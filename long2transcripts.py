#!/usr/bin/env python
desc="""Get transcripts from long reads
- all-vs-all alignment
- clustering and isoform/paralogs separation
- correction with long & short-reads
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 2/11/2018
"""

import glob, os, sys, gzip, resource, subprocess, pysam
import numpy as np
from datetime import datetime

def get_clusters_minimap2_mcl(outdir, fasta, threads=4, log=sys.stderr):
    """Compute all-vs-all matches with minimap2 and cluster with mcl"""
    # init graph with all reads
    faidx = pysam.FastaFile(fasta)
    # minimap2
    paffn = "%s.paf.gz"%fasta
    if not os.path.isfile(paffn):
        #cmd1 = "minimap2 -Hk10 -Xw5 -t%s %s %s 2> %s.log | gzip > %s"%(threads, fasta, fasta, paffn, paffn)
        cmd1 = "minimap2 -x ava-ont -t%s %s %s 2> %s.log | gzip > %s"%(threads, fasta, fasta, paffn, paffn)
        if log: log.write(" running minimap2...\n  %s\n"%cmd1)
        os.system(cmd1)
    # mcl
    mclfn = "%s.mcl"%paffn
    if not os.path.isfile(mclfn):
        cmd2 = "zcat %s | cut -f1,6,10 | mcl - -I 2.0 --abc -te %s -o %s 2> %s.log"%(paffn, threads, mclfn, mclfn)
        if log: log.write(" running mcl...\n  %s\n"%cmd2)
        os.system(cmd2)
    # get clusters
    k = 0    
    clusters, clusterfnames = [], []
    outclusters = open(os.path.join(outdir, "clusters.txt"), "w")
    for i, l in enumerate(open(mclfn), 1):
        cluster = l[:-1].split()
        clusters.append(cluster)
        k += len(cluster)
        fn = "cluster%6s.fa"%i
        fn = fn.replace(" ", "0")
        outclusters.write("%s\t%s\t%s\n"%(fn, len(cluster), " ".join(cluster)))
        clusterfn = os.path.join(outdir, fn)
        with open(clusterfn, "w") as out:
            for q in cluster:
                out.write(">%s\n%s\n"%(q, faidx[q]))
        clusterfnames.append(clusterfn)
    log.write(" %s seqs in %s cluster(s)\n"%(k, i))
    outclusters.close()
    #print(sorted(map(len, clusters), reverse=1)[:20])
    return clusterfnames 
    
def _lastal2sam(fasta, threads=4):
    """Align globally with last and report SAM."""
    # build db
    dbfn = "%s.lastdb"%fasta
    os.system("lastdb %s %s"%(dbfn, fasta)) #if not os.path.isfile(fasta+".suf"):
    # run last -s2 both strands for protein -C1? -s1? "-f TAB", 
    args1 = ["lastal", "-a1", "-b1", "-T1", "-P", str(threads), dbfn, fasta] # "-fTAB", FastQ input (-Q 1)
    args2 = ["maf-convert", "sam"]
    args3 = ["samtools", "view", "-ShT", fasta]
    #print(" | ".join(" ".join(a) for a in (args1, args2, args3)))
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc1.stdout)
    proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc2.stdout)
    return proc3.stdout

def get_cigar2bases(cigartuples):
    """Return cigar2base dictionary"""
    cigar2bases = {}
    for c, b in cigartuples:
        if c in cigar2bases:
            cigar2bases[c].append(b)
        else:
            cigar2bases[c]=[b]
    return cigar2bases
            
def get_hits_last_sam(fasta, identity=0.7, overlap=0.95, threads=4, maxindel=20, log=''):
    """Return hits from SAM file. Works with LAST sam files!!"""
    # get handle
    if not os.path.isfile(fasta+".sam.gz"):
        #handle = _lastal2sam(fasta, threads)
        cmd = "lastdb %s.lastdb %s; "%(fasta, fasta)
        cmd += "lastal -a1 -b1 -T 1 -P4 %s.lastdb %s | maf-convert sam | samtools view -ShT %s | gzip > %s.sam.gz"%(fasta, fasta, fasta, fasta)
        if log: log.write("  %s\n"%cmd)
        os.system(cmd)
    q2hits = {}
    sam = pysam.AlignmentFile(fasta+".sam.gz") 
    t2size = {r: l for r, l in zip(sam.references, sam.lengths)}
    for r in sam: 
        q, t = r.qname, sam.references[r.rname]
        # target always larger
        if t2size[q]>t2size[t] or t2size[q]==t2size[t] and q>t:
            continue
        # unload cigar
        c2b = get_cigar2bases(r.cigar)
        # skip if too much clipped
        if 5 in c2b and max(c2b[5])>100:
            continue
        # skip if indel>maxindel
        if 1 in c2b and max(c2b[1])>maxindel or 2 in c2b and max(c2b[2])>maxindel:
            continue
        qalg = r.qend - r.qstart
        score = matches = sum(c2b[7])
        # filter by identity and overlap
        _identity = 1.0 * matches / qalg
        _overlap = 1.0 * qalg / t2size[q]
        if _identity < identity or _overlap < overlap:
            continue

        # store q2hits
        algdata = (score, qalg, _identity, _overlap)
        if q not in q2hits:
            q2hits[q] = [(t, algdata)]
        else:
            q2hits[q].append((t, algdata))
    # remove lastdb FASTA.lastdb.*
    for fn in glob.glob("%s.lastdb.*"%fasta):
        os.unlink(fn)
    return q2hits, t2size
        
def _get_best_cluster(hits, clusters, t2size):
    """Return best target among clusters"""
    # alg to the longest isoform first
    for t, algdata in sorted(hits, key=lambda x: t2size[x[0]], reverse=1):
        if t in clusters:
            return t

def get_clusters(fasta, q2hits, t2size):
    """Return list of fnames with longest sequence for clusters of sequences based on all-vs-all matches.
    Clusters are stored in FASTA.subclusterXXXX.fa files. 
    """
    clusters = {}
    for q, tsize in sorted(t2size.iteritems(), key=lambda x: x[1], reverse=1):
        # create new cluster
        if q not in q2hits:
            clusters[q] = [q]
            continue
        # add to existing cluster
        t = _get_best_cluster(q2hits[q], clusters, t2size) 
        if t:    
            clusters[t].append(q)
        else:
            clusters[q] = [q]

    fnames = []
    #print(len(clusters))#, len(q2cluster))
    faidx = pysam.FastaFile(fasta)
    for i, (t, cl) in enumerate(sorted(clusters.iteritems(), key=lambda x: t2size[x[0]], reverse=1), 1):
        #print i, t, t2size[t], len(cl)
        # store longest read are representative
        ifn = ".longest%4s.fa"%i
        ifn = ifn.replace(" ", "0")
        with open(fasta+ifn, "w") as out:
            out.write(">%s\n%s\n"%(t, faidx[t]))
        # don't store subcluster & symlink
        if len(cl)<2:
            # symlink as OUTDIR/FASTA.racon.fa -> FASTA
            os.symlink(os.path.basename(fasta+ifn), os.path.join(os.path.dirname(fasta+ifn), os.path.basename(fasta+ifn)+".racon.fa"))
            continue
        # store all reads
        fn = ifn.replace("longest", "subcluster")
        with open(fasta+fn, "w") as out:
            for q in cl:
                out.write(">%s\n%s\n"%(q, faidx[q]))
        fnames.append((fasta+ifn, fasta+fn))
    return fnames #clusters

def correct_racon(fasta, sam, reads, minq=7, log=0):
    """Correct longest read with racon"""
    # <sequences> <overlaps> <target sequences>
    cmd = "racon -u -q %s %s %s %s > %s.racon.fa 2> %s.racon.fa.log"%(minq, reads, sam, fasta, fasta, fasta)
    if log: log.write('  %s\n'%cmd)
    os.system(cmd)
    
def long2transcripts(outdir, fasta, threads, identity, overlap, maxindel, log=sys.stderr):
    """Scaffold de novo transcripts using reference transcripts/peptides"""
    # load fasta index - skip as it takes lots of memory for large files!
    #faidx = pysam.FastaFile(fasta)
    #t2size = {r: l for r, l in zip(faidx.references, faidx.lengths)}    

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    transcriptsfn = os.path.join(outdir, "transcripts.fa")
    if os.path.isfile(transcriptsfn):
        return transcriptsfn
        
    if log: log.write(logger("Generating clusters..."))# for %s seqs..."%len(faidx)))
    #clusters = get_clusters_minimap2(outdir, fasta, threads=threads, log=log)
    clusterfnames = get_clusters_minimap2_mcl(outdir, fasta, threads, log)

    # separate isoforms using all-vs-all LASTal
    if log: log.write(logger("Separating isoforms from %s clusters..."%len(clusterfnames)))
    for fn in clusterfnames:
        if log: log.write(" %s   \r"%fn)
        faidx = pysam.FastaFile(fn)
        if len(faidx)<2:
            # symlink as OUTDIR/FASTA.racon.fa -> FASTA
            os.symlink(os.path.basename(fn), os.path.join(os.path.dirname(fn), os.path.basename(fn)+".racon.fa"))
            continue
        q2hits, t2size = get_hits_last_sam(fn, identity, overlap, threads, maxindel)
        fnames = get_clusters(fn, q2hits, t2size)
        # correct with long reads
        for longestfn, subclusterfn in fnames:
            correct_racon(longestfn, fn+".sam.gz", subclusterfn)

    # concat all transcripts
    cmd = "cat %s/*.racon.fa > %s"%(outdir, transcriptsfn)
    if log: log.write(" %s\n"%cmd)
    os.system(cmd)
        
    return transcriptsfn
    
def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return "".join(p.stdout.readlines())

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in dependencies.items():
        out = _check_executable(cmd)
        if "not found" in out:
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            p = subprocess.Popen([cmd, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = p.stdout.readlines() #"".join(p.stdout.readlines())
            # ".".join(out.split()[-1].strip('"').split('-')[0].split('_')[0].split('.')[:2]) # openjdk version "1.8.0_181" > 1.8
            curver = out[0].split()[-1].split('-')[0] # samtools 1.8-24-g2035cf3 > 1.8
            try:
                curver = float(curver)
            except:
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            if curver<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/pyScaf#dependencies !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)

def get_memusage():
    return "[%5i Mb] "%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, )            

def logger(mssg):
    timestamp = "[%s]"% datetime.ctime(datetime.now())
    return "%s %s %s\n"%(timestamp, mssg, get_memusage())
    
def main():
    import argparse
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    #parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.10a')
    parser.add_argument("-l", "--long", required=1, help="FASTA/Q file with long reads")
    parser.add_argument("-s", "--short", nargs="*", help="FASTQ files with short reads")
    parser.add_argument("-o", "--outdir", default='long2transcripts', help="output dir [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max no. of threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="output log to [stderr]")
    parser.add_argument("-i", "--maxindel", default=20, type=int, help="max allowed indel [%(default)s]")
    parser.add_argument("--identity", default=0.50, type=float, help="min. identity [%(default)s]")
    parser.add_argument("--overlap", default=0.95, type=float, help="min. overlap  [%(default)s]")
    #parser.add_argument("--test", action='store_true', help="test run")
    
    '''# kmer maxclip
    refo = parser.add_argument_group('Reference-based scaffolding options')
    #refo.add_argument("-r", "--ref", "--reference", required=1, help="reference transcripts FastA file")
    refo.add_argument("--identity", default=0.33, type=float, help="min. identity [%(default)s]")
    refo.add_argument("--overlap", default=0.33, type=float, help="min. overlap  [%(default)s]")
    #refo.add_argument("-g", "--maxgap", default=100, type=int, help="max. distance between adjacent contigs []")
    #refo.add_argument("--norearrangements", action='store_true', help="high identity mode (rearrangements not allowed)")
    '''
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    #if o.verbose:
    o.log.write("Options: %s\n"%str(o))
    log = o.log
    
    # check if all executables exists & in correct versions
    dependencies = {'lastal': 959, 'lastdb': 959, 'minimap2': 0, 'mcl': 0, 'racon': 0, 'samtools': 1.8} #, java: 0, }
    _check_dependencies(dependencies)

    # get transcripts fname (and correct with long reads)
    tfn = long2transcripts(o.outdir, o.long, o.threads, o.identity, o.overlap, o.maxindel, o.log)

    if o.short and not os.path.isfile("%s.last.bam.fasta"%tfn):
        if not os.path.isfile("%s.last.bam"%tfn):
            if log: log.write(logger("Correcting transcripts with short reads..."))
            cmd1 = "lastdb %s %s; lastal -P %s -a6 -b6 -T1 -Q1 %s %s | maf-convert sam | "%(tfn, tfn, o.threads, tfn, " ".join(o.short)) + \
                    "samtools view -SbuT %s | samtools sort -o %s.last.bam; samtools index %s.last.bam"%(tfn, tfn, tfn)
            if log: log.write(" %s\n"%cmd1)
            os.system(cmd1)
        cmd2 = "java -jar ~/src/pilon-1.22.jar --threads %s --genome %s --bam %s.last.bam --output %s.last.bam > %s.last.bam.fasta.log"%(o.threads, tfn, tfn, tfn, tfn)
        if log: log.write(" %s\n"%cmd2)
        os.system(cmd2)
    #
    if log: log.write(logger("Done!"))
    
            
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    