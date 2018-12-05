#!/usr/bin/env python
desc="""Get transcripts from long reads
- all-vs-all alignment
- clustering and isoform/paralogs separation
- correction with long & short-reads


TBD: split reads into 1M batches, process with minimap2 + mcl
and then process resulting transcripts with minimap2 (or lastal) + mcl.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 2/11/2018
"""

import glob, os, sys, gzip, resource, subprocess, pysam
import numpy as np
from datetime import datetime
from multiprocessing import Pool

def _minimap2_mcl(fasta, threads=1):
    """Align globally with minimap2 and cluster with mcl."""
    opts = "-Hk13 -X".split() 
    args1 = ["minimap2", "-t %s" %threads, fasta, fasta] + opts # ["-xava-ont", ]
    args2 = ["cut", "-f1,6,10"]
    args3 = ["mcl", "-", "-o", "-", "-I", "2.0", "--abc", "-te", "%s"%threads]
    #print(" | ".join(" ".join(a) for a in (args1, args2, args3)))
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=open(fasta+".minimap2.log", "w"))
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stdin=proc1.stdout, stderr=sys.stderr)
    proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stdin=proc2.stdout, stderr=open(fasta+".mcl.log", "w"))
    return proc3.stdout

def _lastal2sam(fasta, threads=4):
    """Align globally with last and report SAM."""
    # build db
    dbfn = "%s.lastdb"%fasta
    os.system("lastdb -P%s %s %s"%(threads, dbfn, fasta)) #if not os.path.isfile(fasta+".suf"):
    # run last -s2 both strands for protein -C1? -s1? "-f TAB"; slowing down ~6x:  "-a1", "-b1"
    args1 = ["lastal", "-N50", "-a1", "-b1", "-s1", "-T1", "-P%s"%threads, dbfn, fasta] # "-fTAB", FastQ input (-Q 1)
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
    q2hits = {}
    sam = pysam.AlignmentFile(_lastal2sam(fasta, threads)) # has to be processed by samtools
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

def get_clusters(fasta, identity, overlap, threads, maxindel):
    """Save longest sequence for each subcluster based on all-vs-all matches.
    Return number of clusters. 
    """
    faidx = pysam.FastaFile(fasta)
    q2hits, t2size = get_hits_last_sam(fasta, identity, overlap, threads, maxindel)
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
    # save longest isoform for each subcluster
    with open("%s.isoforms.fa"%fasta, "w") as out:
        for i, (t, cl) in enumerate(sorted(clusters.iteritems(), key=lambda x: t2size[x[0]], reverse=1), 1):
            out.write(">%s\n%s\n"%(t, faidx[t]))
    return len(clusters)
            
def worker2(args):
    clusterfn, identity, overlap, threads, maxindel = args
    return get_clusters(clusterfn, identity, overlap, threads, maxindel)
    
def batch2clusters(fasta, identity, overlap, threads, maxindel, log=sys.stderr, minCluster=10):
    """Compute all-vs-all matches with minimap2 and cluster with mcl"""
    faidx = pysam.FastaFile(fasta)
    singletons = set(faidx.references)
    # get clusters
    i = k = isoforms = 0
    clusterfnames = []
    outclusters = open("%s.clusters.txt"%fasta, "w")
    for l in _minimap2_mcl(fasta, threads):
        cluster = l[:-1].split()
        if len(cluster)<minCluster:
            continue
        i += 1
        k += len(cluster)
        singletons.difference_update(cluster)
        fn = "cluster%6s.fa"%i
        fn = fn.replace(" ", "0")
        outclusters.write("%s\t%s\t%s\n"%(fn, len(cluster), " ".join(cluster)))
        clusterfn = "%s.%s"%(fasta, fn) 
        with open(clusterfn, "w") as out:
            for q in cluster:
                out.write(">%s\n%s\n"%(q, faidx[q]))
        clusterfnames.append(clusterfn)
        # get longest
        #isoforms += get_clusters(clusterfn, identity, overlap, threads, maxindel)
    outclusters.close()
        
    # get longest
    p = Pool(threads)
    for _i in p.imap_unordered(worker2, [(fn, identity, overlap, 1, maxindel) for fn in clusterfnames]):
        isoforms += _i
    
    with open("%s.singletons.fa"%fasta, "w") as out:
        for q in singletons:
            out.write(">%s\n%s\n"%(q, faidx[q]))
    log.write(" %s isoforms from %s seqs in %s clusters (with >=%s members) + %s singletons\n"%(isoforms+len(singletons), k, i, minCluster, len(singletons)))
    cmd = "cat %s*.isoforms.fa %s.singletons.fa > %s.transcripts.fa"%(fasta, fasta, fasta)
    os.system(cmd)
    return "%s.transcripts.fa"%fasta

def fasta_parser(handle):
    """Simple fasta parser"""
    data = []
    for l in handle:
        if l.startswith(">"):
            if data:
                yield data
            data = [l, ]
        else:
            data.append(l)
    if data:
        yield data
    
def fastq_parser(handle):
    """Simple fastq parser assuming 4 lines per entry returning header and sequence"""
    data = []
    for l in handle:
        data.append(l)
        if len(data)==4:
            yield data[:2]
            data = []
    if len(data)==4:
        yield data[:2]

def worker(args):
    fn, identity, overlap, threads, maxindel = args
    batch2clusters(fn, identity, overlap, threads, maxindel)
        
def process_in_batches(outdir, fasta, identity, overlap, threads, maxindel, maxseq=1000000, log=sys.stderr):
    """Process reads in batches"""
    handle = open(fasta)
    if fasta.endswith('.gz'):
        handle = gzip.open(fasta)
        fasta = fasta[:-3] # strip .gz
    parser = fasta_parser(handle)
    if fasta.endswith(('.fq', '.fastq')):
        parser = fastq_parser(handle) #

    fastas = []
    for i, r in enumerate(parser):
        if not i%maxseq:
            subdir = "batches/%4s"%(len(fastas)+1, )
            subdir = subdir.replace(' ', '0')
            fn = os.path.join(outdir, subdir, "reads.fa")
            os.makedirs(os.path.dirname(fn))
            fastas.append(fn)
            out = open(fn, "w")
        # write
        out.write(">"+"".join(r)[1:]) # fastq @ fasta >, thus skip first element
    out.close()

    # process all reads
    for fn in fastas:
        batch2clusters(fn, identity, overlap, threads, maxindel)
    '''
    if len(fastas)<2:
        # use all threads locally if only 1 batch
        for fn in fastas:
            batch2clusters(fn, identity, overlap, threads, maxindel)
    else:
        p = Pool(threads)
        for fn in p.imap_unordered(worker, [(fn, identity, overlap, 1, maxindel) for fn in fastas]):
            pass
    '''
    outfn = "%s/batches.transcripts.fa"%outdir
    os.system("cat %s/batches/*/reads.fa.transcripts.fa > %s"%(outdir, outfn))
    return outfn
            
def correct_racon(fasta, reads, sam="", minq=7, threads=4, log=0):
    """Correct longest read with racon"""
    racon = "%s.racon.fa"%fasta    
    if os.path.isfile(racon):
        return racon
    if not sam or not os.path.isfile(sam):
        sam = "%s.paf.gz"%fasta
        cmd0 = "minimap2 -xava-ont -t%s %s %s 2> %s.minimap2.log | gzip > %s"%(threads, fasta, reads, fasta, sam)
        if log: log.write('  %s\n'%cmd0)
        os.system(cmd0)

    # <sequences> <overlaps> <target sequences>
    cmd = "racon -u -t %s -q %s %s %s %s > %s 2> %s.log"%(threads, minq, reads, sam, fasta, racon, racon)
    if log: log.write('  %s\n'%cmd)
    os.system(cmd)
    return racon
    
def long2transcripts(outdir, fasta, threads, identity, overlap, maxindel, maxseq=1000000, log=sys.stderr):
    """Scaffold de novo transcripts using reference transcripts/peptides"""
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    tfn = os.path.join(outdir, "batches.transcripts.fa.transcripts.fa")
    if not os.path.isfile(tfn):        
        if log: log.write(logger("Generating clusters in batches..."))
        batchesfn = process_in_batches(outdir, fasta, identity, overlap, threads, maxindel, maxseq, log)

        # combine all transcripts from batches
        if log: log.write(logger("Generating clusters from batches..."))
        tfn = batch2clusters(batchesfn, identity, overlap, threads, maxindel, log, minCluster=2)

    # correct racon
    if log: log.write(logger("Correcting with long reads..."))
    tfn = correct_racon(tfn, fasta, threads=threads)
    
    return tfn
    
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
    parser.add_argument("--overlap", default=0.95, type=float, help="min. overlap [%(default)s]")
    parser.add_argument("-m", "--maxseq", default=100000, type=float, help="max no. of seq in batch [%(default)s]")
    #parser.add_argument("--test", action='store_true', help="test run")
    
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
    tfn = long2transcripts(o.outdir, o.long, o.threads, o.identity, o.overlap, o.maxindel, o.maxseq, o.log)

    # correct pilon
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
    