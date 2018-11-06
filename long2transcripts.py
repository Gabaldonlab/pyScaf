#!/usr/bin/env python
desc="""Get transcripts from long reads
- all-vs-all alignment
- clustering and isoform/paralogs separation
- correction with long & short-reads
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 13/06/2017
"""

import os, sys, gzip, resource, subprocess, pysam
import numpy as np
from datetime import datetime
#from FastaIndex import FastaIndex
#from pyScaf import percentile, median, mean, pstdev, _check_dependencies

class MyGraph(object):
    """Undirected Graph class."""
    def __init__(self, vertices=[]):
        """Construct a graph with the given vertices, but no lines
        """
        self._neighbours = {}
        self._features = {}
        for vertex in vertices:
            self._neighbours[vertex] = set()
            self._features[vertex] = {}

    def add_vertice(self, vertex, key="", value=""):
        """Add vertice to vertices.
        If already present, just update vertice features."""
        if vertex not in self._neighbours:
            self._neighbours[vertex] = set()
            self._features[vertex] = {}
        if key:
            self._features[vertex][key] = value

    def add_vertices(self, vertices, key="", value=""):
        """Add vertice to vertices.
        If already present, just update vertice features."""
        for vertex in vertices:
            if vertex not in self._neighbours:
                self._neighbours[vertex] = set()
                self._features[vertex] = {}
            if key:
                self._features[vertex][key] = value

    def add_line(self, v1, v2):
        """Add a line from v1 to v2 (with no error checking!)
        """
        self._neighbours[v1].add(v2)
        self._neighbours[v2].add(v1)

    def add_lines(self, v1, vertices):
        """Add a lines from v1 to vertices
        """
        for v2 in vertices:
            self._neighbours[v1].add(v2)
            self._neighbours[v2].add(v1)
                        
    def __str__(self):
        """Produce string representation of the graph
        """
        out = 'Vertices: %s \nLines:\n' % self._neighbours.keys() 
        for vertex, neighbours in self._neighbours.items():
            for neighbour in neighbours:
                if vertex < neighbour:
                    out += '\t%s - %s\n' % (vertex,neighbour)
        return out

    def get_clusters(self, key="", value=""):
        """Return connected nodes having given feature (key: value)."""
        processed = set()
        clusters = []
        for v in self._neighbours:
            if v not in processed:
                clusters.append([v])
                processed.add(v)
            i = 0
            while i<len(clusters[-1]): #for v2 in self._neighbours[v]:
                v = clusters[-1][i]
                for v2 in self._neighbours[v]:
                    if v2 in processed:
                        continue
                    clusters[-1].append(v2)
                    processed.add(v2)
                i += 1
        return clusters

def _minimap2(fasta, threads=3, k=10, log=sys.stderr):
    """Start minimap2"""
    # run minimap -m100 -g10000 -r2000 --max-chain-skip 25
    args1 = ["minimap2", "-Hk%s"%k, "-t%s"%threads, "-Xw5", fasta, fasta]
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=log)
    return proc1.stdout

def get_hits_paf(fasta, identity=0.1, overlap=0.3, threads=4, log=sys.stderr):
    """Return hits from PAF file"""
    t2size, q2hits = {}, {}
    if os.path.isfile(fasta+".paf.gz"):
        if log: log.write(" parsing from file...\n")
        handle = gzip.open(fasta+".paf.gz")
    else:
        if log: log.write(" running minimap2...\n")
        handle = _minimap2(fasta, threads, log=log)
    for l in handle:  
        # unpack
        (q, qsize, qstart, qend, qstrand, t, tsize, tstart, tend, matches, algLengths, mapq) = l.split()[:12]
        (matches, qstart, qend, qsize, tstart, tend, tsize) = map(int, (matches, qstart, qend, qsize, tstart, tend, tsize))
        # target always larger
        if qsize>tsize or qsize==tsize and q>t:
            q, qsize, qstart, qend, t, tsize, tstart, tend = t, tsize, tstart, tend, q, qsize, qstart, qend
        # filter by identity and overlap
        qalg = qend - qstart
        _identity = 1.0 * matches / qalg 
        _overlap = 1.0 * qalg / qsize
        if _identity < identity or _overlap < overlap:
            continue

        # store q2hits
        if q not in q2hits:
            q2hits[q] = set([t])
        else:
            q2hits[q].add(t)
    return q2hits

def _lastal(fasta, ref, threads=4):
    """Start LAST in local global and with FastQ input (-Q 1)."""
    # build db
    dbcmd = "lastdb %s %s"
    if not os.path.isfile(ref+".suf"):
        os.system(dbcmd%(ref, ref))
    # run last -s2 both strands for protein -C1? -s1?
    args1 = ["lastal", "-a1", "-b1", "-T1", "-f TAB", "-P", str(threads), ref, fasta] # "-fTAB", 
    args2 = ["awk '$6>$11 || $6==$11 && $2>$7'"]
    args3 = ["gzip"]
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc1.stdout)
    proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc2.stdout)
    return proc3.stdout

def get_hits(fasta, identity=0.3, overlap=0.8, threads=4):
    """Return hits from TAB file."""
    t2size = {}
    q2hits = {}
    for l in gzip.open(fasta+".C2.global.tab.gz"): #_lastal(fasta, ref, threads): # 
        if l.startswith('#'):
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        # filter smaller
        if qsize>tsize or qsize==tsize and q>t:
            continue
        # filter by identity and overlap - score is +1/-1 for match/mismatch&indel
        mismatches = (qalg - score) / 2.0
        _identity = 1 - 1.0 * mismatches / qalg # 1.0 * (qalg + score) / qalg / 2
        _overlap = 1.0 * qalg / qsize
        if _identity < identity or _overlap < overlap:
            continue

        # store t2size
        if t not in t2size:
            t2size[t] = tsize
        if q not in t2size:
            t2size[q] = qsize
   
        # store q2hits
        algdata = (score, qalg, _identity, _overlap)
        if q not in q2hits:
            q2hits[q] = [(t, algdata)]
        else:
            q2hits[q].append((t, algdata))

    print(len(q2hits), len(t2size))
    return q2hits, t2size

def get_cigar2bases(cigartuples):
    """Return cigar2base dictionary"""
    cigar2bases = {}
    for c, b in cigartuples:
        if c in cigar2bases:
            cigar2bases[c].append(b)
        else:
            cigar2bases[c]=[b]
    return cigar2bases

def _get_dv(tags):
    dv = 1.0
    for t, v in tags:
        if t=='dv': return v
    return dv
    
def get_hits_minimap2_sam(fasta, identity=0.7, overlap=0.95, threads=4, maxindel=20):
    """Return hits from SAM file. Works with minimap2 sam files!! UNFINISHED!"""
    q2hits = {}
    sam = pysam.AlignmentFile(fasta+".k11.sam.gz")
    t2size = {r: l for r, l in zip(sam.references, sam.lengths)}
    for r in sam: 
        q, t = r.qname, sam.references[r.rname]
        # target always larger
        qstart, qend = r.qstart, r.qend
        if t2size[q]>t2size[t] or t2size[q]==t2size[t] and q>t:
            q, qstart, qend, t, = t, r.reference_start, r.reference_end, q
            
        # unload cigar
        c2b = get_cigar2bases(r.cigar)
        # skip if too much clipped
        if 4 in c2b and max(c2b[4])>100:
            continue
        # skip if indel>maxindel
        if 1 in c2b and max(c2b[1])>maxindel or 2 in c2b and max(c2b[2])>maxindel:
            continue
        qalg = qend - qstart
        score = 0
        # filter by identity and overlap
        _identity = 1.0 - _get_dv(r.tags)
        _overlap = 1.0 * qalg / t2size[q]
        if _identity < identity or _overlap < overlap:
            continue

        # store q2hits
        algdata = (score, qalg, _identity, _overlap)
        if q not in q2hits:
            q2hits[q] = [(t, algdata)]
        else:
            q2hits[q].append((t, algdata))

    return q2hits, t2size
    
def get_hits_last_sam(fasta, identity=0.7, overlap=0.95, threads=4, maxindel=20):
    """Return hits from SAM file. Works with LAST sam files!!"""
    q2hits = {}
    sam = pysam.AlignmentFile(fasta+".global.sam.gz")
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
    return q2hits, t2size
        
def _get_best_cluster(hits, clusters, q2cluster, t2size):
    """Return best target among clusters"""
    # alg to the longest isoform first
    for t, algdata in sorted(hits, key=lambda x: t2size[x[0]], reverse=1):
        if t in clusters:
            return t
        elif q2cluster and t in q2cluster:
            return q2cluster[t]

def get_clusters(outdir, fasta, q2hits, t2size, join_by_other_members=False):
    """Return clusters or sequences based on all-vs-all matches.
    Clusters are stored in outdir/clusterXXXXXX.fa files. 

    If join_by_other_members is True, sequences are added to cluster
    if another member of the cluster is matched (more relaxed matches).
    If it's false, then valid match to the longest sequence in the
    cluster is required. 
    """
    clusters = {}
    q2cluster = {}
    for q, tsize in sorted(t2size.iteritems(), key=lambda x: x[1], reverse=1):
        # create new cluster
        if q not in q2hits:
            clusters[q] = [q]
            q2cluster[q] = q
            continue
        # add to existing cluster
        t = _get_best_cluster(q2hits[q], clusters, q2cluster, t2size) 
        if t:    
            clusters[t].append(q)
            if join_by_other_members:
                q2cluster[q] = t
        else:
            clusters[q] = [q]
            if join_by_other_members:
                q2cluster[q] = q
    print(len(clusters), len(q2cluster))
    faidx = pysam.FastaFile(fasta)
    for i, (t, cl) in enumerate(sorted(clusters.iteritems(), key=lambda x: t2size[x[0]], reverse=1), 1):
        print i, t, t2size[t], len(cl)
        fn = "cluster%6s.fa"%i
        fn = fn.replace(" ", "0")
        with open(os.path.join(outdir, fn), "w") as out:
            for q in cl:
                out.write(">%s\n%s\n"%(q, faidx[q]))
    return clusters

def reads2clusters(outdir, fasta, t2size, threads=3, log=sys.stderr):
    """Align reads all-vs-all and get clusters likely representing
    all transcripts for given gene or group of paralogs.
    
    """
    if log: log.write("Generating clusters...\n")
    q2hits = get_hits_paf(fasta, threads=threads, log=log)
    # get clusters
    ## init graph with all reads
    g = MyGraph(t2size.keys())
    # add matches
    for q, hits in q2hits.iteritems():
        g.add_lines(q, hits)
    # get clusters
    clusters = g.get_clusters()
    # & report    
    faidx = pysam.FastaFile(fasta)
    k = 0
    for i, cluster in enumerate(clusters, 1): 
        k += len(cluster)
        #print(i, max(t2size[e] for e in cluster), len(cluster))
        fn = "cluster%6s.fa"%i
        fn = fn.replace(" ", "0")
        with open(os.path.join(outdir, fn), "w") as out:
            for q in cluster:
                out.write(">%s\n%s\n"%(q, faidx[q]))
    log.write(" %s seqs in %s cluster(s)\n"%(k, i))
    return clusters    
    
def scaffold_transcripts(fasta, ref, out, threads, \
                         identity, overlap, maxgap, norearrangements=0, log=sys.stderr, protein=0):
    """Scaffold de novo transcripts using reference transcripts/peptides"""
    # get hits
    t2hits, t2size = get_hits(fasta, ref, identity, overlap, threads, norearrangements, protein)

    # generate scaffolds
    i = k = 0
    faidx = FastaIndex(fasta)
    added = {c: 0 for c in faidx}
    for i, (t, algs) in enumerate(t2hits.iteritems(), 1):
        # recognise many-to-one alignements
        algs = split_many2one(sorted(algs), t2size[t])
        for _algs in algs:
            # skip singletons and overlap below cut-off
            if len(_algs)<2:# or sum(a[1]-a[0] for a in _algs)<overlap*t2size[t]:
                continue
            k += 1
            sys.stderr.write("%s / %s  %s   \r"%(i, k, t))#, t2size[t], _algs
            name = 'scaffold{0:07d}'.format(k) + " %s:%s-%s %s"%(t, _algs[0][0], _algs[-1][1], t2size[t])
            added = report_scaffold(faidx, added, out, name, _algs, protein)

    # store unscaffolded contigs
    for c in faidx:
        if not added[c]:
            out.write(faidx[c])

    sys.stderr.write("Joined %s out of %s contigs into %s scaffolds\n"%(len([c for c in added if added[c]]), len(faidx), k))

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
            out = "".join(p.stdout.readlines())
            curver = out.split()[-1]
            if not curver.isdigit():
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            elif int(curver)<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/pyScaf#dependencies !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)
    
def main():
    import argparse
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    #parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.10a')
    parser.add_argument("-f", "--fasta", required=1, help="reads FASTA file")
    parser.add_argument("-o", "--outdir", default='long2transcripts', help="output dir [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max no. of threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="output log to [stderr]")
    parser.add_argument("-i", "--maxindel", default=20, type=int, help="max allowed indel [%(default)s]")
    parser.add_argument("--test", action='store_true', help="test run")
    
    """ # kmer maxclip
    refo = parser.add_argument_group('Reference-based scaffolding options')
    refo.add_argument("-r", "--ref", "--reference", required=1, help="reference transcripts FastA file")
    refo.add_argument("--identity", default=0.33, type=float, help="min. identity [%(default)s]")
    refo.add_argument("--overlap", default=0.33, type=float, help="min. overlap  [%(default)s]")
    refo.add_argument("-g", "--maxgap", default=100, type=int, help="max. distance between adjacent contigs []")
    refo.add_argument("--norearrangements", action='store_true', help="high identity mode (rearrangements not allowed)")
    """
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    #if o.verbose:
    o.log.write("Options: %s\n"%str(o))

    
    # check if all executables exists & in correct versions
    dependencies = {'lastal': 959, 'lastdb': 959, 'minimap2': 0}
    _check_dependencies(dependencies)

    if not os.path.isdir(o.outdir):
        os.makedirs(o.outdir)

    if o.test:
        """
        q2hits, q2size = get_hits(fasta)
        q2hits, q2size = get_hits_paf(fasta)
        #"""
        q2hits, q2size = get_hits_last_sam(fasta)
        #q2hits, q2size = get_hits_minimap2_sam(fasta)
        clusters = get_clusters(o.outdir, fasta, q2hits, q2size)
        return

    # load fasta index
    faidx = pysam.FastaFile(o.fasta)
    t2size = {r: l for r, l in zip(faidx.references, faidx.lengths)}    

    # reads2clusters
    clusters = reads2clusters(o.outdir, o.fasta, t2size, o.threads, o.log)
        
    # align all clusters precisely in parallel
    #scaffold_transcripts(o.fasta, o.ref, o.out, o.threads, o.identity, o.overlap, o.maxgap, o.norearrangements, o.log)
            
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    