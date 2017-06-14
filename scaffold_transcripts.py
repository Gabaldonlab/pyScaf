#!/usr/bin/env python
desc="""Perform scaffolding of de novo transcripts fragments using reference transcripts/peptides
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 13/06/2017
"""

import os, sys, gzip, commands, resource, subprocess
import numpy as np
from datetime import datetime
from FastaIndex import FastaIndex
from pyScaf import percentile, median, mean, pstdev, _check_dependencies

def _lastal(fasta, ref, norearrangements=0, protein=0):
    """Start LAST in local mode and with FastQ input (-Q 1)."""
    # build db
    dbcmd = "lastdb %s %s"
    if protein:
        dbcmd = "lastdb -p %s %s"
    if not os.path.isfile(ref+".suf"):
        os.system(cmd%(ref, ref))
    # run last -s2 both strands for protein
    args1 = ["lastal", "-s2", "-P", str(threads), ref, fasta] # "-fTAB", 
    args2 = ["last-split"]
    args3 = ["maf-convert", "tab"]
    if norearrangements:
        args1.append("-T1")
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc1.stdout)
    proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stderr=sys.stderr, stdin=proc2.stdout)
    return proc3

def _is_rearranged(algs, strand):
    """Return True if rearrangement between query and target"""
    # return True if inversion
    if strand not in (0, 1):
        return True
    # return True if translocation
    for i in range(1, len(algs)):
        if algs[i][3] < algs[i-1][3]:
            return True

def get_hits(fasta, ref, identity, overlap, norearrangements, protein):
    """Resolve & report scaffolds"""
    ## consider splitting into two functions
    ## to facilitate more input formats
    t2hits = {}
    t2size = {}
    q2hits = {}
    for l in gzip.open(fasta+".tab.gz"): #_lastalal(fasta, ref, norearrangements, protein):
        if l.startswith('#'):
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        # filter by identity and overlap
        _identity = 1.0 * score / qalg
        if _identity < identity:# or _overlap < overlap:
            continue
            
        # store t2size
        if t not in t2size:
            t2size[t] = tsize

        # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
        strand = 0 # forward
        if qstrand == "-":
            # reverse
            strand = 1 
            #qstart = qsize - qstart - qalg
        qend, tend = qstart + qalg, tstart + talg
            
        # store q2hits
        alg = (tstart, tend, q, qstart, qend, qsize, strand, score)
        if q not in q2hits:
            q2hits[q] = {t: [score, [alg]]}
        elif t not in q2hits[q]:
            q2hits[q][t] = [score, [alg]]
        else:
            q2hits[q][t][0] += score
            q2hits[q][t][1].append(alg)
            
    # populate t2hits only with best q2t hits
    for q in q2hits:
        # get t with highest score
        t, (score, algs) = sorted(q2hits[q].iteritems(), key=lambda x: x[1][0], reverse=1)[0]
        # check q overlap
        qsize = algs[0][5]
        _overlap  = 1.0 * sum(x[4]-x[3] for x in algs) / qsize
        if _overlap < overlap:
            continue
        # get strand
        strand = mean([x[6] for x in algs])
        # check if no rearrangements
        algs.sort()
        if _is_rearranged(algs, strand):
            #print "rearrangement in:", strand, t, algs
            continue
        # collapse q2t partial matches
        first, last = algs[0], algs[-1]
        _identity = round(1.*sum(x[-1] for x in algs) / sum(x[4]-x[3] for x in algs), 4)
        alg = (first[0], last[1], first[2], first[3], last[4], first[5], strand, _identity)
        # store target matches
        if t not in t2hits:
            t2hits[t] = []
        t2hits[t].append(alg)

    # process 
    return t2hits, t2size
            
def split_many2one(algs, tsize, maxOverlap=.66):
    """Return one2one algs from many to one"""
    # convert alg intervals into array of 0/1 False/True
    a = np.zeros((len(algs), tsize), dtype='bool')
    for i in range(len(algs)):
        a[i][algs[i][0]:algs[i][1]] = True
    # check overlaps between algs and store each overlap into seperate block
    uniq = []
    added = set()
    for i in range(0, len(algs)-1):
        # tstart, tend, q, qstart, qend, qsize, strand, identity
        if algs[i][2] in added:
            continue
        uniq.append([algs[i]])
        for ii in range(i+1, len(algs)):
            if np.all(a[(i, ii),], axis=0).sum() > maxOverlap*a[(i, ii),].sum(axis=1).min():
                #print "", np.all(a[(i, ii),], axis=0).sum(), a[(i, ii),].sum(axis=1).min(), algs[i], algs[ii]
                # update
                uniq[-1].append(algs[ii])
                added.add(algs[ii][2])
    #print uniq
    # [[(1056, 1818, 'TRINITY_DN29363_c1_g1_i2', 267, 1026, 1405, 1.0, 0.5549),
    #   (1081, 1304, 'TRINITY_DN29363_c0_g1_i1', 5, 228, 246, 0.0, 0.4529),
    #   (1467, 1818, 'TRINITY_DN29363_c1_g1_i1', 18, 369, 748, 1.0, 0.6752)]]            
    _algs = []
    if len(uniq)<2:
        return _algs
    '''
    [[(258, 867, 'TRINITY_DN47477_c0_g1_i1', 16, 622, 622, 0.0, 0.7459)],
     [(853, 1328, 'TRINITY_DN7505_c0_g1_i1', 0, 475, 476, 1.0, 0.7684)],
     [(1298, 1886, 'TRINITY_DN79803_c0_g1_i1', 0, 591, 591, 0.0, 0.6937)],
     [(1898, 2226, 'TRINITY_DN71663_c0_g1_i1', 19, 347, 347, 1.0, 0.811)],
     [(2954, 3268, 'TRINITY_DN39262_c0_g1_i1', 21, 335, 401, 1.0, 0.5096)]]'''
    for u in uniq:
        while len(u)>len(_algs):
            _algs.append([])
        for i, _u in enumerate(sorted(u, key=lambda x: x[-1], reverse=1)):
            _algs[i].append(_u)

    # remove contigs covered completely by bigger conting overhang
    # 12289 2494 ENSDART00000057044.7 2687
    #  [(71, 265, 'TRINITY_DN1115_c0_g2_i1', 390, 584, 584, 0.0, 1.0),
    #   (300, 2138, 'TRINITY_DN8977_c0_g1_i2', 304, 2139, 2246, 0.0, 0.9619)]        
    return _algs
            
def scaffold_transcripts(fasta, ref, out, threads, \
                         identity, overlap, maxgap, norearrangements=0, log=sys.stderr, protein=0):
    """Scaffold de novo transcripts using reference transcripts/peptides"""
    # get hits
    t2hits, t2size = get_hits(fasta, ref, identity, overlap, norearrangements, protein)

    i = k = 0
    for i, (t, algs) in enumerate(t2hits.iteritems(), 1):
        # recognise many-to-one alignements
        algs = split_many2one(sorted(algs), t2size[t])
        for _algs in algs:
            # skip singletons
            if len(_algs)<2:
                continue
            k += 1
            print i, k, t, t2size[t], _algs

def main():
    import argparse
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    #parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.12a')
    parser.add_argument("-f", "--fasta", required=1, help="assembly FASTA file")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType('w'), help="output stream [stdout]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max no. of threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="output log to [stderr]")
    
    refo = parser.add_argument_group('Reference-based scaffolding options')
    refo.add_argument("-r", "--ref", "--reference", required=1, help="reference transcripts FastA file")
    refo.add_argument("--identity", default=0.33, type=float, help="min. identity [%(default)s]")
    refo.add_argument("--overlap", default=0.33, type=float, help="min. overlap  [%(default)s]")
    refo.add_argument("-g", "--maxgap", default=100, type=int, help="max. distance between adjacent contigs []")
    refo.add_argument("--norearrangements", action='store_true', help="high identity mode (rearrangements not allowed)")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    #if o.verbose:
    o.log.write("Options: %s\n"%str(o))

    # check if all executables exists & in correct versions
    dependencies = {'lastal': 700, 'lastdb': 700, }
    _check_dependencies(dependencies)
    
    scaffold_transcripts(o.fasta, o.ref, o.out, o.threads, \
                         o.identity, o.overlap, o.maxgap, o.norearrangements, o.log)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    