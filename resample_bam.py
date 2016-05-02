#!/usr/bin/python

import argparse
import shutil
import sys
import pysam
import numpy as np

def calculate_samples(bamfile, count, ref_prefix=None):
    refcounts = dict()
    for s in pysam.idxstats(bamfile).split("\n"):
        tok = s.rstrip().split("\t")
        if ref_prefix is not None and tok[0].startswith(ref_prefix) == False:
            continue
        refcounts[tok[0]] = int(tok[2])
    coef = (count * 1.0) / sum(refcounts.values())
    return dict((k, (v, int(np.round(v * coef)))) for k, v in refcounts.items())


def get_random_alignments(bfd, obfd, ref, refcounts):
    refsizes = dict(zip(bfd.references, bfd.lengths))
    count = refcounts[ref][1]
    total = refcounts[ref][0]
    #q = np.arange(total)
    #np.random.shuffle(q)
    #q = sorted(q[:count])
    q = sorted(np.random.permutation(np.arange(total))[:count])
    qi = 0
    algni = 0
    if len(q) == 0:
        return
    for algn in bfd.fetch(ref):
        if algni == q[qi]:
            obfd.write(algn)
            qi += 1
            if qi == len(q):
                break
        algni += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-b", "--bamfile", required=True)
    parser.add_argument("-f", "--outfile", required=True)
    parser.add_argument("-c", "--count", required=True, type=int)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    args = parser.parse_args()
    refcounts = calculate_samples(args.bamfile, args.count, "chr")
    bfd = pysam.Samfile(args.bamfile, "rb")
    if bfd.mapped < args.count:
        print "number of reads specified is greater than number of reads in source file"
        exit(1)
    obfd = pysam.Samfile(args.outfile, "wb", template=bfd)
    for ref in sorted(refcounts):
        if args.verbose:
            print "processing %s %d" % (ref, refcounts[ref][1])
            sys.stdout.flush()
        get_random_alignments(bfd, obfd, ref, refcounts)
    obfd.close()
    bfd.close()
    sortedprefix = "%s.sorted" % ".".join(args.outfile.split(".")[:-1])
    if args.verbose:
        print "sorting %s -> %s.bam" % (args.outfile, sortedprefix)
    pysam.sort("-o", "%s.bam" % sortedprefix, args.outfile)
    shutil.move("%s.bam" % sortedprefix, args.outfile)
    if args.verbose:
        print "indexing %s" % args.outfile
    pysam.index(args.outfile)
