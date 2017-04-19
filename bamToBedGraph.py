#!/usr/bin/python

"""Generate a bedGraph file from a given indexed BAM file

Constraints:  BAM file should be sorted by genomic coordinate and indexed

For a full list of options, run the script with the '--help' option
"""

__version__ = 2.0

import sys
import argparse
import re

import BedGraphTable as bgt
import hamrhein.genome.bam as bam

class BamfilterAction(argparse.Action):
    bamfilter = 0

    def __init__(self, **kwargs):
        super(BamfilterAction, self).__init__(**kwargs)

        self.dest = "bamfilter"
        self.nargs = 0
        self.default = bam.BAMQ_ALL

    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == "--all":
            BamfilterAction.bamfilter |= bam.BAMQ_ALL
            BamfilterAction.bamfilter ^= bam.BAMQ_UNMAPPED

        if option_string == "--unique":
            BamfilterAction.bamfilter |= bam.BAMQ_UNIQUES

        if option_string == "--multi":
            BamfilterAction.bamfilter |= bam.BAMQ_MULTIS

        if option_string == "--match":
            BamfilterAction.bamfilter |= bam.BAMQ_MATCHES

        if option_string == "--splice":
            BamfilterAction.bamfilter |= bam.BAMQ_SPLICES

        BamfilterAction.bamfilter |= bam.BAMQ_MAPPED
        setattr(namespace, self.dest, BamfilterAction.bamfilter)


def insertAlignment(tbl, generator, strand):
    count = 0

    for r in generator:
        if strand == "+" and r.is_reverse:
            continue

        if strand == "-" and not r.is_reverse:
            continue

        if strand is None:
            tbl.insertPysamAlignment(r, False, False)

        else:
            tbl.insertPysamAlignment(r, False, True)

        count += 1

    return count


def main(infile, outfile):
    bamfilt = bam.BamFilter(args.bamfilter)

    if args.normalize:
        norm_factor = 1000000.0 / len(infile)
    else:
        norm_factor = 1.0

    if args.region:
        reg_split = re.split("[:-]", args.region)
        ref = reg_split[0]

        if len(reg_split) == 3:
            length = int(reg_split[2]) - int(reg_split[1])

        else:
            ridx = infile.bam.references.index(ref)
            length = infile.bam.lengths[ridx]

        bgtable = bgt.BedGraphTable(ref, length, norm_factor)
        alignments = infile.fetch(bamfilt, args.region)

        count = insertAlignment(bgtable, alignments, args.strand)

        for line in bgtable:
            outfile.write("%s\n" % line)

        if args.verbose:
            print("\t%d alignments" % count)

    else:
        for ref, length in infile.references:
            if args.chromonly and not ref.startswith("chr"):
                continue

            if args.verbose:
                print("Processing %s" % ref)

            bgtable = bgt.BedGraphTable(ref, length, norm_factor)
            alignments = infile.fetch(bamfilt, reference=ref, start=1, end=length)

            count = insertAlignment(bgtable, alignments, args.strand)

            for line in bgtable:
                outfile.write("%s\n" % line)

            if args.verbose:
                print("\t%d alignments" % count)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a bedGraph file from a given indexed BAM file")
    parser.add_argument("--all", "--unique", "--multi", "--match", "--splice",
                        action=BamfilterAction, type=int,
                        help="set various query constraints", default="--all")

    parser.add_argument("--chromonly", action="store_true",
                        help="Only process references beginning with ``chr''")

    parser.add_argument("--normalize", action="store_true",
                        help="Normalize results", default=False)

    parser.add_argument("--region",
                        help="Report only reads within REGION", default=None)

    parser.add_argument("--strand", metavar="STRAND",
                        help="Report only reads from STRAND", default=None)

    parser.add_argument("--verbose", action="store_true",
                        help="Show progress", default=False)

    parser.add_argument("--version", action="store_true",
                        help="Print version and exit", default=False)

    parser.add_argument("infile")
    parser.add_argument("outfile")

    if len(sys.argv) < 3:
        parser.print_help()
        print("")
        sys.exit(-1)

    args = parser.parse_args()

    if args.version:
        print("bamToBed.py version %s" % __version__)
        sys.exit(0)

    try:
        ifs = bam.Bam(args.infile)
    except IOError, value:
        print("BAM file error: %s" % value)
        sys.exit(-1)

    try:
        ofs = open(args.outfile, "w")
    except IOError, value:
        print("bedGraph file error: %s" % value)
        sys.exit(-1)

    main(ifs, ofs)

    ofs.close()
    sys.exit(0)
