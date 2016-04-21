##################################
#                                #
# Last modified 04/26/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
import string

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAM chrom.sizes outfilename' % sys.argv[0]
        sys.exit(1)

    BAM = sys.argv[1]
    chromsizes = sys.argv[2]
    outputfilename = sys.argv[3]

    ChrDict={}
    linelist = open(chromsizes)
    for line in linelist:
        fields = line.strip().split('\t')
        chr = fields[0]
        size = fields[1]
        ChrDict[chr]=size

    outfile = open(outputfilename, "w")
    outline = '@HD\tVN:1.0\tSO:sorted'
    outfile.write(outline+'\n')

    SeenDict={}

    i=0
    samfile = pysam.Samfile(BAM, "rb" )
    for alignedread in samfile.fetch():
        i+=1
        if i % 5000000 == 0:
            print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed', chr, alignedread.pos, ChrDict[chr]
        chr = samfile.getrname(alignedread.tid)
        if SeenDict.has_key(chr):
            continue
        else:
            outline = '@SQ\tSN:' + chr + '\tLN:' + ChrDict[chr]
            print outline
            outfile.write(outline+'\n')
            SeenDict[chr]=''
    outfile.close()

run()

