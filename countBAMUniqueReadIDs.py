##################################
#                                #
# Last modified 03/21/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAM chrom.size outputfilename' % sys.argv[0]
        print '       BAM file has to be indexed'
        sys.exit(1)

    BAM = sys.argv[1]
    outputfilename = sys.argv[3]
    chrominfo=sys.argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))

    ReadDict={}

    i=0
    samfile = pysam.Samfile(BAM, "rb" )
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print str(i/1000000) + 'M alignments processed', chr,start,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                ReadDict[ID]=0
        except:
            print 'skipping', chr, start, end, ' - no alignments found in that range'

    outfile=open(outputfilename, 'w')

    outline=BAM + '\t'+str(len(ReadDict.keys()))
    print outline
    outfile.write(outline+'\n')
             
    outfile.close()

run()
