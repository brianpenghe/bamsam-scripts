##################################
#                                #
# Last modified 05/04/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
import string
from sets import Set

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s BAM chrom.sizes outfilename' % argv[0]
        sys.exit(1)

    BAM = argv[1]
    outputfilename = argv[3]
    chrominfo=argv[2]

    chromInfoList=[]
    ChrDict={}
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
        ChrDict[chr]={}
        ChrDict[chr]['unique']=0
        ChrDict[chr]['multi']=0

    ReadMultiplicity={}

    i=0
    samfile = pysam.Samfile(BAM, "rb" )
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            continue
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print 'examining read multiplicity and inputting reads', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
            ID = alignedread.qname
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if ReadMultiplicity.has_key(ID):
                ReadMultiplicity[ID]+=1
            else:
                ReadMultiplicity[ID]=1

    i=0
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            continue
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print 'counting alignments', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
            ID = alignedread.qname
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if ReadMultiplicity[ID]==1:
                ChrDict[chr]['unique']+=1
            else:
                ChrDict[chr]['multi']+=1

    outfile = open(outputfilename, "w")

    outfile.write('#chr\tlength\tunique\tmulti\n')

    for (chr,start,end) in chromInfoList:
        outline = chr + '\t' + str(end) + '\t' + str(ChrDict[chr]['unique']) + '\t' + str(ChrDict[chr]['multi']) 
        outfile.write(outline + '\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

