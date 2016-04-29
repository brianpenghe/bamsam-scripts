##################################
#                                #
# Last modified 02/17/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

# FLAG field meaning
# 0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair
# 0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
# 0x0004 4 the query sequence itself is unmapped
# 0x0008 8 the mate is unmapped 1
# 0x0010 16 strand of the query (0 for forward; 1 for reverse strand)
# 0x0020 32 strand of the mate 1
# 0x0040 64 the read is the first read in a pair 1,2
# 0x0080 128 the read is the second read in a pair 1,2
# 0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records)
# 0x0200 512 the read fails platform/vendor quality checks
# 0x0400 1024 the read is either a PCR duplicate or an optical duplicate

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024]

    FLAGList=[]

    MaxNumberList=[]
    for i in Numbers:
        if i <= FLAG:
            MaxNumberList.append(i)

    Residual=FLAG
    maxPos = len(MaxNumberList)-1

    while Residual > 0:
        if MaxNumberList[maxPos] <= Residual:
            Residual = Residual - MaxNumberList[maxPos]
            FLAGList.append(MaxNumberList[maxPos])
            maxPos-=1
        else:
            maxPos-=1
  
    return FLAGList

import sys
import pysam
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s BAM chrom.sizes' % argv[0]
        print '       BAM file has to be indexed'
        sys.exit(1)

    SAM = argv[1]

    CoverageDict={}
    CoverageDict['+'] = {}
    CoverageDict['-'] = {}

    chrominfo=argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
        CoverageDict['+'][chr]={}
        CoverageDict['-'][chr]={}

    ReadMultiplicityDict={}

    print 'examining read multiplicty'

    samfile = pysam.Samfile(SAM, "rb" )
    for (chr,start,end) in chromInfoList:
        try:
            has_reads = False
            for alignedread in samfile.fetch(chr, start, end):
                has_reads = True
                if has_reads:
                    break
        except:
            continue
        for alignedread in samfile.fetch(chr, start, end):
            fields=str(alignedread).split('\t')
            if len(fields[3].split('),')) > 1:
                continue
            ID=fields[0]
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if ReadMultiplicityDict.has_key(ID):
                pass
            else:
                ReadMultiplicityDict[ID]=0
            ReadMultiplicityDict[ID]+=1

    print 'counting distinct reads'

    for (chr,start,end) in chromInfoList:
        try:
            has_reads = False
            for alignedread in samfile.fetch(chr, start, end):
                has_reads = True
                if has_reads:
                    break
        except:
            continue
        for alignedread in samfile.fetch(chr, start, end):
            fields=str(alignedread).split('\t')
            if 16 in FLAG(alignedread.flag):
                strand = '-'
            else:
                strand = '+'
            ID=fields[0]
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if ReadMultiplicityDict.has_key(ID):
                if ReadMultiplicityDict[ID] > 1:
                    continue
                else:
                    pass
            else:
                continue
            start=alignedread.pos
            if CoverageDict[strand][chr].has_key(start):
                pass
            else:
                CoverageDict[strand][chr][start]=0
            CoverageDict[strand][chr][start]+=1
 
    TotalUniqueReads = 0
    DistinctUniqueReads = 0

    for (chr,start,end) in chromInfoList:
        for pos in CoverageDict['+'][chr].keys():
            DistinctUniqueReads +=1
            TotalUniqueReads += CoverageDict['+'][chr][pos]
        for pos in CoverageDict['-'][chr].keys():
            DistinctUniqueReads +=1
            TotalUniqueReads += CoverageDict['-'][chr][pos]

    outline = '#File\tTotal Unique Reads\tDistinct Unique Reads\tComplexity'
    print outline
    outline = SAM + '\t' + str(TotalUniqueReads) + '\t' + str(DistinctUniqueReads) + '\t' + str(DistinctUniqueReads/(TotalUniqueReads + 0.0))[0:4]
    print outline

if __name__ == '__main__':
    main(sys.argv)
