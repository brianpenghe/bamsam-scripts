##################################
#                                #
# Last modified 12/01/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam

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

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s BAMfilename chrom.sizes outputfilename [-nomulti] [-firstN number_pairs] [-chr chr1,...,chrN] [-regions file chrFiledID leftFieldID rightFieldID]' % argv[0]
        print '\Note: the -regions option and the -chr option will be integrated if both run, i.e. only the regions within the wanted chromosomes will be used'
        sys.exit(1)

    doChr = False
    if '-chr' in argv:
        doChr = True
        chromosomes = argv[argv.index('-chr') + 1].split(',')
        WantedDict = {}
        for chr in chromosomes:
            WantedDict[chr] = ''

    doRegions = False
    if '-regions' in argv:
        doRegions = True
        regionsFile = argv[argv.index('-regions') + 1]
        regionsChr = int(argv[argv.index('-regions') + 2])
        regionsLeft = int(argv[argv.index('-regions') + 3])
        regionsRight = int(argv[argv.index('-regions') + 4])
        linelist = open(regionsFile)
        Regions = []
        for line in linelist:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chr = fields[regionsChr]
            left = int(fields[regionsLeft])
            right = int(fields[regionsRight])
            if doChr:
                if WantedDict.has_key(chr):
                    Regions.append((chr,left,right))
            else:
                Regions.append((chr,left,right))
    
    BAM = argv[1]
    chrominfo=argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        if doChr:
            if WantedDict.has_key(chr):
                chromInfoList.append((chr,start,end))
        else:
            chromInfoList.append((chr,start,end))
    outfilename = argv[3]

    noMulti = False
    if '-nomulti' in argv:
        noMulti = True

    doFirstN = False
    if '-firstN' in argv:
        doFirstN = True
        FN = int(argv[argv.index('-firstN') + 1])

    InsertLengthDistribution = {}
    InsertLengthDistribution['singleton'] = 0

    samfile = pysam.Samfile(BAM, "rb" )

    RN=0
    PN=0
    if doRegions:
        for (chr,start,end) in Regions:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            currentPos=0
            if doFirstN and PN >= FN:
                break
            for alignedread in samfile.fetch(chr, start, end):
                RN+=1
                if RN % 5000000 == 0:
                    print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
                if doFirstN and PN >= FN:
                    break
                try:
                    multiplicity = alignedread.opt('NH')
                except:
                    print 'no NH: tags in BAM file, exiting'
                    sys.exit(1)
                if noMulti and multiplicity > 1:
                    continue
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                pos = alignedread.pos
                if 8 in FLAGfields:
                    InsertLengthDistribution['singleton'] += 1
                    PN+=1
                    continue
                matepos = alignedread.pnext
                if matepos > pos:
                    continue
                IL = pos - matepos + len(alignedread.query)
                if InsertLengthDistribution.has_key(IL):
                    pass
                else:
                    InsertLengthDistribution[IL] = 0
                InsertLengthDistribution[IL] += 1
                PN+=1
    else:
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            currentPos=0
            if doFirstN and PN >= FN:
                break
            for alignedread in samfile.fetch(chr, start, end):
                RN+=1
                if RN % 5000000 == 0:
                    print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
                if doFirstN and PN >= FN:
                    break
                try:
                    multiplicity = alignedread.opt('NH')
                except:
                    print 'no NH: tags in BAM file, exiting'
                    sys.exit(1)
                if noMulti and multiplicity > 1:
                    continue
                fields=str(alignedread).split('\t')
                FLAGfields = FLAG(int(fields[1]))
                pos = alignedread.pos
                if 8 in FLAGfields:
                    InsertLengthDistribution['singleton'] += 1
                    PN+=1
                    continue
                matepos = alignedread.pnext
                if matepos > pos:
                    continue
                IL = pos - matepos + len(alignedread.query)
                if InsertLengthDistribution.has_key(IL):
                    pass
                else:
                    InsertLengthDistribution[IL] = 0
                InsertLengthDistribution[IL] += 1
                PN+=1

    outfile = open(outfilename, 'w')

    outline = '#Length\tNumberPairs'
    outfile.write(outline + '\n')

    keys = InsertLengthDistribution.keys()
    keys.sort()
    for IL in keys:
        outline = str(IL) + '\t' + str(InsertLengthDistribution[IL])
        outfile.write(outline + '\n')

    outfile.close()
            
if __name__ == '__main__':
    main(sys.argv)
