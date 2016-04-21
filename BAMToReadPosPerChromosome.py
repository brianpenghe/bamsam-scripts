##################################
#                                #
# Last modified 02/23/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import os

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

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAMfilename chrom.sizes outputfile_prefix [-nomulti] [-mismatchesMD M] [-mismatches M] [-readLength min max] [-chr chrN1(,chrN2....)] [-uniqueBAM] [-noNH samtools]' % sys.argv[0]
        print '\tUse the -mismatches option to specify the maximum number of mismatches allowed for an alignment to be considered; use the -mimatchesMD option is mismatches are specified with the MD special tag'
        print '\tuse the -noNH option and supply a path to samtools in order to have the file converted to one that has NH tags'
        print '\tNote: Right now the script does not work with gapped alignments'
        sys.exit(1)
    
    BAM = sys.argv[1]
    chrominfo=sys.argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
    outfileprefix = sys.argv[3]

    doUniqueBAM = False
    if '-uniqueBAM' in sys.argv:
        print 'will treat all alignments as unique'
        doUniqueBAM = True
        TotalReads = 0
        pass

    samfile = pysam.Samfile(BAM, "rb" )
    try:
        print 'testing for NH tags presence'
        for alignedread in samfile.fetch():
            multiplicity = alignedread.opt('NH')
            print 'file has NH tags'
            break
    except:
        if '-noNH' in sys.argv:
            print 'no NH: tags in BAM file, will replace with a new BAM file with NH tags'
            samtools = sys.argv[sys.argv.index('-noNH')+1]
            BAMpreporcessingScript = sys.argv[0].rpartition('/')[0] + '/bamPreprocessing.py'
            cmd = 'python ' + BAMpreporcessingScript + ' ' + BAM + ' ' + BAM + '.NH'
            os.system(cmd)
            cmd = 'rm ' + BAM
            os.system(cmd)
            cmd = 'mv ' + BAM + '.NH' + ' ' + BAM
            os.system(cmd)
            cmd = samtools + ' index ' + BAM
            os.system(cmd)
        elif doUniqueBAM:
            pass
        else:
            print 'no NH: tags in BAM file, exiting'
            sys.exit(1)

    doReadLength=False
    if '-readLength' in sys.argv:
        doReadLength=True
        minRL = int(sys.argv[sys.argv.index('-readLength')+1])
        maxRL = int(sys.argv[sys.argv.index('-readLength')+2])
        print 'will only consider reads between', minRL, 'and', maxRL, 'bp length'

    doMaxMMMD=False
    if '-mismatchesMD' in sys.argv:
        doMaxMMMD=True
        maxMM = int(sys.argv[sys.argv.index('-mismatchesMD')+1])
        print 'Will only consider alignments with', maxMM, 'or less mismatches'

    doMaxMM=False
    if '-mismatches' in sys.argv:
        doMaxMM=True
        maxMM = int(sys.argv[sys.argv.index('-mismatches')+1])
        print 'Will only consider alignments with', maxMM, 'or less mismatches'

    doChrSubset=False
    if '-chr' in sys.argv:
        doChrSubset=True
        WantedChrDict={}
        for chr in sys.argv[sys.argv.index('-chr')+1].split(','):
            WantedChrDict[chr]=''

    noMulti=False
    if '-nomulti' in sys.argv:
        print 'will only consider unique alignments'
        noMulti=True

    doRPM=False
    if '-RPM' in sys.argv:
        doRPM=True

    doStranded=False
    if '-stranded' in sys.argv:
        doStranded=True
        strand=sys.argv[sys.argv.index('-stranded')+1]
        print 'will only consider', strand, 'strand reads'

    samfile = pysam.Samfile(BAM, "rb" )

    RN=0
    for (chr,start,end) in chromInfoList:
        if doChrSubset:
            if WantedChrDict.has_key(chr):
                pass
            else:
                continue
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            print 'region', chr,start,end, 'not found in bam file, skipping'
            continue
        outfileplus = open(outfileprefix + '.' + chr + '.plus', 'w')
        outfileminus = open(outfileprefix + '.' + chr + '.minus', 'w')
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 5000000 == 0:
                print 'counting total number of reads', str(RN/1000000) + 'M alignments processed', chr, currentPos, end
            if doReadLength:
                if len(alignedread.seq) > maxRL or len(alignedread.seq) < minRL:
                    continue
            if doMaxMM:
                mismatches = 0
                for (m,bp) in alignedread.cigar:
                    if m == 8:
                        mismatches+=1
                if mismatches > maxMM:
                    continue
            if doMaxMMMD:
                MM = alignedread.opt('MD')
                mismatches = 0
                if MM.isdigit():
                    pass
                else:
                    for s in range(len(MM)):
                        if MM[s].isalpha():
                            mismatches+=1
                if mismatches > maxMM:
                    continue
            if doUniqueBAM:
                multiplicity = 1
                continue
            else:
                multiplicity = alignedread.opt('NH')
            if noMulti and multiplicity > 1:
                continue
            fields=str(alignedread).split('\t')
            ID=fields[0]
            FLAGfields = FLAG(int(fields[1]))
            if 16 in FLAGfields:
                s = '-'
            else:
                s = '+'
            if doStranded and s!=strand:
                continue
            currentPos=alignedread.pos
            endPos = currentPos
            for (m,bp) in alignedread.cigar:
                if m == 0:
                    endPos+=bp
            if s == '+':
                outfileplus.write(str(currentPos) + '\n')
            if s == '-':
                outfileminus.write(str(endPos) + '\n')
        outfileplus.close()
        outfileminus.close()

run()
