##################################
#                                #
# Last modified 09/05/2013       # 
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

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s BAMfilename chrom.sizes outputfilename [-uniqueBAM] [-nomulti] [-noNH samtools]' % argv[0]
        sys.exit(1)

    BAM = argv[1]
    chrominfo=argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
    outfilename = argv[3]

    noMulti=False
    if '-nomulti' in argv:
        print 'will only consider unique alignments'
        noMulti=True

    doUniqueBAM = False
    if '-uniqueBAM' in argv:
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
        if '-noNH' in argv:
            print 'no NH: tags in BAM file, will replace with a new BAM file with NH tags'
            samtools = argv[argv.index('-noNH')+1]
            BAMpreporcessingScript = argv[0].rpartition('/')[0] + '/bamPreprocessing.py'
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

    outfile = open(outfilename, 'w')

    RN=0
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            print 'region', chr,start,end, 'not found in bam file, skipping'
            continue
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 5000000 == 0:
                print str(RN/1000000) + 'M alignments processed', chr, pos, end
            fields=str(alignedread).split('\t')
            FLAGfields = FLAG(int(fields[1]))
            ID = fields[0]
            if doUniqueBAM:
                multiplicity = 1
            else:
                multiplicity = alignedread.opt('NH')
            if noMulti and multiplicity > 1:
                continue
            scroe = 1000.0/multiplicity
            FLAGfields = FLAG(int(fields[1]))
            if alignedread.is_reverse:
                strand = '-'
            else:
                strand = '+'
            pos = alignedread.pos
            seq = alignedread.seq
            qual = alignedread.qual
            outline = chr + '\t' + str(pos) + '\t' + str(pos + len(seq)) + '\t' + seq + '\t' + qual + '\t' + strand
            outfile.write(outline + '\n')

    outfile.close()
            
if __name__ == '__main__':
    main(sys.argv)
