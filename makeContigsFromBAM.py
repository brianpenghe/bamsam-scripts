##################################
#                                #
# Last modified 08/29/2011       # 
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

def WriteContig(contig,outfile):

    outline = contig[0] + '\t' + str(contig[1]) + '\t' + str(contig[2]) + '\t' + contig[3]  + '\t' + str(contig[2]-contig[1]) + '\t' + str(contig[4])
    outfile.write(outline+'\n')

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

    if len(argv) < 3:
        print 'usage: python %s BAM chrom.sizes outputfilename [-stranded] [-splices] [-withmulti]' % argv[0]
        print '       BAM file has to be indexed'
        print '       paired reads are not taken advantage of right now'
        print '       use the -splices option if you want spliced reads to be counted too; the spliced read will be assigned to the leftmost exon it covers'
        print '       use the -withmulti option if you want multi reads to be counted too'
        print '       the scripts needs an NH attribute in all alignments'
        sys.exit(1)

    SAM = argv[1]
    outputfilename = argv[3]
    chrominfo=argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))

    doSpliced=False
    doMulti=False
    doStranded=False

    if '-stranded' in argv:
        doStranded=True

    if '-withmulti' in argv:
        doMulti=True

    if '-splices' in argv:
        doSpliced=True

    outfile=open(outputfilename, 'w')

    outline='#chr\tstart\tend\tstrand\tlength\treads\t'
    outfile.write(outline+'\n')

    i=0
    samfile = pysam.Samfile(SAM, "rb" )
    for (chr,start,end) in chromInfoList:
        print chr,start,end
        contig=[chr,0,0,'.',0]
        if doStranded:
            pluscontig=[chr,0,0,'+',0]
            minuscontig=[chr,0,0,'-',0]
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print str(i/1000000) + 'M alignments processed'
            fields=str(alignedread).split('\t')
            if "'NH'," in fields[8]:
                scaleby = float(fields[8].split("'NH', ")[1].split(')')[0])
            else:
                print 'multireads not specified with the NH tag, exiting'
                sys.exit(1)
            if doMulti:
                pass
            else:
                if scaleby > 1:
                    continue
            chr=samfile.getrname(alignedread.rname)
            start=alignedread.pos
            if doSpliced:
                end=start + alignedread.cigar[0][1]
            else:
                if len(alignedread.cigar) > 1:
                    continue
                end=start + len(fields[6])
            if doStranded:
                FLAGfields = FLAG(int(fields[1]))
                if 16 in FLAGfields:
                    strand = '-'
                    if pluscontig == [chr,0,0,'+',0]:
                        pluscontig = [chr,start,end,strand,scaleby]
                        continue
                    if start >= pluscontig[1] and start <= pluscontig[2]:
                        pluscontig[2] = end
                        pluscontig[4] = pluscontig[4] + 1
                    else:
                        WriteContig(pluscontig,outfile)
                        pluscontig=[chr,start,end,strand,1]
                else:
                    strand = '+'
                    if minuscontig == [chr,0,0,'-',0]:
                        minuscontig = [chr,start,end,strand,scaleby]
                        continue
                    if start >= minuscontig[1] and start <= minuscontig[2]:
                        minuscontig[2] = end
                        minuscontig[4] = minuscontig[4] + 1
                    else:
                        WriteContig(minuscontig,outfile)
                        minuscontig=[chr,start,end,strand,1]
            else:
                strand = '.'
                if contig == [chr,0,0,'.',0]:
                    contig = [chr,start,end,strand,scaleby]
                    continue
                if start >= contig[1] and start <= contig[2]:
                    contig[2] = end
                    contig[4] = contig[4] + 1
                else:
                    WriteContig(contig,outfile)
                    contig=[chr,start,end,strand,1]
        if doStranded:
            WriteContig(pluscontig,outfile)
            WriteContig(minuscontig,outfile)
        else:
            WriteContig(contig,outfile)

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
