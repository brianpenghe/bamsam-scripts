##################################
#                                #
# Last modified 05/21/2012       # 
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

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s title BAMfilename chrom.sizes outputfilename [-stranded + | -] [-nomulti] [-RPM] [-notitle] [-singlebasepair] [-mismatchesMD M] [-mismatches M] [-chr chrN1(,chrN2....)]' % sys.argv[0]
        print '       Use the -mismatches option to specify the maximum number of mismatches allowed for an alignment to be considered; use the -mimatchesMD option is mismatches are specified with the MD special tag'
        sys.exit(1)
    
    doSingleBP=False
    if '-singlebasepair' in sys.argv:
        doSingleBP=True

    doTitle=True
    if '-notitle' in sys.argv:
        doTitle=False

    title = sys.argv[1]
    BAM = sys.argv[2]
    chrominfo=sys.argv[3]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
    outfilename = sys.argv[4]

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

    CorrectionDict={}
    i=0
    unique=0
    multi=0
    samfile = pysam.Samfile(BAM, "rb" )
    for (chr,start,end) in chromInfoList:
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            print 'region', chr,start,end, 'not found in bam file, skipping'
            continue
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed', len(CorrectionDict.keys()), 'reads found', chr, start, alignedread.pos, end
            fields=str(alignedread).split('\t')
            ID=fields[0]
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
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            if CorrectionDict.has_key(ID):
                CorrectionDict[ID]+=1
            else:
                CorrectionDict[ID]=1

    readNumber = len(CorrectionDict.keys())
    print 'found', readNumber, 'reads'
    normFactor = readNumber/1000000.
    print 'RPM normalization Factor =', normFactor

    outfile = open(outfilename, 'w')
    if doTitle:
        outline='track type=bedGraph name="' + title + '"'
        outfile.write(outline+'\n')

    RN=0
    for (chr,start,end) in chromInfoList:
        coverageDict={}
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
        currentPos=0
        for alignedread in samfile.fetch(chr, start, end):
            RN+=1
            if RN % 5000000 == 0:
                print str(RN/1000000) + 'M alignments processed', chr, currentPos, end
            fields=str(alignedread).split('\t')
            ID = fields[0]
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
            if alignedread.is_read1:
                ID = ID + '/1'
            if alignedread.is_read2:
                ID = ID + '/2'
            multiplicity = CorrectionDict[ID]
            if noMulti and multiplicity > 1:
                continue
            scaleby=1.0/multiplicity
            if doStranded:
                FLAGfields = FLAG(int(fields[1]))
                if 16 in FLAGfields:
                    s = '-'
                else:
                    s = '+'
                if s!=strand:
                    continue
            currentPos=alignedread.pos
            for (m,bp) in alignedread.cigar:
                if m == 0:
                    for j in range(currentPos,currentPos+bp):
                        if coverageDict.has_key(j+1):
                            coverageDict[j+1]+=scaleby
                        else:
                            coverageDict[j+1]=scaleby
                elif m == 2:
                    pass
                elif m == 3:
                    pass
                else:
                    continue
                currentPos=currentPos+bp
        posKeys=coverageDict.keys()
        posKeys.sort()
        if len(posKeys) == 0:
            continue
        initial=(posKeys[0],coverageDict[posKeys[0]])
        previous=(posKeys[0],coverageDict[posKeys[0]])
        written=['']
        if doSingleBP:
            for i in range(1,max(posKeys)+1):
                if coverageDict.has_key(i):
                    if doStranded and strand == '-':
                        if doRPM:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t-' + str(coverageDict[i]/normFactor)
                        else:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t-' + str(coverageDict[i])
                    else:
                        if doRPM:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(coverageDict[i]/normFactor)
                        else:
                            outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(coverageDict[i])
                    outfile.write(outline+'\n')
                else:
                    outline = chr + '\t' + str(i-1) + '\t' + str(i+1-1) + '\t' + str(0)
                    outfile.write(outline+'\n')
        else:
            for i in posKeys[1:len(posKeys)]:
                if previous[0]+1 == i and previous[1]==coverageDict[i]:
                     previous=(i,coverageDict[i])
                else:
                     if written[0]==initial[0]:
                         print '####', written, initial, previous
                     if doStranded and strand == '-':
                         if doRPM:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
                         else:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t-'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
                     else:
                         if doRPM:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]/normFactor).split('.')[0] + '.' + str(initial[1]/normFactor).split('.')[1][0:4]
                         else:
                             outline=chr+'\t'+str(initial[0]-1)+'\t'+str(previous[0]+1-1)+'\t'+str(initial[1]).split('.')[0] + '.' + str(initial[1]).split('.')[1][0:4]
                     written=(initial[0],previous[0]+1)
                     outfile.write(outline+'\n')
                     initial=(i,coverageDict[i])
                     previous=(i,coverageDict[i])

    outfile.close()
            
run()
