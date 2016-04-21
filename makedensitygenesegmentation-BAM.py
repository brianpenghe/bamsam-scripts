##################################
#                                #
# Last modified 06/04/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import pysam

def countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict):

    reads = 0
    for alignedread in samfile.fetch(chr, rmin, rmax):
        ID=str(alignedread).split('\t')[0]
        if alignedread.is_read1:
            ID = ID + '/1'
        if alignedread.is_read2:
            ID = ID + '/2'
        if noMulti and CorrectionDict[ID] > 1:
            continue
        reads += 1.0/CorrectionDict[ID]
    return reads

def run():

    if len(sys.argv) < 9:
        print 'usage: python %s refGene.txt labelFields chrFieldID leftFieldID rightFieldID strandFieldID BAM chrom.sizes otufilename [-TSSupstream bp bins] [-TSSdownstream bp bins] [-genebodybins bins] [-3UTR bpdownstream bins]' % sys.argv[0]
        print "\tNote: labelFields comma-separated"
        print "\tNote: Genes that are shorter than the -TSSdownstream distance plus the number of genebodybins times the length of TSSdownstream bins will not be considered"
        print "\tDefault values: "
        print "\t\t\tTSSupstreambp = 2000"
        print "\t\t\tTSSupstreambins = 40"
        print "\t\t\tTSSdownstreambp = 1000"
        print "\t\t\tTSSdownstreambins = 20"
        print "\t\t\tgenebodybins = 100"
        print "\t\t\tUTRdownstreambp = 5000"
        print "\t\t\tUTRdownstreambins = 50"

        sys.exit(1)

    genes = sys.argv[1]
    fields=sys.argv[2].split(',')
    LabelFields=[]
    for ID in fields:
        LabelFields.append(int(ID))
    chrFieldID = int(sys.argv[3])
    leftFieldID = int(sys.argv[4])
    rightFieldID = int(sys.argv[5])
    strandFieldID = int(sys.argv[6])
    BAM = sys.argv[7]
    chrominfo=sys.argv[8]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))
    outfilename = sys.argv[9]

    TSSupstreambp = 2000
    TSSupstreambins = 40
    TSSdownstreambp = 1000
    TSSdownstreambins = 20
    genebodybins = 100
    UTRdownstreambp = 5000
    UTRdownstreambins = 50

    if '-TSSupstream' in sys.argv:
        TSSupstreambp = int(sys.argv[sys.argv.index('-TSSupstream') + 1])
        TSSupstreambins = int(sys.argv[sys.argv.index('-TSSupstream') + 2])
    if '-TSSdownstream' in sys.argv:
        TSSdownstreambp = int(sys.argv[sys.argv.index('-TSSdownstream') + 1])
        TSSdownstreambins = int(sys.argv[sys.argv.index('-TSSdownstream') + 2])
    if '-genebodybins' in sys.argv:
        genebodybins = int(sys.argv[sys.argv.index('-genebodybins')+1])
    if '-zeros' in sys.argv:
        doZeros=True
        zeros = float(sys.argv[sys.argv.index('-zeros')+1])
    if '-3UTR' in sys.argv:
        UTRdownstreambp = int(sys.argv[sys.argv.index('-3UTR')+1])
        UTRdownstreambins = int(sys.argv[sys.argv.index('-3UTR')+2])

    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti = True

    ReadDict={}

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
    normalizeBy = readNumber/1000000.
    print 'RPM normalization Factor =', normalizeBy
   
    outfile = open(outfilename, 'w')

    linelist= open(genes)
    j=1
    for line in linelist:
        if line.startswith('#'):
            continue
        if j % 1000 == 0:
            print j
        j+=1
        fields = line.strip().split('\t')
        chr = fields[chrFieldID]
        left = int(fields[leftFieldID])
        right = int(fields[rightFieldID])
        try:
            for alignedread in samfile.fetch(chr, 0, 100):
                a='b'
        except:
            continue
        strand = fields[strandFieldID]
        if right - left < TSSdownstreambp + TSSdownstreambins*genebodybins:
            continue
        values = []
        if strand == 'F' or strand == '+':
            current=left-TSSupstreambp
            slide=float(TSSupstreambp)/TSSupstreambins
            for i in range(TSSupstreambins-1):
                rmin=int(current)
                rmax=int(current+slide)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current+slide
            current=left
            slide=float(TSSdownstreambp)/TSSdownstreambins
            for i in range(TSSdownstreambins-1):
                rmin=int(current)
                rmax=int(current+slide)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current+slide
            current=left+TSSdownstreambp
            slide=float(right-current)/genebodybins
            for i in range(genebodybins-1):
                rmin=int(current)
                rmax=int(current+slide)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current+slide
            current=right
            slide=float(UTRdownstreambp)/UTRdownstreambins
            for i in range(UTRdownstreambins-1):
                rmin=int(current)
                rmax=int(current+slide)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current+slide
        if strand == 'F' or strand == '-':
            start=right+TSSupstreambp
            current=start
            slide=float(TSSupstreambp)/TSSupstreambins
            for i in range(TSSupstreambins-1):
                rmin=int(current-slide)
                rmax=int(current)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current-slide
            current=right
            slide=float(TSSdownstreambp)/TSSdownstreambins
            for i in range(TSSdownstreambins-1):
                rmin=int(current-slide)
                rmax=int(current)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current-slide
            current=right-TSSdownstreambp
            slide=float(current-left)/genebodybins
            for i in range(genebodybins-1):
                rmin=int(current-slide)
                rmax=int(current)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current-slide
            current=left
            slide=float(UTRdownstreambp)/UTRdownstreambins
            for i in range(UTRdownstreambins-1):
                rmin=int(current-slide)
                rmax=int(current)
                values.append(countReads(samfile,chr,rmin,rmax,noMulti,CorrectionDict)/(normalizeBy*slide/1000))
                current=current-slide
        outline = ''
        for geneID in LabelFields:
            outline = outline + fields[geneID] + '\t'
        outline = outline + chr + '\t' + str(left) + '\t' + str(right) + '\t' + strand
        for value in values:
             outline = outline + '\t' + str(value)
        outfile.write(outline + '\n')

    outfile.close()

run()
