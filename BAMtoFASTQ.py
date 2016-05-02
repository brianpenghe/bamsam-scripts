##################################
#                                #
# Last modified 12/12/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
import string
from sets import Set
import os

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for i in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-i-1]]
    return sequence

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s BAM chrom.sizes outfilename [-ignoreChromSizes] [-OneEndOnly 1 2] [-alignedOnly] [-chr chrN1[,....,chrNN]] [-withinRegions regions_file chrFieldID leftFieldID rightFieldID] [-noNH samtools]' % argv[0]
        print '\tNote: the script assumes no duplicate readIDs'
        print '\t-withinRegions option does not work together with the -ignoreChromSizes option'
        print '\tuse the -ignoreChromSizes option to also get the unaligned reads'
        print '\tthe -chr option will output all alignments on the set of chromosomes indicate and will override the -withinRegions option'
        print '\tuse the -noNH option and supply a path to samtools in order to have the file converted to one that has NH tags'
        print '\tuse - for outfilename if you want to print to standard output'
        sys.exit(1)

    BAM = argv[1]
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

    doStdOut = False
    if outputfilename == '-':
        doStdOut = True

    doEnd = False
    if '-OneEndOnly' in argv:
        doEnd = True
        End = int(argv[argv.index('-OneEndOnly')+1])

    doChr=False
    if '-chr' in argv:
        doChr=True
        WantedChrDict={}
        for chr in argv[argv.index('-chr')+1].split(','):
            WantedChrDict[chr]=''
        if not doStdOut:
            print 'will output all reads aligning to the following chromosomes:'
            print WantedChrDict.keys()

    doIgnoreChromSizes=False
    if '-ignoreChromSizes' in argv:
        doIgnoreChromSizes=True
        if not doStdOut:
            print 'will ignore chrom.sizes file'

    doAlignedOnly=False
    if '-alignedOnly' in argv:
        doAlignedOnly=True
        if not doStdOut:
            print 'will only output aligned reads'

    if not doStdOut:
        print 'will write reads into: ', outputfilename

    ReadDict={}

    samfile = pysam.Samfile(BAM, "rb" )
    try:
        for read in samfile.fetch():
            multiplicity = read.opt('NH')
            if not doStdOut:
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
        else:
            print 'no NH: tags in BAM file, exiting'
            sys.exit(1)

    if not doStdOut:
        outfile = open(outputfilename, "w")

    i=0
    samfile = pysam.Samfile(BAM, "rb" )
    if doIgnoreChromSizes:
        for read in samfile.fetch(until_eof=True):
            i+=1
            if i % 5000000 == 0:
                if not doStdOut:
                    print 'examining read multiplicity and inputting reads', str(i/1000000) + 'M alignments processed processed'
            fields=str(read).split('\t')
            if doAlignedOnly:
                try:
                    chr = samfile.getrname(read.tid)
                except:
                    print fields
                    sys.exit(1)
            ID = read.qname
            if read.is_read1:
                ID = ID + '/1'
                if doEnd and End == 2:
                    continue
            if read.is_read2:
                ID = ID + '/2'
                if doEnd and End == 1:
                    continue
            sequence = read.seq
            quality = read.qual
            if read.is_reverse:
                sequence = getReverseComplement(sequence)
                quality=quality[::-1]
            if read.is_unmapped:
                if not doStdOut:
                    outfile.write('@' + ID + '\n')
                    outfile.write(sequence + '\n')
                    outfile.write('+\n')
                    outfile.write(quality + '\n')
                else:
                    print '@' + ID
                    print sequence
                    print '+'
                    print quality
                continue
            if read.opt('NH') == 1:
                if not doStdOut:
                    outfile.write('@' + ID + '\n')
                    outfile.write(sequence + '\n')
                    outfile.write('+\n')
                    outfile.write(quality + '\n')
                else:
                    print '@' + ID
                    print sequence
                    print '+'
                    print quality
                continue
            if ReadDict.has_key(ID):
                if ReadDict[ID] != (sequence,quality):
                    sequence = getReverseComplement(sequence)
                    quality=quality[::-1]
                    if ReadDict[ID] != (sequence,quality):
                        print ID, sequence, quality
                        print ReadDict[ID]
                        print 'duplicate read ID detected, exiting'
                        sys.exit(1)
                else:
                    pass
            else:
                ReadDict[ID] = (sequence,quality)
    elif doChr:
        for (chr,start,end) in chromInfoList:
            if WantedChrDict.has_key(chr):
                pass
            else:
                continue
            try:
                for read in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                continue
            for read in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    if not doStdOut:
                        print str(i/1000000) + 'M alignments processed processed', chr,start,read.pos,end
                fields=str(read).split('\t')
                ID = read.qname
                if read.is_read1:
                    ID = ID + '/1'
                    if doEnd and End == 2:
                        continue
                if read.is_read2:
                    ID = ID + '/2'
                    if doEnd and End == 2:
                        continue
                sequence = read.seq
                quality = read.qual
                if read.is_reverse:
                    sequence = getReverseComplement(sequence)
                    quality=quality[::-1]
                if read.opt('NH') == 1:
                    if not doStdOut:
                        outfile.write('@' + ID + '\n')
                        outfile.write(sequence + '\n')
                        outfile.write('+\n')
                        outfile.write(quality + '\n')
                    else:
                        print '@' + ID
                        print sequence
                        print '+'
                        print quality
                    continue
                if ReadDict.has_key(ID):
                    if ReadDict[ID] != (sequence,quality):
                        print ID, sequence, qualtiy
                        print ReadDict[ID]
                        print 'duplicate read ID detected, exiting'
                        sys.exit(1)
                    else:
                        pass
                else:
                    ReadDict[ID] = (sequence,quality)
    elif '-withinRegions' in argv:
        regionDict={}
        if not doStdOut:
            print 'will only output reads with alignments within regions defined in', argv[argv.index('-withinRegions')+1]
        linelist = open(argv[argv.index('-withinRegions')+1])
        chrFieldID = int(argv[argv.index('-withinRegions')+2])
        leftFieldID = int(argv[argv.index('-withinRegions')+3])
        rightFieldID = int(argv[argv.index('-withinRegions')+4])
        for line in linelist:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chr = fields[chrFieldID]
            if regionDict.has_key(chr):
                pass
            else:
                regionDict[chr]={}
            left = int(fields[leftFieldID])
            right = int(fields[rightFieldID])
            test = 0
            try:
                for read in samfile.fetch(chr, left, right):
                    test+=1
                    if test == 1:
                        break
            except:
                continue
            for read in samfile.fetch(chr, left, right):
                i+=1
                if i % 5000000 == 0:
                    if not doStdOut:
                        print 'examining read multiplicity and inputting reads', str(i/1000000) + 'M alignments processed processed', chr,start,read.pos,end
                fields=str(read).split('\t')
                ID = read.qname
                if read.is_read1:
                    ID = ID + '/1'
                    if doEnd and End == 2:
                        continue
                if read.is_read2:
                    ID = ID + '/2'
                    if doEnd and End == 2:
                        continue
                sequence = read.seq
                quality = read.qual
                if read.is_reverse:
                    sequence = getReverseComplement(sequence)
                    quality=quality[::-1]
                if read.opt('NH') == 1:
                    if not doStdOut:
                        outfile.write('@' + ID + '\n')
                        outfile.write(sequence + '\n')
                        outfile.write('+\n')
                        outfile.write(quality + '\n')
                    else:
                        print '@' + ID
                        print sequence
                        print '+'
                        print quality
                    continue
                if ReadDict.has_key(ID):
                    if ReadDict[ID] != (sequence,quality):
                        print ID, sequence, qualtiy
                        print ReadDict[ID]
                        print 'duplicate read ID detected, exiting'
                        sys.exit(1)
                    else:
                        pass
                else:
                    ReadDict[ID] = (sequence,quality)
    else:
        for (chr,start,end) in chromInfoList:
            try:
                for read in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                continue
            for read in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    if not doStdOut:
                        print 'examining read multiplicity and inputting reads', str(i/1000000) + 'M alignments processed processed', chr,start,read.pos,end
                fields=str(read).split('\t')
                ID = read.qname
                if read.is_read1:
                    ID = ID + '/1'
                    if doEnd and End == 2:
                        continue
                if read.is_read2:
                    ID = ID + '/2'
                    if doEnd and End == 2:
                        continue
                sequence = read.seq
                quality = read.qual
                if read.is_reverse:
                    sequence = getReverseComplement(sequence)
                    quality=quality[::-1]
                if read.opt('NH') == 1:
                    if not doStdOut:
                        outfile.write('@' + ID + '\n')
                        outfile.write(sequence + '\n')
                        outfile.write('+\n')
                        outfile.write(quality + '\n')
                    else:
                        print '@' + ID
                        print sequence
                        print '+'
                        print quality
                    continue
                if ReadDict.has_key(ID):
                    if ReadDict[ID] != (sequence,quality):
                        print ID, sequence, qualtiy
                        print ReadDict[ID]
                        print 'duplicate read ID detected, exiting'
                        sys.exit(1)
                    else:
                        pass
                else:
                    ReadDict[ID] = (sequence,quality)

    if not doStdOut:
        print i

    if not doStdOut:
        for readID in ReadDict.keys():
            (sequence,quality) = ReadDict[readID]
            outfile.write('@' + readID + '\n')
            outfile.write(sequence + '\n')
            outfile.write('+\n')
            outfile.write(quality + '\n')
        outfile.close()
    else:
        for readID in ReadDict.keys():
            (sequence,quality) = ReadDict[readID]
            print '@' + readID
            print sequence
            print '+'
            print quality

if __name__ == '__main__':
    main(sys.argv)

