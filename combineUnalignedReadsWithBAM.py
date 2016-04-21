##################################
#                                #
# Last modified 11/08/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
import string
import os
import subprocess

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s samtools BAM unaligned_reads [-compression bzip2 | gunzip | none]' % sys.argv[0]
        print '\tNote1: BAM file has to be indexed'
        print '\tNote2: the script will check whether the unaligned reads are compressed or not, but it does so relying on their name; if they end with .gz, gunzip will be used, if they end with .bz2, bzip2 will be used, in any other case, an uncompressed file will be assumed; use the -compression option if you want to specify the compression method directly'
        print '\tNote3: the script will print SAM aligments into standard output; couple it with a samtools run into a new BAM file'
        sys.exit(1)

    samtools = sys.argv[1]
    BAM = sys.argv[2]
    unaligned = sys.argv[3]
#    outputfilename = sys.argv[3]

    doCompression = False
    if '-compression' in sys.argv:
        doCompression = True
        compression = sys.argv[sys.argv.index('-compression')+1]

#    outfile = open(outputfilename,'w')

#    samfile = pysam.Samfile(BAM, "rb" )
#    outfile = pysam.Samfile(outputfilename, "wb", template=samfile)
#    i=0
#    for read in samfile.fetch(until_eof=True):
#        i+=1
#        if i % 5000000 == 0:
#            print str(i/1000000) + 'M alignments processed'
#        outfile.write(read)


    cmd1 = samtools + ' view -H ' + BAM
    p1 = os.popen(cmd1, "r")
    line = 0
    while line != '':
        line = p1.readline()
        if line != '':
            print line.strip()

    cmd1 = samtools + ' view ' + BAM
    p1 = os.popen(cmd1, "r")
    line = 0
    while line != '':
        line = p1.readline()
        if line != '':
            print line.strip()

    if doCompression:
        if compression == 'gunzip':
            cmd1 = 'gunzip -c ' + unaligned
        if compression == 'bzip2':
            cmd1 = 'bzip2 -cd ' + unaligned
        if compression == 'none':
            cmd1 = 'cat ' + unaligned
    else:
        if unaligned.endswith('.gz'):
            cmd1 = 'gunzip -c' + unaligned
        elif unaligned.endswith('.bz2'):
            cmd1 = 'bzip2 -cd' + unaligned
        else:
            cmd1 = 'cat ' + unaligned

    p1 = os.popen(cmd1, "r")
    i=0
    i=0
    pos=1
    line = 0
    while line != '':
        i+=1
        line = p1.readline()
        scoresNext=False
        seqNext=False
        if pos==1:
            if line.startswith('@'):
                ID = line.strip()[0:-1]
                pos=2
                continue
            else:
                if line == '':
                    continue
                else:
                    print 'invalid read', line, 'line', i
                    break
        if pos==2:
            if i % 10000000 == 0:
                print str(i/4000000) + 'M reads processed'
            sequence = line.strip()
            pos=3
            continue
        if pos==3:
            if line.startswith('+'):
                pos=4
                continue
            else:
                print 'invalid read', line, 'line', i
                break
        if pos==4:
            quality = line.strip()
            pos=1
            line = ID + '\t4\t*\t0\t0\t*\t*\t0\t0\t' + sequence + '\t' + quality + '\tXM:i:0\tNH:i:1'
#            a = pysam.AlignedRead()
#            a.qname = ID
#            a.seq = sequence
#            a.flag = 4
#            a.rname = None
#            a.pos = 0
#            a.mapq = 0
#            a.cigar = ( (0,10), (2,1), (0,25) )
#            a.cigar = ((8,len(sequence)))
#            a.mrnm = 0
#            a.mpos = 0
#            a.isize = 0
#            a.qual = quality
#            a.tags = (  )
            print line
#            outfile.write(a)
            continue
        
#    outfile.close()
        
run()