##################################
#                                #
# Last modified 06/17/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import pysam
import string

def run():

    if len(sys.argv) < 4:
        print 'usage: python %s BAM samtools max_multiplcity outfilename' % sys.argv[0]
        sys.exit(1)

    BAM = sys.argv[1]
    samtools=sys.argv[2]
    M = int(sys.argv[3])
    outputfilename = sys.argv[4]
#    chromInfoList=[]
#    linelist=open(chrominfo)
#    for line in linelist:
#        fields=line.strip().split('\t')
#        chr=fields[0]
#        start=0
#        end=int(fields[1])
#        chromInfoList.append((chr,start,end))

#    doIgnoreChromSizes=False
#    if '-ignoreChromSizes' in sys.argv:
#        doIgnoreChromSizes=True
#        print 'will ignore chrom.sizes file'

    print 'will write alignments into: ', outputfilename

    samfile = pysam.Samfile(BAM, "rb" )
    try:
        for alignedread in samfile.fetch():
            chr = samfile.getrname(alignedread.tid)
            print chr, alignedread
            if chr == '*':
                continue
            else:
                multiplicity = alignedread.opt('NH')
                print 'file has NH tags'
                break
    except:
        print 'no NH: tags in BAM file, will replace with a new BAM file with NH tags'
        BAMpreporcessingScript = sys.argv[0].rpartition('/')[0] + '/bamPreprocessing.py'
        cmd = 'python ' + BAMpreporcessingScript + ' ' + BAM + ' ' + BAM + '.NH'
        os.system(cmd)
        cmd = 'rm ' + BAM
        os.system(cmd)
        cmd = 'mv ' + BAM + '.NH' + ' ' + BAM
        os.system(cmd)
        cmd = samtools + ' index ' + BAM
        os.system(cmd)

#    MultiplicityDict={}
#    i=0
    samfile = pysam.Samfile(BAM, "rb" )
#    for alignedread in samfile.fetch():
#        i+=1
#        if i % 5000000 == 0:
#            print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed'
#        fields=str(alignedread).split('\t')
#        ID=fields[0]
#        if alignedread.is_read1:
#            ID = ID + '/1'
#        if alignedread.is_read2:
#            ID = ID + '/2'
#        if MultiplicityDict.has_key(ID):
#            MultiplicityDict[ID]+=1
#        else:
#            MultiplicityDict[ID]=1
    outfile = pysam.Samfile(outputfilename, "wb", template=samfile)
    o=0
    i=0
    for alignedread in samfile.fetch():
        i+=1
        if i % 5000000 == 0:
            print 'outputting alignments', str(i/1000000) + 'M alignments processed processed'
#        fields=str(alignedread).split('\t')
#        ID=fields[0]
#        if alignedread.is_read1:
#            ID = ID + '/1'
#        if alignedread.is_read2:
#            ID = ID + '/2'
#        if MultiplicityDict[ID] != alignedread.opt('NH'):
#            print alignedread
#            print MultiplicityDict[ID], alignedread.opt('NH')
        if alignedread.opt('NH') <= M:
            o+=1
            outfile.write(alignedread)
#        if MultiplicityDict[ID] == 1:
#            del MultiplicityDict[ID]
#    else:
#        MultiplicityDict={}
#        i=0
#        samfile = pysam.Samfile(BAM, "rb" )
#        for (chr,start,end) in chromInfoList:
#            try:
#                for alignedread in samfile.fetch(chr, 0, 100):
#                    a='b'
#            except:
#                continue
#            for alignedread in samfile.fetch(chr, start, end):
#                i+=1
#                if i % 5000000 == 0:
#                    print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
#                fields=str(alignedread).split('\t')
#                ID=fields[0]
#                if alignedread.is_read1:
#                    ID = ID + '/1'
#                if alignedread.is_read2:
#                    ID = ID + '/2'
#                if MultiplicityDict.has_key(ID):
#                    MultiplicityDict[ID]+=1
#                else:
#                    MultiplicityDict[ID]=1
#        outfile = pysam.Samfile(outputfilename, "wb", template=samfile)
#        o=0
#        i=0
#        for (chr,start,end) in chromInfoList:
#            try:
#                for alignedread in samfile.fetch(chr, 0, 100):
#                    a='b'
#            except:
#                continue
#            for alignedread in samfile.fetch(chr, start, end):
#                i+=1
#                if i % 5000000 == 0:
#                    print 'outputting alignments', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
#                fields=str(alignedread).split('\t')
#                ID=fields[0]
#                if alignedread.is_read1:
#                    ID = ID + '/1'
#                if alignedread.is_read2:
#                    ID = ID + '/2'
#                if MultiplicityDict[ID] != alignedread.opt('NH'):
#                    print alignedread
#                    print MultiplicityDict[ID], alignedread.opt('NH')
#                if alignedread.opt('NH') <= M:
#                    o+=1
#                    outfile.write(alignedread)
#                if MultiplicityDict[ID] == 1:
#                    del MultiplicityDict[ID]

#    MultiplicityDict={}

    print 'outputted', o, 'alignments'

    outfile.close()

run()

