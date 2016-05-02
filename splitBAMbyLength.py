##################################
#                                #
# Last modified 06/06/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set
import pysam
import string

def main(argv):

    if len(argv) < 4:
        print 'usage: python %s BAM chrom.sizes sizes outfileprefix [-nomulti]' % argv[0]
        print '\tsizes should be comma-separated and can be in from-to format as well as individual'
        sys.exit(1)

    BAM = argv[1]
    outfileprefix= argv[4]
    chrominfo=argv[2]
    chromInfoList=[]
    linelist=open(chrominfo)
    for line in linelist:
        fields=line.strip().split('\t')
        chr=fields[0]
        start=0
        end=int(fields[1])
        chromInfoList.append((chr,start,end))

    doMulti=True
    if '-nomulti' in argv:
        doMulti=False
        print 'will ignore multimappers'

    samfile = pysam.Samfile(BAM, "rb" )
    sizes = argv[3].split(',')
    sizeDict={}
    print sizes
    for size in sizes:
        if '-' in size:
            s1 = int(size.split('-')[0])
            s2 = int(size.split('-')[1])
            for s in range(s1,s2+1):
                sizeDict[s]= pysam.Samfile(outfileprefix + '.' + size + 'mers.bam', "wb", template=samfile)
        else:
            sizeDict[int(size)]= pysam.Samfile(outfileprefix + '.' + size + 'mers.bam', "wb", template=samfile)

    if doMulti:
        i=0
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                continue
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
                length = len(alignedread.seq)
                if sizeDict.has_key(length):
                    outfile = sizeDict[length]
                    outfile.write(alignedread)
    else:
        MultiplicityDict={}
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                continue
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                if MultiplicityDict.has_key(ID):
                    MultiplicityDict[ID]+=1
                else:
                    MultiplicityDict[ID]=1
        i=0
        for (chr,start,end) in chromInfoList:
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                continue
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print 'outputting alignments', str(i/1000000) + 'M alignments processed processed', chr,start,alignedread.pos,end
                fields=str(alignedread).split('\t')
                ID=fields[0]
                if alignedread.is_read1:
                    ID = ID + '/1'
                if alignedread.is_read2:
                    ID = ID + '/2'
                if nomulti and MultiplicityDict[ID] > 1:
                    continue
                else:
                    length = len(alignedread.seq)
                    if sizeDict.has_key(length):
                        outfile = sizeDict[length]
                        outfile.write(alignedread)

    outfilelist = []
    for s in sizeDict.keys():
        outfilelist.append(sizeDict[s])
    outfilelist = list(Set(outfilelist))
    for outfile in outfilelist:
        print outfile
        outfile.close()

if __name__ == '__main__':
    main(sys.argv)

