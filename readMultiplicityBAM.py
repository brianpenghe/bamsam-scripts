##################################
#                                #
# Last modified 09/14/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import pysam
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAMfilename samtools outputfilename' % sys.argv[0]
        print 'BAM file need not be indexed or sorted'
        sys.exit(1)

    SAM = sys.argv[1]
    samtools = sys.argv[2]
    outputfilename = sys.argv[3]

    ReadMultiplicityDict={}

    print 'examining read multiplicty'

    i=0
    cmd = samtools + ' view ' + SAM
    p = os.popen(cmd, "r")
    line = '...'
    while line != '':
        line = p.readline()
        if line == '':
            continue
        i+=1
        if i % 5000000 == 0:
            print str(i/1000000) + 'M alignments processed'
        fields = line.strip().split('\t')
        ID=fields[0]
        if fields[2] == '*':
            continue
        if ReadMultiplicityDict.has_key(ID):
            pass
        else:
            ReadMultiplicityDict[ID]=0
        ReadMultiplicityDict[ID]+=1

    IDList = ReadMultiplicityDict.keys()
    IDList.sort()

    outfile=open(outputfilename,'w')
    for ID in IDList:
        outfile.write(ID + '\t' + str(ReadMultiplicityDict[ID]) + '\n')
    outfile.close()

run()