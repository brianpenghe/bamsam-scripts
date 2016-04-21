##################################
#                                #
# Last modified 02/03/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import pysam
import random
import string
from sets import Set

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s SAMfilename <Number Reads> outputfilename' % sys.argv[0]
        print '	Note1: do not run this code on SAM files where identical read IDs correspond to different reads'
        print '	      i.e. if for example reads come from lane1 of flowcell A and lane1 of flowcell B'
        print '	Note2: if reads are paired-end, the program will return the specified number of pairs, not individual reads'
        sys.exit(1)

    SAM = sys.argv[1]
    numReads = int(sys.argv[2])
    outputfilename = sys.argv[3]

    doPaired=False
    if '-paired' in sys.argv:
        doPaired=True 

    i=0
    print 'Getting Read IDs'
    ReadIDList=[]
    lineslist = open(SAM)
    for line in lineslist:
        if line.startswith('@'):
            continue
        i+=1
        if i % 10000000 == 0:
            print str(i/1000000) + 'M alignments'
        fields = line.strip().split('\t')
        if fields[0].endswith('/1') or fields[0].endswith('/2'):
            ID=fields[0][0:-2]
        else:
            ID=fields[0]
        ReadIDList.append(ID)

    outfile=open(outputfilename, 'w')

    ReadIDList=list(Set(ReadIDList))

    print 'found', len(ReadIDList), 'fragments'

    SubSampleIDs=random.sample(ReadIDList,numReads)
    SubSampleIDDict={}
    for ID in SubSampleIDs:
        SubSampleIDDict[ID]=''

    print 'sampled', (len(SubSampleIDs)), 'IDs'

    ReadIDList=[]
    SubSampleIDs=[]

    print 'Outputting sub-sampled alignments'
    lineslist = open(SAM)
    i=0
    for line in lineslist:
        if line.startswith('@'):
#            outfile.write(line)
            continue
        i+=1
        if i % 10000000 == 0:
            print str(i/1000000) + 'M alignments'
        fields = line.strip().split('\t')
        if fields[0].endswith('/1') or fields[0].endswith('/2'):
            ID=fields[0][0:-2]
        else:
            ID=fields[0]
        if SubSampleIDDict.has_key(ID):
            outfile.write(line)

    print 'outputted', len(SubSampleIDDict.keys()), 'out of', len(ReadIDList), 'fragments'
 
    outfile.close()

run()
