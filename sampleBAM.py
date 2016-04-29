##################################
#                                #
# Last modified 02/07/2013       # 
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

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s SAMfilename SAMTOOLS <Number Reads>' % argv[0]
        print '	Note1: do not run this code on SAM files where identical read IDs correspond to different reads'
        print '	      i.e. if for example reads come from lane1 of flowcell A and lane1 of flowcell B'
        print '	Note2: if reads are paired-end, the program will return the specified number of pairs, not individual reads'
        print '	Note3: the script will print alignments ot standard output - use samtools to capture them into a bam file'
        sys.exit(1)

    SAM = argv[1]
    samtools = argv[2]
    numReads = int(argv[3])

    doPaired=False
    if '-paired' in argv:
        doPaired=True 

    i=0
    ReadDict={}
    cmd = samtools + ' view ' + SAM
    p = os.popen(cmd, "r")
    line = p.readline()
    fields = line.strip().split('\t')
    if fields[0].endswith('/1') or fields[0].endswith('/2'):
        ID=fields[0][0:-2]
    else:
        ID=fields[0]
    ReadDict[ID]=0
    while line != '':
        line = p.readline()
        if line == '':
            continue
        fields = line.strip().split('\t')
        if fields[0].endswith('/1') or fields[0].endswith('/2'):
            ID=fields[0][0:-2]
        else:
            ID=fields[0]
        ReadDict[ID]=0

    ReadIDList=ReadDict.keys()

    SubSampleIDs=random.sample(ReadIDList,numReads)
    SubSampleIDDict={}
    for ID in SubSampleIDs:
        SubSampleIDDict[ID]=''

    ReadIDList=[]
    ReadDict={}
    SubSampleIDs=[]

    cmd = samtools + ' view -H ' + SAM
    p = os.popen(cmd, "r")
    line = p.readline()
    print line.strip()
    while line != '':
        line = p.readline()
        if line == '':
            continue
        print line.strip()

    cmd = samtools + ' view ' + SAM
    p = os.popen(cmd, "r")
    line = p.readline()
    fields = line.strip().split('\t')
    if fields[0].endswith('/1') or fields[0].endswith('/2'):
        ID=fields[0][0:-2]
    else:
        ID=fields[0]
    if SubSampleIDDict.has_key(ID):
        print line.strip()
    while line != '':
        line = p.readline()
        if line == '':
            continue
        fields = line.strip().split('\t')
        if fields[0].endswith('/1') or fields[0].endswith('/2'):
            ID=fields[0][0:-2]
        else:
            ID=fields[0]
        if SubSampleIDDict.has_key(ID):
            print line.strip()

if __name__ == '__main__':
    main(sys.argv)
