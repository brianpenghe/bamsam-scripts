##################################
#                                #
# Last modified 01/17/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
from commoncode import *

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s bedfilename SAMfilename outputfilename [-chr chrN1,chrN2...,chrNx]' % sys.argv[0]
        sys.exit(1)
    
    bed = sys.argv[1]
    SAM = sys.argv[2]
    outfilename = sys.argv[3]

    wantedDict={}

    doChr=False
    if '-chr' in sys.argv:
        doChr=True
        ChrList=sys.argv[sys.argv.index('-chr')+1].split(',')
        WantedChrDict={}
        for chr in ChrList:
            WantedChrDict[chr]=''
        print 'will output all alignments to', ChrList

    if doChr:
         pass
    else:
        lineslist = open(bed)
        for line in lineslist:
            if line[0]=='#':
                continue
            if line.startswith('track'):
                continue
            fields=line.strip().split('\t')
            chr=fields[0]
            start=int(fields[1])
            stop=int(fields[2])
            if wantedDict.has_key(chr):
                pass
            else:
                wantedDict[chr]={}
            for i in range(start,stop):
                wantedDict[chr][i]=''

    outfile = open(outfilename, 'w')

    lineslist = open(SAM)
    i=0
    TotalScore=0
    for line in lineslist:
        i+=1
        if i % 1000000 == 0:
            print i, 'lines processed'
        fields=line.strip().split('\t')        
        chr=fields[2]
        if doChr:
            if WantedChrDict.has_key(chr):
                outfile.write(line)
            continue            
        start=int(fields[3])
        if wantedDict.has_key(chr):
            if wantedDict[chr].has_key(start):
                outfile.write(line)
           
    outfile.close()
   
run()
