##################################
#                                #
# Last modified 06/26/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import random
import pysam
import string
from sets import Set

try:
	import psyco
	psyco.full()
except:
	pass

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s SAMfilename outputfilename [-bam]' % sys.argv[0]
        sys.exit(1)

    SAM = sys.argv[1]
    outputfilename = sys.argv[2]

    doBAM = False
    if '-bam' in sys.argv:
        doBAM = True

    outfile=open(outputfilename, 'w')

    ChrDict={}

    if doBAM:
        samfile = pysam.Samfile(SAM, "rb" )
        for read in samfile.fetch(until_eof=True):
            chr = samfile.getrname(read.tid)
            if ChrDict.has_key(chr):
                continue
            else:
                ChrDict[chr]=''
                print chr
                outfile.write(chr + '\n')
    else:
        lineslist = open(SAM)
        for line in lineslist:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            chr = fields[2]
            if ChrDict.has_key(chr):
                continue
            else:
                ChrDict[chr]=''
                print chr
                outfile.write(chr + '\n')

    outfile.close()

run()
