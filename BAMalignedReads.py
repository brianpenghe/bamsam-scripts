##################################
#                                #
# Last modified 05/27/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import pysam
import string
from sets import Set
import os
import subprocess

def run():

    if len(sys.argv) < 2:
        print 'usage: python %s BAM samtools' % sys.argv[0]
        print '\tthe script will print to standard output - capture it with samtools and direct it to a bam file'
        print '\tspaces in optional fields will be treated as tabs'
        sys.exit(1)

    BAM = sys.argv[1]
    samtools = sys.argv[2]

    cmd = samtools + ' view -H ' + BAM
    p = os.popen(cmd, "r")
    currentLine = p.readline()
    line = currentLine
    print line.strip()
    while line != '':
        line = p.readline()
        if line == '':
            continue
        print line.strip()

    cmd = samtools + ' view ' + BAM
    p = os.popen(cmd, "r")
    currentLine = p.readline()
    fields = currentLine.strip().split('\t')
    aligned = True
    if fields[2] == '*':
        aligned = False
    line = currentLine
    if aligned:
        print line.strip()
    while line != '':
        line = p.readline()
        if line == '':
            continue
        fields = line.strip().split('\t')
        aligned = True
        if fields[2] == '*':
            aligned = False
        if aligned:
            print line.strip()
    
run()