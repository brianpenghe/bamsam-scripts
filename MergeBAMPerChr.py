##################################
#                                #
# Last modified 12/14/2013       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import math
import os
import subprocess

def run():

    if len(sys.argv) < 5:
        print 'usage: python %s samtools BAM1 chromosome1 BAM2 chromosome2'
        print '\tchromosomes should be comma-separated'
        print '\tthe script will print to stdout by default'
        sys.exit(1)

    samtools = sys.argv[1]
    BAM1 = sys.argv[2]
    chrDict1 = {}
    for chr in sys.argv[3].split(','):
        chrDict1[chr]=1
    BAM2 = sys.argv[4]
    chrDict2 = {}
    for chr in sys.argv[5].split(','):
        chrDict2[chr]=1

    cmd1 = samtools + ' view ' + BAM1
    p1 = os.popen(cmd1, "r")
    line = 'line1'
    while line != '':
        line = p1.readline().strip()
        if line == '':
            continue
        fields = line.strip().split('\t')
        chr = fields[2]
        if chrDict1.has_key(chr):
            print line.strip()

    cmd2 = samtools + ' view ' + BAM2
    p2 = os.popen(cmd2, "r")
    line = 'line2'
    while line != '':
        line = p2.readline().strip()
        if line == '':
            continue
        fields = line.strip().split('\t')
        chr = fields[2]
        if chrDict2.has_key(chr):
            print line.strip()
        
run()