##################################
#                                #
# Last modified 01/28/2013       # 
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

    if len(sys.argv) < 3:
        print 'usage: python %s BAM samtools fields_to_be_retained' % sys.argv[0]
        print '\tfields to be retained format: comma separated, or from:to (including)'
        print '\tthe script will print to standard output - capture it with samtools and direct it to a bam file'
        print '\tspaces in optional fields will be treated as tabs'
        sys.exit(1)

    BAM = sys.argv[1]
    samtools = sys.argv[2]
    fieldIDs = []
    fields = sys.argv[3].split(',')
    for IDs in fields:
        if ':' in IDs:
            start = int(IDs.split(':')[0])
            end = int(IDs.split(':')[1]) + 1
            for i in range(start,end):
                fieldIDs.append(i)
        else:
            fieldIDs.append(int(IDs))
    fieldIDs.sort()

    cmd = samtools + ' view ' + BAM
    p = os.popen(cmd, "r")
    currentLine = p.readline()
    fields = currentLine.strip().split('\t')
    newfields = fields[0:12]
    for ID in range(12,len(fields)):
        newfields += fields[ID].split(' ')
    outline = ''
    for ID in fieldIDs:
        outline = outline + newfields[ID] + '\t'
    print outline.strip()
    line = currentLine
    while line != '':
        line = p.readline()
        if line == '':
            continue
        fields = line.strip().split('\t')
        if len(fields) == 11:
            print outline.strip()
            continue
        newfields = fields[0:12]
        for ID in range(12,len(fields)):
            newfields += fields[ID].split(' ')
        outline = ''
        for ID in fieldIDs:
            outline = outline + newfields[ID] + '\t'
        print outline.strip()
    
run()