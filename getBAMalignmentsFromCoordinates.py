##################################
#                                #
# Last modified 04/13/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import pysam
import sys
import string
import math

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s bedfilename BAMfilename outputfilename' % argv[0]
        sys.exit(1)
    
    bed = argv[1]
    SAM = argv[2]
    outfilename = argv[3]

    samfile = pysam.Samfile(SAM, "rb" )
    outfile = pysam.Samfile(outfilename, "wb", template=samfile)
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
        for alignedread in samfile.fetch(chr, start, stop):
            outfile.write(alignedread)

    outfile.close()
   
if __name__ == '__main__':
    main(sys.argv)
