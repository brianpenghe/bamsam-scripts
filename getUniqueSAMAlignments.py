##################################
#                                #
# Last modified 02/02/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
from sets import Set

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s SAMfilename outputfilename' % argv[0]
        sys.exit(1)

    SAM = argv[1]
    outputfilename = argv[2]

    outfile=open(outputfilename, 'w')

    i=0
    lineslist = open(SAM)
    for line in lineslist:
        if line.startswith('@'):
            continue
        i+=1
        if i % 1000000 == 0:
            print i, 'lines processed'
        fields = line.strip().split('\t')
        if fields[4]=='255':
            outfile.write(line)

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
