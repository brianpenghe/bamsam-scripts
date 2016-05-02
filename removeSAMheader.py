##################################
#                                #
# Last modified 04/24/2012       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys

try:
	import psyco
	psyco.full()
except:
	pass

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s input outfilename' % argv[0]

    input = argv[1]
    outputfilename = argv[2]

    outfile = open(outputfilename, 'w')

    listoflines = open(input)
    i=0
    for line in listoflines:
        i+=1
        if i % 5000000 == 0:
            print str(i/1000000) + 'M alignments processed'
        if line.startswith('@'):
            continue
        outfile.write(line)

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)

