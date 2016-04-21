##################################
#                                #
# Last modified 02/27/2012       # 
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

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s SAM <max multiplicity> outfilename ' % sys.argv[0]
        print '       ONLY USE THIS SCRIPT FOR SINGLE-END DATA!!!!'
        sys.exit(1)

    inputfilename = sys.argv[1]
    M = int(sys.argv[2])
    outputfilename = sys.argv[3]

    ReadDict={}

    lineslist = open(inputfilename)
    i=0
    print 'examining read multiplicity'
    for line in lineslist:
        i+=1
        if i % 5000000 == 0:
            print str(i/1000000) + 'M alignments processed'
        fields = line.strip().split('\t')
        if line.startswith('@') and len(fields) < 10:
            continue
        if fields[2] == '*':
            continue
        ID = fields[0]
        if ReadDict.has_key(ID):
            ReadDict[ID]+=1
        else:
            ReadDict[ID]=1

    outfile = open(outputfilename, 'w')
    lineslist = open(inputfilename)
    i=0
    kept=0
    print 'outputting alignments'
    for line in lineslist:
        i+=1
        if i % 5000000 == 0:
            print str(i/1000000) + 'M alignments processed'
        fields = line.strip().split('\t')
        if line.startswith('@') and len(fields) < 10:
            outfile.write(line)
        if fields[2] == '*':
            continue
        ID = fields[0]
        if ReadDict.has_key(ID):
            if ReadDict[ID] <= M:
                outfile.write(line)
                kept+=1
            else:
                del ReadDict[ID] 

    print 'retained', kept, 'alignments'

    outfile.close()


run()

