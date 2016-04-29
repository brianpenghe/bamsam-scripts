##################################
#                                #
# Last modified 03/18/2014       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import os

# FLAG field meaning
# 0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair
# 0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
# 0x0004 4 the query sequence itself is unmapped
# 0x0008 8 the mate is unmapped 1
# 0x0010 16 strand of the query (0 for forward; 1 for reverse strand)
# 0x0020 32 strand of the mate 1
# 0x0040 64 the read is the first read in a pair 1,2
# 0x0080 128 the read is the second read in a pair 1,2
# 0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records)
# 0x0200 512 the read fails platform/vendor quality checks
# 0x0400 1024 the read is either a PCR duplicate or an optical duplicate

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024]

    FLAGList=[]

    MaxNumberList=[]
    for i in Numbers:
        if i <= FLAG:
            MaxNumberList.append(i)

    Residual=FLAG
    maxPos = len(MaxNumberList)-1

    while Residual > 0:
        if MaxNumberList[maxPos] <= Residual:
            Residual = Residual - MaxNumberList[maxPos]
            FLAGList.append(MaxNumberList[maxPos])
            maxPos-=1
        else:
            maxPos-=1
  
    return FLAGList

def main(argv):

    if len(argv) < 1:
        print 'usage: python %s 1 | 2' % argv[0]
        print '\tthe script will read from stdin (from samtools) and print to stdout'
        print '\t1 or 2 - read1 or read2'
        sys.exit(1)

    read1or2 = argv[1]

    input_stream = sys.stdin
    for line in input_stream:
        fields=line.strip().split('\t')
        FLAGfields = FLAG(int(fields[1]))
        if 128 in FLAGfields and read1or2 == '2':
            print line.strip()
        if 64 in FLAGfields and read1or2 == '1':
            print line.strip()
            
if __name__ == '__main__':
    main(sys.argv)
