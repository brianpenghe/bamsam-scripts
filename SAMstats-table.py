##################################
#                                #
# Last modified 07/28/2010       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys

def main(argv):

    if len(argv) < 2:
        print 'usage: python %s <list of files> outputfilename' % argv[0]
        print '	Format of input file: <lable> <tab> <filename>'
        sys.exit(1)

    input = argv[1]
    outputfilename = argv[2]

    DataDict={}
    LabelList=[]

    lineslist = open(input)
    for line in lineslist:
        if line.strip() == '':
            continue
        fields = line.strip().split('\t')
        label=fields[0]
        file=fields[1]
        listoflines=open(file)
        LabelList.append(label)
        for line in listoflines:
            fields = line.strip().split(':\t')
            category=fields[0]
            if category == 'MultiSplices':
                category = 'Multi Splices'
            if category == 'UniqueSplices':
                category = 'Unique Splices'
            counts=fields[1]
            if DataDict.has_key(category):
                DataDict[category][label]=counts
            else:
                DataDict[category]={}
                DataDict[category][label]=counts

    outfile=open(outputfilename, 'w')

    CatKeys=DataDict.keys()
    CatKeys.sort()
    outline='#'
    for Category in CatKeys:
        outline=outline+'\t'+Category
    outfile.write(outline+'\n')
    LabelList.sort()
    for label in LabelList:
        outline=label
        for Category in CatKeys:
            outline=outline+'\t'+DataDict[Category][label]
        outfile.write(outline+'\n')
        
    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
