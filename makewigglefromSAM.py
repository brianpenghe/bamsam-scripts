##################################
#                                #
# Last modified 12/20/2011       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s title SAMfilename outputfilename [-stranded + | -] [-nomulti] [-noNH] [-multiReadCorrection] [-RPM] [-notitle] [-bam chrom.sizes] [-singlebasepair] [-chr chrN1(,chrN2....)]' % sys.argv[0]
        print '       Note: the script is designed to work with constant length reads; if the read length varies, then the weighted average of the length will be used'
        print '       Note: only run it on BAM files with NH filed indicating multiread multiplicity present; the -noNH option does not work on BAM files'
        print '       Note: -multiReadCorrection option deprecated and correcting for multireads is the default now'
        print '       Note: do not use the -chr subset option together with the RPM option - the script will not consider all reads in this case'
        sys.exit(1)
    
    cachePages = 2000000

    doSingleBP=False
    if '-singlebasepair' in sys.argv:
        doSingleBP=True

    doTitle=True
    if '-notitle' in sys.argv:
        doTitle=False

    title = sys.argv[1]
    SAM = sys.argv[2]
    outfilename = sys.argv[3]

    doChrSubset=False
    if '-chr' in sys.argv:
        doChrSubset=True
        WantedChrDict={}
        for chr in sys.argv[sys.argv.index('-chr')+1].split(','):
            WantedChrDict[chr]=''

    doBAM=False
    if '-bam' in sys.argv:
        doBAM=True
        chrominfo=sys.argv[sys.argv.index('-bam')+1]
        chromInfoList=[]
        linelist=open(chrominfo)
        for line in linelist:
            fields=line.strip().split('\t')
            chr=fields[0]
            start=0
            end=int(fields[1])
            chromInfoList.append((chr,start,end))

    noMulti=False
    if '-nomulti' in sys.argv:
        noMulti=True

    doRPM=False
    if '-RPM' in sys.argv:
        doRPM=True

    multireads=0

    noNH=False
    if '-noNH' in sys.argv:
        noNH=True
        CorrectionDict={}
        CorrectionDict['paired']={}
        CorrectionDict['single']={}
        i=0
        if doBAM:
            samfile = pysam.Samfile(SAM, "rb" )
            for (chr,start,end) in chromInfoList:
                try:
                    for alignedread in samfile.fetch(chr, 0, 100):
                        a='b'
                except:
                    print 'region', chr,start,end, 'not found in bam file, skipping'
                    continue
                print chr,start,end
                for alignedread in samfile.fetch(chr, start, end):
                    i+=1
                    if i % 1000000 == 0:
                        print 'examining read multiplicity', str(i/1000000) + 'M alignments processed processed', len(CorrectionDict['paired'].keys()), 'paired multireads', len(CorrectionDict['single'].keys()), 'reads', multireads, 'multiread alignments'
                    fields=str(alignedread).split('\t')
                    ID=fields[0]
                    if alignedread.is_read1:
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        ID = ID + '/2'
                    if CorrectionDict['single'].has_key(ID):
                        multireads+=1
                        pass
                    else:
                        CorrectionDict['single'][ID]=0
                    CorrectionDict['single'][ID]+=1
        else:
            i=0
            linelist = open(SAM)
            for line in linelist:
                i+=1
                if i % 5000000 == 0:
                    print 'examining read multiplicity', str(i/1000000) + 'M alignments processed', len(CorrectionDict['paired'].keys()), 'paired multireads', len(CorrectionDict['single'].keys()), 'unpaired multireads'
                if line.startswith('@'):
                    continue
                fields=line.strip().split('\t')
                ID=fields[0]
                if fields[4]=='255':
                    continue
                if fields[7]!='0':
                    if CorrectionDict['paired'].has_key(ID):
                        pass
                    else:
                        CorrectionDict['paired'][ID]=0
                    CorrectionDict['paired'][ID]+=1
                else:
                    if CorrectionDict['single'].has_key(ID):
                        pass
                    else:
                        CorrectionDict['single'][ID]=0
                    CorrectionDict['single'][ID]+=1

    doStranded=False
    if '-stranded' in sys.argv:
         doStranded=True
         strand=sys.argv[sys.argv.index('-stranded')+1]
         print 'will only consider', strand, 'strand reads'

    coverageDict={}
    cumSumCoverage=0.0
    i=0
    r=0
    if doBAM:
        samfile = pysam.Samfile(SAM, "rb" )
        for (chr,start,end) in chromInfoList:
            if doChrSubset:
                if WantedChrDict.has_key(chr):
                    pass
                else:
                    continue
            try:
                for alignedread in samfile.fetch(chr, 0, 100):
                    a='b'
            except:
                print 'region', chr,start,end, 'not found in bam file, skipping'
                continue
            print chr,start,end
            for alignedread in samfile.fetch(chr, start, end):
                i+=1
                if i % 5000000 == 0:
                    print str(i/1000000) + 'M alignments processed'
                fields=str(alignedread).split('\t')
                if noNH:
                    if alignedread.is_read1:
                        ID = ID + '/1'
                    if alignedread.is_read2:
                        ID = ID + '/2'
                    scaleby = float(CorrectionDict['single'][ID])
                else:
                    try:
                        scaleby = alignedread.opt('NH')
                    except:
                        print 'multireads not specified with the NH tag, exiting'
                        sys.exit(1)
                if noMulti and scaleby > 1:
                    continue
                r+=1.0/scaleby
                if doStranded:
                    if 'XS' in fields[8]:
                        s = fields[8].split("'XS', '")[1].split("')")[0]
                        if s!=strand:
                            length=0
                            for (m,bp) in alignedread.cigar:
                                if m == 0:
                                    cumSumCoverage+=bp/scaleby
                            continue
                    else:
                        print 'BAM file does not have XS tags for all reads, exiting'
                        sys.exit(1)
                chr=samfile.getrname(alignedread.rname)
                if coverageDict.has_key(chr):
                    pass
                else:
                    coverageDict[chr]={} 
                start=alignedread.pos
                currentPos=start
                for (m,bp) in alignedread.cigar:
                    if m == 0:
                        cumSumCoverage+=bp/scaleby
                        for j in range(currentPos,currentPos+bp):
                            if coverageDict[chr].has_key(j+1):
                                coverageDict[chr][j+1]+=1/scaleby
                            else:
                                coverageDict[chr][j+1]=1/scaleby
                    elif m == 2:
                        pass
                    elif m == 3:
                        pass
                    else:
                        continue
                    currentPos=currentPos+bp
    else:
        linelist = open(SAM)
        for line in linelist:
            if line.startswith('@'):
                continue
            fields=line.strip().split('\t')
            chr=fields[2]
            if doChrSubset:
                if WantedChrDict.has_key(chr):
                    pass
                else:
                    continue
            start=int(fields[3])
            i+=1
            if i % 5000000 == 0:
                print str(i/1000000) + 'M alignments processed'
            if noNH:
                ID=fields[0]
                if CorrectionDict['single'].has_key(ID):
                    scaleby = float(CorrectionDict['single'][ID])
                elif CorrectionDict['paired'].has_key(ID):
                    scaleby = CorrectionDict['paired'][ID]/2.0
                else:
                    scaleby = 1.0
            else:
                if 'NH:i:' not in line:
                    print 'multireads not specified with the NH tag, exiting'
                    sys.exit(1)
                scaleby = float(line.split('NH:i:')[1].split('\t')[0].strip())
            if noMulti and scaleby > 1:
                continue
            r+=1.0/scaleby
            if doStranded:
                if 'XS:A:' not in line:
                    print 'SAM file does not have strand assignment that can be parsed by this script, '
                    sys.exit(1)
                s = line.split('XS:A:')[1].split('\t')[0].strip()
                if s != strand:
                    length=0
                    alignFields=fields[5].replace('D','N').split('N')
                    for a in alignFields:
                        length+=int(a.split('M')[0])
                    cumSumCoverage+=length/scaleby
                    continue
            if coverageDict.has_key(chr):
                pass
            else:
                coverageDict[chr]={} 
            alignFields=fields[5].replace('D','N').split('N')
            current=start
            for AF in alignFields:
                ExonMatch=int(AF.split('M')[0])
                stop=current+ExonMatch
                for j in range(current,stop):
                    cumSumCoverage+=1.0/scaleby
                    if coverageDict[chr].has_key(j):
                        coverageDict[chr][j]+=1.0/scaleby
                    else:
                        coverageDict[chr][j]=1.0/scaleby
                current=current+ExonMatch
                if AF.split('M')[1]!='':
                    intron=int(AF.split('M')[1])
                    current=current+intron

    readNumber=r

    if noNH:
        readNumber=i
        for ID in CorrectionDict['paired']:
            readNumber = readNumber - CorrectionDict['paired'][ID]/2.0 +1
        for ID in CorrectionDict['single']:
            readNumber = readNumber - CorrectionDict['single'][ID] + 1

    print 'read number', r
    print 'cumulative coverage', cumSumCoverage

    readLength=cumSumCoverage/readNumber

    print 'Estimated read length', readLength

    if doRPM:
        normFactor = readNumber/1000000.0

    outfile = open(outfilename, 'w')
    if doTitle:
        outline='track type=bedGraph name="' + title + '"'
        outfile.write(outline+'\n')

    chromosomes=coverageDict.keys()
    chromosomes.sort()
    for chr in chromosomes:
        print chr
        posKeys=coverageDict[chr].keys()
        posKeys.sort()
        initial=(posKeys[0],coverageDict[chr][posKeys[0]])
        previous=(posKeys[0],coverageDict[chr][posKeys[0]])
        written=['']
        if doSingleBP:
            for i in range(1,max(posKeys)+1):
                if coverageDict[chr].has_key(i):
                    outline = chr + '\t' + str(i) + '\t' + str(i+1) + str(coverageDict[chr][i])
                    outfile.write(outline+'\n')
                else:
                    outline = chr + '\t' + str(i) + '\t' + str(i+1) + str(0)
                    outfile.write(outline+'\n')
        else:
            for i in posKeys[1:len(posKeys)]:
                if previous[0]+1 == i and previous[1]==coverageDict[chr][i]:
                     previous=(i,coverageDict[chr][i])
                else:
                     if written[0]==initial[0]:
                         print written, initial, previous
                     if doStranded and strand == '-':
                         if doRPM:
                             outline=chr+'\t'+str(initial[0])+'\t'+str(previous[0]+1)+'\t-'+str(initial[1]/normFactor)
                         else:
                             outline=chr+'\t'+str(initial[0])+'\t'+str(previous[0]+1)+'\t-'+str(initial[1])
                     else:
                         if doRPM:
                             outline=chr+'\t'+str(initial[0])+'\t'+str(previous[0]+1)+'\t'+str(initial[1]/normFactor)
                         else:
                             outline=chr+'\t'+str(initial[0])+'\t'+str(previous[0]+1)+'\t'+str(initial[1])
                     written=(initial[0],previous[0]+1)
                     outfile.write(outline+'\n')
                     initial=(i,coverageDict[chr][i])
                     previous=(i,coverageDict[chr][i])

    outfile.close()
            
run()
