#!/usr/bin/env python

'''
Written by Lauren Chircus
lchircus@stanford.edu
Stanford Univeristy
August 1, 2014

This script accepts a fastq file and returns a fastq file with the reverse complement of the reads in the original file and the reverse of the qScores in the original file.
'''

''' Import Modules '''
# import necessary for python
import os
import re
import sys
import bz2
import gzip
import string
import datetime
from argparse import ArgumentParser


''' Functions '''

# Reverse complement
complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

# Writes a single line with a newline to the specified file
# assumes that file is already open
def writeLine(myfile, line):
    myfile.write(line+'\n')

# Prints a line to the console and writes it to the specified file
def printLine(myfile, line):
    print line
    writeLine(myfile, line)


def main():
    #### OPTIONS ####
    # define options
    opts = ArgumentParser(description="This script accepts a fastq file and returns a fastq file with the reverse complement of the reads in the original file and the reverse of the qScores in the original file.")

    opts.add_argument("-i", dest="inFile", help="Input file; accepts fastq or fastq.gz", required=True)
    opts.add_argument("-u",  dest="u", action="store_true", default=False, help="Print uncompressed output file")
    opts.add_argument("-L", dest='lenTarget', help="length of the region of interest in read", action='store', type=int, required=True)
    opts.add_argument("-S", dest='startTarget', help="start position of the  region of interest in read", action='store', type=int, required=True)

    options = opts.parse_args()

    # return usage information if no argvs given
    if len(sys.argv)==1:
        os.system(sys.argv[0]+" --help")
        sys.exit()

    ##### INPUTS AND OUTPUTS #####
    # name input and outputs
    inputFile = options.inFile
    startTarget = options.startTarget -1
    endTarget = startTarget + options.lenTarget

    # name outputs and print to working dir
    fastq_file = inputFile.split('/')[-1]

    #check for file type and open input file
    fileType = inputFile.split('.')[-1]
    if fileType == "fastq":
        fastqFile = open(inputFile,'r')
        fastq_out = re.sub(".fastq", "_rc.fastq", fastq_file)
    elif fileType == "fq":
        fastqFile = open(inputFile,'r')
        fastq_out = re.sub(".fq", "_rc.fastq", fastq_file)
    elif fileType == "gz":
        fastqFile = gzip.open(inputFile,'r')
        fastq_out = re.sub(".fastq.gz", "_rc.fastq", fastq_file)
    elif fileType == "bz2":
        fastqFile = bz2.BZ2File(inputFile,'r')
        fastq_out = re.sub(".fastq.bz2", "_rc.fastq", fastq_file)
    else:
        sys.exit("ERROR! The input file must be a .fastq or .fastq.gz")

    ## Write to log
    logFileName = os.path.basename(__file__).rstrip(".py")+".log" # log for script
    logFile = open(logFileName,'w')
    printLine(logFile, 'Analysis performed ' + str(datetime.datetime.now())  + ' using ' + os.path.basename(__file__))
    printLine(logFile, 'Input fastq: ' + options.inFile)
    printLine(logFile, 'Output printed to: ' + fastq_out)
    printLine(logFile, 'Region of read reverse complemented began at position: ' + str(options.startTarget))
    printLine(logFile, 'Length of region reverse complemented: ' + str(options.lenTarget))

    # initialize variables
    nSequences=0
    count=1

    # initilize write files
    if options.u == False:
        outputFile = gzip.open(fastq_out+'.gz', 'wb')
    elif options.u == True:
        outputFile = open(fastq_out, 'w')

    while 1:
        # read line
        fastq_line = fastqFile.readline()

        # break if at end of file
        if not fastq_line:
            break

        # load fastq chunk
        if count ==1:
            sequenceHeader = fastq_line
        elif count ==2:
            seq = fastq_line.rstrip()
            seq = seq[startTarget:endTarget]
        elif count ==3:
            qualityScoreHeader = fastq_line
        elif count ==4:
            qual = fastq_line.rstrip()
            qual = qual[startTarget:endTarget]

            # reverse complement reads
            nSequences = nSequences+1  # total reads
            rc_seq = reverse_complement(seq)

            #reverse the quality scores
            qual = qual[::-1]

            # print data
            outputFile.write(sequenceHeader)
            outputFile.write(rc_seq+"\n")
            outputFile.write(qualityScoreHeader)
            outputFile.write(qual+"\n")

        # increment count
        count = count + 1
        if count == 5:
            count = 1
        else:
            count = count

    # close files to write the file
    outputFile.close()
    fastqFile.close()

    # give summary
    printLine(logFile, str(nSequences)+" total sequences reverse complemented")
    printLine(logFile, 'Analysis cpmpleted ' + str(datetime.datetime.now()))

if __name__ == '__main__':
    main()
