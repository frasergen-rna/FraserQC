import sys
import os

if len(sys.argv) != 3:
	print "Usage: python %prog fq1 fq2"
	print "Function: check paired seq name for two fastq files"
	sys.exit(1)

lineNum = 0
seqNum = 0
with open(sys.argv[1]) as fp1, open(sys.argv[2]) as fp2:
    for line1 in fp1:
	lineNum += 1
        line2 = fp2.readline()

	if lineNum % 4 == 1:
		readName1 = line1.split()[0]
		readName2 = line2.split()[0] 
		if readName1 == readName2:
			seqNum += 1
		else:
			print "read name are not identical for reads", seqNum, " in the ", lineNum, " line of the two fastqc file."
			sys.exit(0)

print "identify ", seqNum, "reads were all identical from two fastq files"	

