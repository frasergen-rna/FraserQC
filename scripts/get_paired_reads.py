import sys
import os

if len(sys.argv) != 3:
	print "Usage: python %prog fq1 fq2"
	print "Function: get paired seq name for two fastq files"
	sys.exit(1)

lineNum = 0
readName1Dict = {}
fIn1 = open(sys.argv[1], 'r')
for line in fIn1:
	lineNum += 1
	if lineNum % 4 == 1:
		readName = line.split()[0]
		if readName[0] == "@":
			readName1Dict[readName] = ""
		else:
			print "wrong read name in the  ", lineNum, " of the fastq file: ", sys.argv[1]
			sys.exit(1)
fIn1.close()

lineNum = 0
readNameCommonDict = {}
fIn2 = open(sys.argv[2], 'r')
for line in fIn2:
	lineNum += 1
        if lineNum % 4 == 1:
                readName = line.split()[0]
                if readName[0] == "@" and  readName in readName1Dict:
                        readNameCommonDict[readsName] = ""
                else:
                        print "wrong read name in the  ", lineNum, " of the fastq file: ", sys.argv[2]
                        sys.exit(1)
fIn2.close()

flag = 1
lineNum = 0
fIn1 = open(sys.argv[1], 'r')
fOut1 = open(sys.argv[1]+"PE_1.fq", 'w')
for line in fIn1:
        lineNum += 1
        if lineNum % 4 == 1:
                readName = line.split()[0]
              	if readsName not in readNameCommonDict:
			flag = 0
		else:
			flag = 1	
	if flag == 1:
		print >> fOut1, line, 
fIn1.close()

flag = 1
lineNum = 0
fIn2 = open(sys.argv[2], 'r')
fOut2 = open(sys.argv[2]+"PE_2.fq", 'w')
for line in fIn2:
        lineNum += 1
        if lineNum % 4 == 1:
                readName = line.split()[0]
                if readsName not in readNameCommonDict:
                        flag = 0
                else:
                        flag = 1        
        if flag == 1:
                print >> fOut2, line, 
fIn2.close()
fOut2.close()
