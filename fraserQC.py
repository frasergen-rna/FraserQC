#!/usr/bin/python
#-*-coding:utf-8-*- 

#fraserQC for fastq data quality control and trim
#frasergen reserved all right regrading to the pipelines

import sys
import os
import re

from optparse import OptionParser
from optparse import OptionGroup

import multiprocessing

def rm_folder(folder):
	if os.path.exists(folder):
		os.system("rm -r "+folder)
	else:
		print "Warning: cannot find ", folder

def rm_file(filE):
        if os.path.isfile(filE):
                os.system("rm -r "+filE)

def create_folder(folder):
	if not os.path.exists(folder):
		os.system("mkdir "+folder)

##access fastqc for sample
def fastqc_run(sampleName, fileList, outputFolder, threads):
	#clean output files
	rm_folder(sampleName+"_fastqc")
	rm_folder(sampleName+"_htqc")

	fileType = fileList[0].split('.')[-1]
	statDict = {"readNum":" ", "baseNum":" ", "readLen":" ", "Q20Rate":" ", "Q30Rate":" ", "GCRate":" ", "passOverRepTest":" "}

	print "############################\nperform quality control for "+ sampleName + " ..."
	create_folder(outputFolder+"/"+sampleName)
	create_folder(sampleName+"_fastqc")
		
	#############################################################pick the first 100000 lines for fastqc
	#rm_file(sampleName+"_fastqc/tmp.fq")

	#pickNum = 4000000
	#annotated for single output
	#for fq in fileList:
	#	if fq.split('.')[-1] == "gz":
	#		cmd = "gunzip -c "+fq+" | head -n " + str(pickNum) +" >> "+sampleName+"_fastqc/tmp.fq"
	#	else:
	#		cmd = "head -n " + str(pickNum) + " " + fq + " >> " + sampleName+"_fastqc/tmp.fq"
	#	os.system(cmd)

	#now use paired output for paired-end reads
	#for i in xrange(1, len(fileList)+1):
	#	fq = fileList[i-1]
	#	if fq.split('.')[-1] == "gz":
        #                cmd = "gunzip -c "+fq+" | head -n " + str(pickNum) +" > "+sampleName+"_fastqc/tmp"+str(i)+".fq"
	#		os.system(cmd)
	#		cmd = progDict["fastqc"] + " -t " + str(threads) + " -q " + sampleName+"_fastqc/tmp"+str(i)+".fq -o "+sampleName+"_fastqc"
	#		os.system(cmd)
        #       else:
        #                cmd = "head -n " + str(pickNum) + " " + fq + " > " + sampleName+"_fastqc/tmp"+str(i)+".fq"
        #        	os.system(cmd)
	#		cmd = progDict["fastqc"] + " -t " + str(threads) + " -q " + sampleName+"_fastqc/tmp"+str(i)+".fq -o "+sampleName+"_fastqc"
	#		os.system(cmd)

	#fastqc run
	for fqFile in fileList:
		cmd = progDict["fastqc"] + " -t " + str(threads) + " -q " + fqFile +" -o "+sampleName+"_fastqc"
		#cmd = progDict["fastqc"] + " -t " + str(threads) + " " + fqFile +" -o "+sampleName+"_fastqc"	
		print "Performing command: " + cmd
		os.system(cmd)

	#collect fastqc result
	for i in xrange(1, len(fileList)+1):
		#############################################################extract fastqc output
		#rm_folder(sampleName+"_fastqc/"+sampleName+"_fastqc")
		###modified by liuzhenhua 20180529
		if "fastq.gz" in fileList[i-1]:
			zipFile = fileList[i-1].split('/')[-1].split('.fastq')[0]
		elif "fastq" in fileList[i-1]:
			zipFile = fileList[i-1].split('/')[-1].split('.fastq')[0]
		else:
			zipFile = fileList[i-1].split('/')[-1].split('.fq')[0]
		#print zipFile
		cmd = "unzip -q "+sampleName+"_fastqc/"+zipFile+"_fastqc.zip -d "+sampleName+"_fastqc/"
		print "Performing command: " + cmd
		os.system(cmd)

		#############################################################parse fastqc output
		#load in fastqc result
		FastqcOut = open(sampleName+"_fastqc/"+zipFile+"_fastqc/fastqc_data.txt", 'r')

		#output important info for plot
		flag = 0
		fOut = ""
		for line in FastqcOut:
			if '_'.join(line.split()[0:2]) == ">>Overrepresented_sequences":
				if statDict["passOverRepTest"] == " ":
					if line.replace('\n','').split()[2] == "pass":
						statDict["passOverRepTest"] = "Y"
					else:
						statDict["passOverRepTest"] = "N"
				else:
					if line.replace('\n','').split()[2] == "pass":
						statDict["passOverRepTest"] = statDict["passOverRepTest"]+";Y"
					else:
						statDict["passOverRepTest"] = statDict["passOverRepTest"]+";N"
		
			if '_'.join(line.split()[0:4]) == ">>Per_base_sequence_quality":
				flag = 1
				fOut = open(outputFolder+"/"+sampleName+"/per_base_sequence_quality_"+str(i), 'w')
				#print "open", fOut
				continue
	
			if '_'.join(line.split()[0:4]) == ">>Per_base_sequence_content":
				flag = 1
       	 			fOut = open(outputFolder+"/"+sampleName+"/per_base_sequence_content_"+str(i), 'w')
				#print "open", fOut
       		 		continue
	
			if '_'.join(line.split()[0:4]) == ">>Per_sequence_GC_content":
       			 	flag = 1
       		 		fOut = open(outputFolder+"/"+sampleName+"/per_sequence_GC_content_"+str(i), 'w')
				#print "open", fOut
       	 			continue

			if '_'.join(line.split()[0:3]) == ">>Sequence_Length_Distribution":
        			flag = 1
        			fOut = open(outputFolder+"/"+sampleName+"/sequence_length_distribution_"+str(i), 'w')
				#print "open", fOut
        			continue

			elif line.split()[0] == ">>END_MODULE":
				if fOut:
					#print "close", fOut
					fOut.close()
				flag = 0

			if flag == 1:
				print >> fOut, line,
	
		FastqcOut.close()
	
		cmd="cp "+sampleName+"_fastqc/"+zipFile+"_fastqc/fastqc_data.txt "+outputFolder+"/"+sampleName+"/fastqc_data_"+str(i)+".txt"
		print "Performing command: " + cmd
		os.system(cmd)

        #htqc run
        if len(fileList) == 1:
                if fileType == "gz":
                        cmd = progDict["htStat"] + " -q -z -i " +  fileList[0] + " -S -o " + sampleName + "_htqc"
                else:
                        cmd = progDict["htStat"] + " -q -i " +  fileList[0] + " -S -o " + sampleName + "_htqc"
        else:
                if fileType == "gz":
                        cmd = progDict["htStat"] + " -q -z -i " +  ' '.join(fileList) + " -P -o " + sampleName + "_htqc"
                else:
                        cmd = progDict["htStat"] + " -q -i " +  ' '.join(fileList) + " -P -o " + sampleName + "_htqc"
		print "Performing htstat ..."
		print "Performing command: " + cmd
        os.system(cmd)

	#collect htqc results
	fIn = open(sampleName + "_htqc/info.tab", 'r')
	for line in fIn:
		tmp = line.replace('\n', '').split()
                if tmp[0] == "amount":
              		statDict["readNum"] = tmp[1]
                	break
        fIn.close()
	
	readLen = ""
        baseNum = ""
        Q20Rate = ""
        Q30Rate = ""
	GCRate = ""

	for i in xrange(1, len(fileList)+1):
		all = 0
		Q20 = 0
		Q30 = 0
		GCNum = 0

		fIn = open(sampleName + "_htqc/cycle_quality_"+str(i)+".tab", 'r')
		lineNum = 0
		for line in fIn:
			lineNum += 1
			if line[0] != "#" and lineNum > 1:
				tmp = [int(elem) for elem in line.replace('\n', '').split()[1:]]
				all += sum(tmp)
				Q20 += sum(tmp[20:])
				Q30 += sum(tmp[30:])
		fIn.close()
		
		fIn = open(sampleName + "_htqc/cycle_composition_"+str(i)+".tab", 'r')
		lineNum = 0
		for line in fIn:
			lineNum += 1
                        if line[0] != "#" and lineNum > 1:
				tmp = [int(elem) for elem in line.replace('\n', '').split()[1:]]
				GCNum += sum(tmp[2:4])
		fIn.close()
		#print "GCNum=", GCNum
	
		if baseNum == "":
                        readLen = str(all/int(statDict["readNum"].split('x')[0]))
                        baseNum = str(all)
                        Q20Rate = str(round(float(Q20)/all, 3)*100)
                        Q30Rate = str(round(float(Q30)/all, 3)*100)
			GCRate = str(round(float(GCNum)/all, 3)*100)
                else:
                        readLen += ";"+str(all/int(statDict["readNum"].split('x')[0]))
                        baseNum += ";"+str(all)
                        Q20Rate += ";"+str(round(float(Q20)/all, 3)*100)
                        Q30Rate += ";"+str(round(float(Q30)/all, 3)*100)		
			GCRate += ";"+str(round(float(GCNum)/all, 3)*100)
		
	statDict["baseNum"] = baseNum
	statDict["readLen"] = readLen
	statDict["Q20Rate"] = Q20Rate
	statDict["Q30Rate"] = Q30Rate
	statDict["GCRate"] = GCRate

	#output quality results
	fOut = open(outputFolder+"/"+sampleName+"/quality_statistics.dat", 'w')
	print >> fOut, "SampleName\tReadNum\tBaseCount(bp)\tReadLength(bp)\tQ20(%)\tQ30(%)\tGC Content(%)\tpassOverRepTest"
	print >> fOut, sampleName+"\t"+statDict["readNum"]+"\t"+statDict["baseNum"]+"\t"+statDict["readLen"]+"\t"+statDict["Q20Rate"]+"\t"+statDict["Q30Rate"]+"\t"+statDict["GCRate"]+"\t"+statDict["passOverRepTest"]
	fOut.close()	

	#rm_folder(sampleName+"_fastqc")
	#rm_folder(sampleName+"_htqc")

##plot for samples from fastqc result
def qc_plot(sampleName, fileList, outputFolder, water):
	for i in xrange(1, len(fileList)+1):
		if os.path.isfile(outputFolder+"/"+sampleName+"/per_base_sequence_quality_"+str(i)):
			cmd = progDict["Rscript"] + " " + progDict["scriptPath"] + "/base_quality_box_plot.R " + outputFolder+"/"+sampleName+"/per_base_sequence_quality_"+str(i)+" "+str(water)
			os.system(cmd)
	
		if os.path.isfile(outputFolder+"/"+sampleName+"/per_base_sequence_content_"+str(i)):
       		        cmd = progDict["Rscript"] + " " + progDict["scriptPath"] + "/base_content_line_plot.R " + outputFolder+"/"+sampleName+"/per_base_sequence_content_"+str(i)+" "+str(water)
                	os.system(cmd)
		
		if os.path.isfile(outputFolder+"/"+sampleName+"/per_sequence_GC_content_"+str(i)):
			cmd = progDict["Rscript"] + " " + progDict["scriptPath"] + "/sequence_content_distribution_plot.R " + outputFolder+"/"+sampleName+"/per_sequence_GC_content_"+str(i)+" "+str(water)
			os.system(cmd)

		if os.path.isfile(outputFolder+"/"+sampleName+"/sequence_length_distribution_"+str(i)):
			cmd = progDict["Rscript"] + " " + progDict["scriptPath"] + "/sequence_length_distribution_plot.R "+outputFolder+"/"+sampleName+"/sequence_length_distribution_"+str(i)+" "+str(water)
			os.system(cmd)

#access htqc for samples
def htqc_run(sampleName, fileList, htqc_W, htqc_C, htqc_L):
	print "#################################\nquality trim for sample ", sampleName, " ..."
	#print fileList

	#check the data from that sample is PE or SE
	PE = 1
	if len(fileList) == 1:
		PE = 0
	elif len(fileList) == 2:
		PE = 1
	else:
		print "no or more than two fq files for sample "+sampleName			
		sys.exit(1)

	#quality trim for data according to PE or SE
	fileType = fileList[0].split('.')[-1]
	print "PE = ", PE
	print "File list = ", fileList
	if PE == 1:
		fastq1 = fileList[0]
		fastq2 = fileList[1]
		if fastq1.split('.')[-1] != fastq2.split('.')[-1]: 
			print "unconsistent file type for two reads", fastq1, fastq2
		else:
			if fileType == "gz":
				cmd = progDict["htTrim"] + " -q -i " + fastq1 + " -z -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz"
				#print "Performing command 1: " + cmd
				os.system(cmd)

				cmd = progDict["htTrim"] + " -q -i " + fastq2 + " -z -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq2.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz"
				#print "Performing command 1: " + cmd
				os.system(cmd)

				cmd = progDict["htFilter"] + " -q -i " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz " + fastq2.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz -z -P -F quality " + str(htqc_C) + " -L " + str(htqc_L) + " -o temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""
				os.system(cmd)
				#print "Performing command 1: " + cmd
				cmd = progDict["htFilter"] + " -q -i temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq.gz temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq.gz -z -P -F length -L " + str(htqc_L) + " -o " + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""
				#print "Performing command 1: " + cmd
				os.system(cmd)

				cmd = "rm " + fastq1.split('/')[-1]+"_Q" + str(htqc_C) + ".fastq.gz; rm " + fastq2.split('/')[-1]+"_Q" + str(htqc_C) + ".fastq.gz"
				os.system(cmd)	
				cmd = "rm temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq.gz temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq.gz"
				os.system(cmd)	
			
			elif fileType == "fastq" or fileType == "fq":
				cmd = progDict["htTrim"] + " -q -i " + fastq1 + " -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq"
				#print "Performing command 2: " + cmd #print cmd
				os.system(cmd)

				cmd = progDict["htTrim"] + " -q -i " + fastq2 + " -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq2.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq"
				#print "Performing command 2: " + cmd #print cmd
				os.system(cmd)

				cmd = progDict["htFilter"] + " -q -i " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq " + fastq2.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq -P -F length -L " + str(htqc_L) + " -o temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""     
				#print "Performing command 2: " + cmd #print cmd
				os.system(cmd)
				cmd = progDict["htFilter"] + " -q -i temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq -z -P -F length -L " + str(htqc_L) + " -o " + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""

				cmd = "rm " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq; rm " + fastq2.split('/')[-1] +"_Q" + str(htqc_C) + ".fastq"			
				os.system(cmd)
				cmd = "rm temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq temp" + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq"
				os.system(cmd)

			else:
				print "unrecognized format of " + fastq1.split('.')[-1] + " for raw data for sample " + sampleName
	else:
		fastq1 = fileList[0]
		if fileType == "gz":
			cmd = progDict["htTrim"] + " -q -i " + fastq1 + " -z -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz"
			#print "Performing command 3: " + cmd
			os.system(cmd)

			cmd = progDict["htFilter"] + " -q -i " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz -z -S -F length -L " + str(htqc_L) + " -o " + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""
			#print "Performing command 3: " + cmd
			os.system(cmd)

			cmd = "rm " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq.gz"
			#print "Performing command 3: " + cmd
			os.system(cmd)
		
		elif fileType == "fastq" or fileType == "fq":
			cmd = progDict["htTrim"] + " -q -i " + fastq1 + " -S both -C " + str(htqc_C) + " -W " + str(htqc_W) + " -o " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq"
			#print "Performing command 4: " + cmd
			os.system(cmd)

			cmd = progDict["htFilter"] + " -q -i " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq -S -F length -L " + str(htqc_L) + " -o " + sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ""     
			#print "Performing command 4: " + cmd
			os.system(cmd)
			cmd = "rm " + fastq1.split('/')[-1] + "_Q" + str(htqc_C) + ".fastq"
			#print "Performing command 4: " + cmd
			os.system(cmd)

		else:
			print "unrecognized format of " + fastq1.split('.')[-1] + " for raw data for sample " + sampleName


#check all required arguments were supplied	
def checkRequiredArguments(opts, parser):
	missing_options = []
	for option in parser.option_list:
        	if re.match(r'^\[REQUIRED\]', option.help) and eval('opts.' + option.dest) == None:
			missing_options.extend(option._long_opts)
    	if len(missing_options) > 0:
        	parser.error('Missing REQUIRED parameters: ' + str(missing_options))


##################parameters and help information
usage = '\n python %prog [options] -i fastqFile -s sampleName -o outputFolder'

parser = OptionParser(usage, version='%prog 1.0 by shijun xiao at 20170413')

##parameters for sampleName
parser.add_option('-s','--sampleName',
                  action='store',dest='sampleName',
                  help='[REQUIRED] sample name for fastq files')

parser.add_option('-i','--fqList',
                  action='store',dest='fqList',
                  help='[REQUIRED] input fastq files ended with fq/fastq/gz (PEs were separated by ",")')

##parameters for output folder
parser.add_option('-o','--outputFolder',
                  action='store',dest='outputFolder',
                  metavar='FOLDER',help='[REQUIRED] In this folder, the program will create sample-specific folder, containning all qc result')

##paramete to control quality trim or not
parser.add_option("-f", '--filter', 
		  action="store_true", dest="filter", 
		  help='perform quality trim by htqc [default to skip this step]')

parser.add_option("-t", '--threads',
                  action="store", dest="threads", default=1,
                  help='number of threads running [default to %default]')

parser.add_option("-m", '--watermarker',
                  action="store_false", dest="watermarker",
                  help='print out water marker [default to not print water marker]')

##parameters for htqc
group = OptionGroup(parser,'Advanced Options',
                    'Parameters for htqc')
group.add_option('-w',action='store', metavar='INT', dest='htqc_W', default=5, help='window size for quality trim [default: %default]')
group.add_option('-C',action='store', metavar='INT', dest='htqc_C', default=20, help='quality threshold for trim [default: %default]')
group.add_option('-L',action='store', metavar='INT', dest='htqc_L', default=50, help='length threshold for trim [default: %default]')

parser.add_option_group(group)
(options,args) = parser.parse_args()

if len(sys.argv) == 1:
	print
	parser.print_help()
	print
	sys.exit(0)
else:
	checkRequiredArguments(options, parser)
	

################ modify this part for program path
progDict={
"scriptPath": os.path.join(sys.path[0], "scripts"),
"fastqc": "/local_data1/pipeline/Transcriptome/PB_Isoseq_noref_V1.0/software/FastQC/fastqc",
"htTrim": "/local_data1/pipeline/Transcriptome/PB_Isoseq_noref_V1.0/software/htqc-0.90.8-Source/bin/ht-trim",
"htFilter": "/local_data1/pipeline/Transcriptome/PB_Isoseq_noref_V1.0/software/htqc-0.90.8-Source/bin/ht-filter",
"htStat": "/local_data1/pipeline/Transcriptome/PB_Isoseq_noref_V1.0/software/htqc-0.90.8-Source/bin/ht-stat",
"Rscript": "/local_data1/pipeline/Transcriptome/PB_Isoseq_noref_V1.0/software/anaconda2/bin/Rscript",
}

################check programs are executable
print "#############################\nchecking programs ..."
for prog in progDict:
        if not os.access(progDict[prog],os.X_OK):
                print progDict[prog], "for", prog, "was not executable"
                sys.exit(1)

#######################main function of the pipeline
sampleName = options.sampleName
fqList = options.fqList
outputFolder = options.outputFolder
htqc_W = options.htqc_W
htqc_C = options.htqc_C
htqc_L = options.htqc_L
water = 0

#waterproof or not
if options.watermarker:
	water = 1

#current absolute path
CWD = os.getcwd()

#create output folder
create_folder(outputFolder)

#check if all fastq files were there
print "#############################\nchecking fastq data ..."
fileList = fqList.split(",")
print "Input File List: ", fileList
for fq in fileList:
	if not os.path.isfile(fq):
		print "Cannot find fastq file", fq
		sys.exit(1)

#perform quality control from raw data
if len(fileList) == 2:
	if fileList[0].split('.')[-1] != fileList[1].split('.')[-1]:
		print "file type are not consistent for ", fileList[0], fileList[1]
		sys.exit(1)

fastqc_run(sampleName, fileList, outputFolder, options.threads)	
qc_plot(sampleName, fileList, outputFolder, water)
#print "File list after fqstqc_run and qc_plot", fileList
#perform quality trim and quality control again if needed
if options.filter:
	rm_flag = 0

	#load in sample info
	fileList = fqList.split(",")
	#print "File List for htqc_run: ", fileList
	#print "Threads = ===", options.threads, "^^^"
	#non-parallel mode
	if options.threads == 1:
		print "Running ht_qc with single process"
		htqc_run(sampleName, fileList, htqc_W, htqc_C, htqc_L)
	else:
		print "number of threads should be large than zero."
		sys.exit(1)
	#parallel mode
	"""
	elif options.threads > 1:
		print "Running ht_qc with multiple process"
		#if the file is compressed, extract fastq file and split into small parts
		rm_flag = 0
		if fileList[0].split('.')[-1] == "gz":
			rm_flag = 1
                	if len(fileList) == 2:
                	        cmd = "gunzip -c " + fileList[0] + " > " + fileList[0].split('/')[-1].strip(".gz") + ".fastq"
                       	 	#print cmd
                        	os.system(cmd)
                        	cmd = "gunzip -c " + fileList[1] + " > " + fileList[1].split('/')[-1].strip(".gz") + ".fastq"
                        	#print cmd
                        	os.system(cmd)
                	else:
                        	cmd = "gunzip -c " + fileList[0] + " > " + fileList[0].split('/')[-1].strip(".gz") + ".fastq"
                        	os.system(cmd)

                	fileList = [elem.split('/')[-1].strip(".gz")+".fastq" for elem in fileList]
		
		#split fastq files into parts
		cmd = "perl " + progDict["scriptPath"] + "/fastq-splitter.pl -n " + str(options.threads) + " " + ' '.join(fileList)
		#print cmd
		os.system(cmd)

		#deleted extracted fastq files if the original files were compressed
		if rm_flag == 1:
			for fqFile in fileList:
				rm_file(fqFile)

		#multi-processes to htqc
		pool = multiprocessing.Pool(processes=int(options.threads))
		for i in xrange(1, int(options.threads)+1):
			pool.apply_async(htqc_run, (sampleName+'_'+str(i), [elem.split('/')[-1]+".part-"+str(i)+".fastq" for elem in fileList], htqc_W, htqc_C, htqc_L, ))	
			
		pool.close()
		pool.join()

		#delete splite fastq parts
		for i in xrange(1, int(options.threads)+1):
			for elem in fileList:
				cmd = "rm " + elem.split('/')[-1] + ".part-" + str(i) + ".fastq"
				os.system(cmd)

		#merge htqc results
		if len(fileList) == 2:
			#rm_file(sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq")
			#rm_file(sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq")
			for i in xrange(1, int(options.threads)+1):
				cmd = "cat " + sampleName +'_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq >> " + sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq; "+"cat " + sampleName +'_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq >> " + sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq; " +"cat " + sampleName +'_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_s.fastq >> " + sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_s.fastq" 
				os.system(cmd)
				cmd = "rm " + sampleName + '_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq " + sampleName + '_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq " + sampleName + '_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + "_s.fastq"
				os.system(cmd)

		elif len(fileList) == 1:
			#rm_file(sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq")
			for i in xrange(1, int(options.threads)+1):
				cmd = "cat " + sampleName + '_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq >> " + sampleName + "_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq"
				os.system(cmd)
				cmd = "rm " + sampleName + '_' + str(i) + "_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq"
				os.system(cmd)
		else:
			print "no or more than two fastq files for sample", sampleName
			sys.exit(1)
	"""


	#generate fileList for trimed fastq files
	fqTrimList = ""
	if len(fileList) == 2:
		if fileList[0].split('.')[-1] == "gz" and options.threads == 1:
			fqTrimList=CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq.gz,"+CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq.gz"
		else:
			fqTrimList=CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_1.fastq,"+CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + "_2.fastq"

	elif len(fileList) == 1:
		if fileList[0].split('.')[-1] == "gz" and options.threads == 1:
			fqTrimList=CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq.gz"
		else:
			fqTrimList=CWD+'/'+sampleName +"_Q" + str(htqc_C) + "L" + str(htqc_L) + ".fastq"
	else:
		print "no or more than two fastq files for sample", sampleName
		sys.exit(1)
	#print fqTrimList	

	#load in trimed sample file for quality control
	fileList = fqTrimList.split(",")
	#print "File List prior to fastqc run after trimming: ", fileList
	#quality check anagin after triming
	fastqc_run(sampleName+"_clean", fileList, outputFolder, options.threads)
	qc_plot(sampleName+"_clean", fileList, outputFolder, water)
