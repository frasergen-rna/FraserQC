# FraserQC
A quality control tool for high throughput sequence data

## fraserQC.py
Data quality control and data statistics for high throughput sequences

Usage:
`python fraserQC.py [options] -i fastqFile -s sampleName -o outputFolder`

Options:

    --version                               show program's version number and exit
    -h, --help                              show this help message and exit
    -s SAMPLENAME, --sampleName=SAMPLENAME  [REQUIRED] sample name for fastq files
    -i FQLIST, --fqList=FQLIST              [REQUIRED] input fastq files ended with fq/fastq/gz (PEs were separated by ",")
    -o FOLDER, --outputFolder=FOLDER        [REQUIRED] In this folder, the program will create sample-specific folder, containning all qc result 
    -f, --filter                            perform quality trim by htqc [default to skip this step]
    -m, --watermarker                       print out water marker [default to not print water marker]

  Advanced Options:

    -w INT              window size for quality trim [default: 5]
    -C INT              quality threshold for trim [default: 20]
    -L INT              length threshold for trim [default: 50]

Example:
`python fraserQC.py -f -L 100 -C 30 -w 10 -i test_1.fq.gz,test_2.fq.gz -t 2 -s test -o outdir`

## Dependencies:
- FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- htqc (https://sourceforge.net/projects/htqc/)
- R

## How to use this tool:
1. Install FastQC, htqc and R
2. Modify fraserQC.py with the install path of FastQC and htqc. Specifically modify the script at line 417 and 421

