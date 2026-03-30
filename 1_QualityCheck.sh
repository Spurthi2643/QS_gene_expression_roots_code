#!/bin/bash

######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
## 	Use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Download from SRA (completed in script 0_downloadSRA): Input Data: NA
## 			Output: Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
## FASTQC 	InPut: Downloaded SRA files .fastq
##		Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##			The last line of this script will make a tarball of the output directory to bring back to your computer
## 	After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
## 	then run the script by using "run_script [script name]"
## 	suggested paramenters are below to submit this script.
## 		queue: class
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 1gb
##		run on asax
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## These are represented in the code by [#] replace these according to the examples provided
MyID=aubnsa001         ## Example: MyID=aubrmg001

  ## Make variable that represent YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/aubnsa001/PracticeRNAseq/RawData_QSroot   			## Example: DD=/scratch/${MyID}/PracticeRNAseq/RawData
WD=/scratch/aubnsa001/PracticeRNAseq			## Example: WD=/scratch/${MyID}/PracticeRNAseq
RDQ=RawDataQSQuality
 
## move to the Data Directory
cd ${DD}

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
##Note, if the parent/child directories already exists we don't use the -p 
##CHECK if the parent/child directories exist - if they do not you will need to add -p

mkdir ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd ${WD}/${RDQ}
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.