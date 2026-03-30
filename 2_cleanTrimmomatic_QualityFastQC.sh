#!/bin/bash

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## FASTQC output is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##		Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##			The last line of this script will make a tarball of the output directory to bring back to your computer
## For running the script on the Alabama Super Computer.
                ##For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
        ## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
        ## then run the script by using "run_script [script name]"
        ## suggested paramenters are below to submit this script.
                ## queue: medium 
                ## core: 6
                ## time limit (HH:MM:SS): 02:00:00  (may need to increase, if so run on medium queue)
                ## Memory: 12gb
                ## run on asax
###############################################

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

# Modules
	#  load the module
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

## STOP. You need to replace the [number] with YOUR paths to 
##       make variables for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubnsa001                       ## Example: MyID=aubrmg001

# Variables: raw data directory (DD), working directory(WD), Quality after cleaning (PCQ), name of file containing the adpaters.
WD=/scratch/aubnsa001/PracticeRNAseq                          ## Example: WD=/scratch/$MyID/PracticeRNAseq
DD=/scratch/aubnsa001/PracticeRNAseq/RawData_QSroot                         ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
CD=/scratch/aubnsa001/PracticeRNAseq/CleanData                         ## Example: CD=/scratch/$MyID/PracticeRNAseq/CleanData
PCQ=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
				## In the future, for your data, you will likely need to edit this for other projects based on how your libraries 
				## were made to search for the correct adapters for your project
				
## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
## Note, if the parent directory exists do not use -p
mkdir ${CD}
mkdir ${WD}/${PCQ}

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls *.fastq | sed 's/\.fastq$//' | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
cp /home/aubnsa001/scripts/AdaptersToTrim_All.fa .

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in single end (SE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
	#From Brittany
                 #/apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar 
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
	SE -threads 6 -phred33 \
        "${i}.fastq" \
        ${CD}/${i}_trimmed.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
        ## Trim single-end reads for adapters and low-quality bases.
        ## SE = single-end mode, with one input FASTQ and one output FASTQ.
        ## -threads 6 uses 6 CPU threads, and -phred33 specifies Phred+33 quality encoding.
        ## ILLUMINACLIP removes adapter contamination.
        ## HEADCROP:10 removes the first 10 bases of each read.
        ## LEADING:30 and TRAILING:30 trim bases below Q30 from the ends.
        ## SLIDINGWINDOW:6:30 trims when the average quality in a 6-base window drops below Q30.
        ## MINLEN:36 removes reads shorter than 36 bp after trimming.
	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

fastqc ${CD}/${i}_trimmed.fastq --outdir=${WD}/${PCQ}

done<list			# This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd ${WD}/${PCQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*

## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate the quality of the data.
