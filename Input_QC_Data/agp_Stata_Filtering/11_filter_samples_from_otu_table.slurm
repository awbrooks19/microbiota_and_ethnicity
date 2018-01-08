#!/bin/bash
#SBATCH --mail-user=andrew.w.brooks@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00-2:00:00
#SBATCH --mem=4G
#SBATCH --output=11_filter_samples_from_otu_table.out

##########################################################################################
############################   HELPFUL INFORMATION   #####################################

#############################################
### TO RUN THE SCRIPT #######################
# sbatch accre_submission.slurm

#############################################
### TO MONITOR SUBMITTED JOBS ###############

### HELP ### 
# man squeue

### SHOW ### 
# squeue

### SHOW w/ OPTIONS ###
# squeue –user=<vunetid>
# squeue -user=brooksaw --Format=jobid,username,account,statecompact,starttime,timelimit,numcpus

#############################################
### TO CANCEL THE SCRIPT ####################
# scancel jobID

#############################################
### TO CHECK JOB USAGE STATISTICS ###########
# sacct
# sacct --jobs=jobID
# sacct --starttime=mm.dd.yy
# sacct --jobs=27 --format=User,JobID,account,Timelimit,elapsed,ReqMem,MaxRss,ExitCode
# rtracejob <jobID>

##########################################################################################
#############################   RUNNING CODE   ###########################################
#############################################
### DESCRIPTION ###
#
# Filter samples from OTU table based on those that remain
# after removing basic metadata (BMI,RACE,AGE_YEARS,SEX,LOCATION,SAMPLE_SOURCE) from the American Gut dataset
# 7/11/16
#
#
#############################################

filter_samples_from_otu_table.py -i 01_otu_table.biom -o 11_filter_samples_from_otu_table.biom --sample_id_fp=11_samples_to_keep.txt
