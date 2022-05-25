#!/bin/bash

#SBATCH --job-name=dynamic_protein_1ns                   #This is the name of your job
#SBATCH --cpus-per-task=4                  #This is the number of cores reserved
#SBATCH --mem=10G              #This is the memory reserved per core.
#Total memory reserved: 10GB

#SBATCH --time=10:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=output3_errors     #This is the joined STDOUT and STDERR file
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=manon.albecq2a@hotmail.fr        #You will be notified via email when your task ends or fails

##SBATCH --partition=rtx8000
#SBATCH --gres=gpu:1        # --gres=gpu:2 for two GPU, etc




#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#module load CUDA/11.3.1



#load your required modules below

module purge

export LMOD_DISABLE_SAME_NAME_AUTOSWAP='no'

#module load OpenMM/7.1.1-foss-2018b-Python-3.6.6
#module load --force Python/3.6.6-fosscuda-2018b 

module load nvptx-tools/20180301-GCCcore-7.3.0-CUDA-10.0.130 #CUDA/11.3.1 

#module load numpy

source ~/.bash_profile

conda activate py37

#add your command lines below

python prot_dyna.py
