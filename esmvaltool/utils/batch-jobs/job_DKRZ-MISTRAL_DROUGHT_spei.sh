#!/bin/bash -e
###############################################################################
### BATCH SCRIPT TO RUN THE ESMVALTOOL AT DKRZ MISTRAL
### Author: Mattia Righi (DLR)
###############################################################################
#SBATCH --partition=compute
#SBATCH --ntasks=8
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --account=bd1083
#SBATCH --output=/work/bd0854/b380745/student/ESMValTool/job_%j.out.log
#SBATCH --error=/work/bd0854/b380745/student/ESMValTool/job_%j.out.log
###############################################################################

# Submit job with: sbatch job_DKRZ-MISTRAL_DROUGHT_spei.sh

# Input arguments
RECIPE=recipe_speikeah2.yml
CONFIG=config-user_keah.yml

# Set environment
CONDAPATH=/pf/b/b380745/anaconda3  # e.g. /home/soft/miniconda3/
CONDAENV=$CONDAPATH/envs/esmvaltool/bin   # e.g. $CONDAPATH/envs/esmvaltool/bin
ESMVALPATH=/work/bd0854/b380745/student/ESMValTool/esmvaltool # e.g. /home/ESMValTool/esmvaltool

# Changes below this line should not be required
export PATH=$PATH:$CONDAPATH/bin/
conda info --envs
$CONDAENV/esmvaltool $ESMVALPATH/recipes/$RECIPE -c $ESMVALPATH/$CONFIG
