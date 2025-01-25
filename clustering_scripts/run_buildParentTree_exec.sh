#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/group/corces/user/cpetersen/
#$ -e /wynton/group/corces/user/cpetersen/
#$ -cwd
#$ -pe smp 16
#$ -r y
#$ -j y
#$ -l mem_free=16G
#$ -l scratch=25G
#$ -l h_rt=96:00:00

module load CBI r

Rscript /gladstone/corces/lab/users/cpetersen/Scripts/run_buildParentTree.R

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
