#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/group/corces/user/cpetersen/
#$ -e /wynton/group/corces/user/cpetersen/
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=60G
#$ -l scratch=50G
#$ -l h_rt=48:00:00

module load CBI r

Rscript /gladstone/corces/lab/users/cpetersen/Scripts/run_combineTrees.R

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
