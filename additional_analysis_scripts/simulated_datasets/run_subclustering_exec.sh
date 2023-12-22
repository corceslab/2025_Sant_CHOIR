#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/group/corces/user/cpetersen/
#$ -e /wynton/group/corces/user/cpetersen/
#$ -cwd
#$ -pe smp 16
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l scratch=25G
#$ -l h_rt=06:00:00

cd /wynton/group/muckelab/user/cpetersen/cluster_benchmarking

. bin/activate

module load CBI r/4.2.3

f=$(basename $1)

Rscript /gladstone/corces/lab/users/cpetersen/Scripts/run_subclustering.R ${f}

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
