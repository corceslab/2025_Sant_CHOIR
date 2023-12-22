#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/group/muckelab/user/cpetersen/cluster_benchmarking/temp
#$ -e /wynton/group/muckelab/user/cpetersen/cluster_benchmarking/temp
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=4G
#$ -l scratch=25G
#$ -l h_rt=04:00:00

module load CBI r/4.2.3

f=$(basename $1)

Rscript generate_simulated_data.R $f

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
