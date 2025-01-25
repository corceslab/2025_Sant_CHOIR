#!/bin/bash
#
#$ -S /bin/bash
#$ -o /wynton/group/corces/user/cpetersen
#$ -e /wynton/group/corces/user/cpetersen
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=50G
#$ -l scratch=100G
#$ -l h_rt=36:00:00

cd /wynton/group/muckelab/user/cpetersen/test_sctriangulate

. bin/activate

python3 /gladstone/corces/lab/users/cpetersen/Scripts/run_scTriangulate.py

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

