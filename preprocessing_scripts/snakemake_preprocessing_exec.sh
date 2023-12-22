#!/bin/bash
#Cathrine Petersen 11/20/22
#
#
#

#PIPELINE VARIABLES:

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" #retrieve the source directory of this pipeline script
IN_FLAG=0 #flag to track if input directory was provided. If not, pipeline will fail because this is a required argument
OUT_FLAG=0 #flag to track if an output directory was provided. If not, pipeline will fail because this is a required argument
METADATA_FLAG=0 #flag to track if a metadata file was provided. If not, pipeline will fail because this is a required argument
TEMP_DIR=/wynton/group/muckelab/user/cpetersen/cluster_benchmarking/temp
FORCE=1 #flag to track if temporary directory checks can be ignored (risks overwriting)
REMOVE=0 #flag to track whether the temporary directory should be removed at the end of the pipeline run
EMAIL=""
EMAIL_FLAG=0

#COMMAND LINE OPTION INTERPRETATION:
#Handle command line inputs using getopts
while getopts ":hFNRi:s:o:d:m:e:" opt; do
	case $opt in
	h)
		echo -e "\nUsage: bash /path/to/snakemake_preprocessing_exec.sh -m <Input_Manifest> -o <Output_Directory> ... <other options>\n"
		echo -e "Error, incorrect usage."
		exit 0
		;;
	i)
		if [[ -d ${OPTARG} ]]; then
			INPUT_DIR=${OPTARG}
			IN_FLAG=1
			echo "-i flag observed. INPUT_DIR set to ${OPTARG}." >&2
		else
			echo "ERROR --- -i flag observed but suggested INPUT_DIR does not exist: ${OPTARG}" >&2
			exit 1
		fi
		;;
	s)
		if [[ -d ${OPTARG} ]]; then
			SCRIPT_DIR=${OPTARG}
			echo "-s flag observed. Replacing default SCRIPT_DIR with ${OPTARG}." >&2
		else
			echo "ERROR --- -s flag observed but suggested SCRIPT_DIR does not exist: ${OPTARG}" >&2
			exit 1
		fi
		;;
	o)
		OUTPUT_DIR=${OPTARG}
		mkdir -p ${OUTPUT_DIR}
		OUT_FLAG=1
		echo "-o flag observed. OUTPUT_DIR set to ${OPTARG}." >&2
		;;
	d)
		TEMP_DIR=${OPTARG}
		echo "-d flag observed. Replacing default TEMP_DIR with ${OPTARG}." >&2
		;;
	m)
		METADATA_FILE=${OPTARG}
		METADATA_FLAG=1
		echo "-m flag observed. METADATA_FILE set to ${OPTARG}." >&2
		;;
	e)
		EMAIL=${OPTARG}
		EMAIL_FLAG=1
		echo "-e flag observed. Email notification will go to ${OPTARG}." >&2
		;;
	F)
		echo "-F flag observed. Temporary directory may be overwritten!" >&2
		FORCE=0
		;;
	N)
		echo "-N flag observed. Temporary directory will not be deleted at the end of the run!" >&2
		RSYNC=1
		;;
	\?)
		echo "Invalid option: -$OPTARG. Check the README.md file for detailed usage instructions!" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument. Check the README.md file for detailed usage instructions!" >&2
		exit 1
		;;
	*)
		echo "Unimplemented option: -$OPTARG. Check the README.md file for detailed usage instructions!" >&2
		exit 1
		;;
	esac
done

#Submission strings for snakemake --cluster parameter
#If email was entered
if [ "$EMAIL_FLAG" == 1 ]; then
	echo "Email was provided."
	QSUB="--jobs 200 --cluster "\"'qsub -cwd -pe smp {threads} -l mem_free={resources.mem_qsub} -l h_rt={resources.job_time} -j yes -V -o ${LOG_DIR} -e ${LOG_DIR} -m bea -M ${EMAIL}'\"
else
	QSUB="--jobs 200 --cluster "\"'qsub -cwd -pe smp {threads} -l mem_free={resources.mem_qsub} -l h_rt={resources.job_time} -j yes -V -o ${LOG_DIR} -e ${LOG_DIR}'\"
fi
CLUSTER_SUB=${QSUB}

#Check for presence of input and output directory from command line input
if [ "$IN_FLAG" != 1 ]; then
	echo "No valid input directory supplied. Check -i argument!"
	exit 1
fi
if [ "$OUT_FLAG" != 1 ]; then
	echo "No output directory provided. -o argument is required!"
	exit 1
fi
if [ "$METADATA_FLAG" != 1 ]; then
	echo "No metadata file provided. -m argument is required!"
	exit 1
fi

#Check for validity of TEMP_DIR
#If FORCE is FALSE, check if directory exists. If so quit.
if [[ "$FORCE" -eq 1 ]]; then
	if [[ -d ${TEMP_DIR} ]]; then
		echo "ERROR --- TEMP_DIR already exist: ${TEMP_DIR}" >&2
		echo "Pipeline expects that the temporary directory does not exist to prevent overwriting of data or collision of multiple pipeline runs" >&2
		echo "Either (i) use the -F flag to override this, (ii) delete the suggested directory, or (iii) provide a path to a directory that does not already exist using the -d argument." >&2
		exit 1
	else
		mkdir -p ${TEMP_DIR}
	fi
else
	mkdir -p ${TEMP_DIR}
fi

#Fix all provided directory paths to remove trailing forward slashes ("/"). This prevents snakemake from complaining about double slashes ("//") in paths.
INPUT_DIR=$(echo "$INPUT_DIR" | sed 's:/*$::')
TEMP_DIR=$(echo "$TEMP_DIR" | sed 's:/*$::')
OUTPUT_DIR=$(echo "$OUTPUT_DIR" | sed 's:/*$::')

echo -e "\n\n"
echo "Input Directory: ${INPUT_DIR}"
echo "Temporary Directory: ${TEMP_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"

#Pause for user to see pipeline messages
sleep 2

#Designate a directory within TEMP_DIR to store log files
LOG_DIR=${TEMP_DIR}/logs
mkdir -p ${LOG_DIR}

#Create a directory within OUTPUT_DIR for intermediate files (if it doesn't already exist)
INTERMEDIATE_DIR=${OUTPUT_DIR}/intermediate_files
mkdir -p ${INTERMEDIATE_DIR}

#Create a directory within OUTPUT_DIR for final output files (if it doesn't already exist)
FINAL_DIR=${OUTPUT_DIR}/final_output
mkdir -p ${FINAL_DIR}

#change directory to the temporary directory
cd ${TEMP_DIR}

#Run snakemake command
COMMAND="snakemake --rerun-incomplete --snakefile ${SCRIPT_DIR}/snakemake_preprocessing.py --config in_dir=${INPUT_DIR} out_dir=${OUTPUT_DIR} \
temp_dir=${TEMP_DIR} scripts=${SCRIPT_DIR} metadata=${METADATA_FILE} --keep-going ${CLUSTER_SUB} --latency-wait 120"

echo $COMMAND

eval ${COMMAND}

if [ $? -eq 0 ]; then
	if [[ "$REMOVE" -eq 0 ]]; then
		echo -e "\n\n##########################################################################################################"
		echo Snakemake completed successfully! Removing temporary directory.
		echo -e "##########################################################################################################\n\n"
		rm -r ${TEMP_DIR}
	else
		echo -e "\n\n##########################################################################################################"
		echo Snakemake completed successfully! -N flag observed so temporary directory will not be removed.
		echo -e "##########################################################################################################\n\n"
	fi
else
	echo -e "\n\n##########################################################################################################"
	echo Snakemake failed. Temporary directory is left in place.
	echo -e "##########################################################################################################\n\n"
fi












