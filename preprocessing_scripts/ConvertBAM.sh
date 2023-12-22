#!/bin/bash
#Cathrine Petersen 04/26/22
#ConvertBAM.sh

#VARIABLES:
IN_FLAG=0 #flag to track if input directory was provided. If not, script will fail, because this is a required argument.
SAMPLE_FLAG=0 #flag to track if a sample ID was provided. If not, script will fail because this is a required argument

#COMMAND LINE OPTION INTERPRETATION:
#Handle command line inputs using getopts
while getopts "i:s:" opt; do
	case $opt in
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
		SAMPLE=${OPTARG}
		SAMPLE_FLAG=1
		echo "-s flag observed. SAMPLE set to ${OPTARG}." >&2
		;;
	esac
done

#Check for presence of required arguments from command line input
if [ "$IN_FLAG" != 1 ]; then
	echo "No valid input directory supplied. Check -i argument!"
	exit 1
fi
if [ "$SAMPLE_FLAG" != 1 ]; then
	echo "No sample ID provided. -s argument is required!"
	exit 1
fi

# Set up input and outputs
bam_file=${INPUT_DIR}/bam_files/${SAMPLE}.bam
out_dir=${SAMPLE}_fastq

#RUN COMMAND

module load CBI
module load cellranger

cellranger bamtofastq --nthreads=16 $bam_file $out_dir

# Move & renamefiles
fastq_files=$(find ${SAMPLE}_fastq -name "*bamtofastq*")
for f in ${fastq_files} ; do mv -- "$f" "${INPUT_DIR}/fastq_files/${SAMPLE}${f##*bamtofastq}" ; done

# Write file names to text file
find ${INPUT_DIR}/fastq_files -type f -maxdepth 1 -name "*${SAMPLE}*" > "${INPUT_DIR}/quality_checks/${SAMPLE}_fastq_file_list.txt"
