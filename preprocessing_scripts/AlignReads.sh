#!/bin/bash
#Cathrine Petersen 04/26/22
#AlignReads.sh

export PATH=/wynton/group/muckelab/shared/tools/cellranger-arc-2.0.1:$PATH
export PATH=/wynton/group/muckelab/shared/tools/cellranger-atac-2.1.0:$PATH

#VARIABLES:
IN_FLAG=0 #flag to track if input directory was provided. If not, script will fail, because this is a required argument.
SPECIES_FLAG=0 #flag to track if a species was provided. If not, script will fail because this is a required argument
CELLS_FLAG=0 #flag to track if cells/nuclei info was provided. If not, script will fail because this is a required argument
SAMPLE_FLAG=0 #flag to track if a sample ID was provided. If not, script will fail because this is a required argument
TECH_FLAG=0 #flag to track if a tech was provided. If not, script will fail because this is a required argument

#COMMAND LINE OPTION INTERPRETATION:
#Handle command line inputs using getopts
while getopts "i:a:c:s:t:" opt; do
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
	a)
		SPECIES=${OPTARG}
		SPECIES_FLAG=1
		echo "-a flag observed. SPECIES set to ${OPTARG}." >&2
		;;
	c)
		CELLS_NUCLEI=${OPTARG}
		CELLS_FLAG=1
		echo "-c flag observed. CELLS_NUCLEI set to ${OPTARG}." >&2
		;;
	s)
		SAMPLE=${OPTARG}
		SAMPLE_FLAG=1
		echo "-s flag observed. SAMPLE set to ${OPTARG}." >&2
		;;
	t)
		TECH=${OPTARG}
		TECH_FLAG=1
		echo "-t flag observed. TECH set to ${OPTARG}." >&2
		;;
	esac
done

#Check for presence of required arguments from command line input
if [ "$IN_FLAG" != 1 ]; then
	echo "No valid input directory supplied. Check -i argument!"
	exit 1
fi
if [ "$SPECIES_FLAG" != 1 ]; then
	echo "No species provided. -a argument is required!"
	exit 1
fi
if [ "$CELLS_FLAG" != 1 ]; then
	echo "No cells/nuclei info provided. -c argument is required!"
	exit 1
fi
if [ "$SAMPLE_FLAG" != 1 ]; then
	echo "No sample ID provided. -s argument is required!"
	exit 1
fi

# Based on species & tech, designate a reference genome
if [[ "$TECH" = "RNA" ]]; then
	HUMAN_REFERENCE=/wynton/group/muckelab/shared/genomes/hg38/cellranger/refdata-gex-GRCh38-2020-A
	MOUSE_REFERENCE=/wynton/group/muckelab/shared/genomes/mm10/cellranger/refdata-gex-mm10-2020-A
elif [[ "$TECH" = "ATAC" ]]; then
	HUMAN_REFERENCE=/wynton/group/muckelab/shared/genomes/hg38/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
	MOUSE_REFERENCE=/wynton/group/muckelab/shared/genomes/mm10/cellranger/refdata-cellranger-arc-mm10-2020-A-2.0.0
elif [[ "$TECH" = "RNA_ATAC" ]]; then
	HUMAN_REFERENCE=/wynton/group/muckelab/shared/genomes/hg38/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
	MOUSE_REFERENCE=/wynton/group/muckelab/shared/genomes/mm10/cellranger/refdata-cellranger-arc-mm10-2020-A-2.0.0
fi

if [[ "$SPECIES" = "human" ]]; then
	REFERENCE=$HUMAN_REFERENCE
elif [[ "$SPECIES" = "mouse" ]]; then
	REFERENCE=$MOUSE_REFERENCE
fi


echo "REFERENCE = $REFERENCE"

#RUN COMMAND

module load CBI
module load cellranger

if [[ "$CELLS_NUCLEI" = "nuclei" ]] && [[ "$TECH" = "RNA" ]]; then
	echo "CELLS_NUCLEI = $CELLS_NUCLEI"
	echo "TECH = $TECH"
	echo "Running cellranger for RNA nuclei"
	cellranger count --id=$SAMPLE \
		--fastqs=${INPUT_DIR}/fastq_files \
		--sample=$SAMPLE \
		--include-introns=true \
		--transcriptome=$REFERENCE \
		--localcores=16 \
		--localmem=50
elif [[ "$CELLS_NUCLEI" = "cells" ]] && [[ "$TECH" = "RNA" ]]; then
	echo "CELLS_NUCLEI = $CELLS_NUCLEI"
	echo "TECH = $TECH"
	echo "Running cellranger for RNA cells"
	cellranger count --id=$SAMPLE \
		--fastqs=${INPUT_DIR}/fastq_files \
		--sample=$SAMPLE \
		--transcriptome=$REFERENCE \
		--localcores=16 \
		--localmem=50
elif [[ "$TECH" = "ATAC" ]]; then
	echo "TECH = $TECH"
	echo "Running cellranger for ATAC nuclei"
	/wynton/group/muckelab/shared/tools/cellranger-atac-2.1.0/cellranger-atac count --id=$SAMPLE \
		--reference=$REFERENCE \
		--fastqs=${INPUT_DIR}/fastq_files \
		--sample=$SAMPLE \
		--localcores=16 \
		--localmem=50
elif [[ "$CELLS_NUCLEI" = "cells" ]] && [[ "$TECH" = "RNA_ATAC" ]]; then
	echo "CELLS_NUCLEI = $CELLS_NUCLEI"
	echo "TECH = $TECH"
	echo "Running cellranger for ATAC + RNA cells"
	/wynton/group/muckelab/shared/tools/cellranger-arc-2.0.1/cellranger-arc count --id=$SAMPLE \
		--reference=$REFERENCE \
		--libraries=${INPUT_DIR}/${SAMPLE}_libraries.csv \
		--localcores=16 \
		--localmem=50
		--gex-exclude-introns
elif [[ "$CELLS_NUCLEI" = "nuclei" ]] && [[ "$TECH" = "RNA_ATAC" ]]; then
	echo "CELLS_NUCLEI = $CELLS_NUCLEI"
	echo "TECH = $TECH"
	echo "Running cellranger for ATAC + RNA nuclei"
	/wynton/group/muckelab/shared/tools/cellranger-arc-2.0.1/cellranger-arc count --id=$SAMPLE \
		--reference=$REFERENCE \
		--libraries=${INPUT_DIR}/${SAMPLE}_libraries.csv \
		--localcores=16 \
		--localmem=50
fi

if [[ "$TECH" = "RNA" ]]; then
	#Copy matrix files for raw count matrix
	cp ${SAMPLE}/outs/raw_feature_bc_matrix/* ${INPUT_DIR}/aligned_reads/${SAMPLE}
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/features.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/matrix.mtx.gz
	# Copy just cell ids for filtered count matrix
	cp ${SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	# Copy cellranger clustering assignments
	cp ${SAMPLE}/outs/analysis/clustering/gene_expression_graphclust/clusters.csv ${INPUT_DIR}/aligned_reads/${SAMPLE}/cellranger_clusters.csv
elif [[ "$TECH" = "ATAC" ]]; then
	#Copy matrix files for raw peak matrix
	cp ${SAMPLE}/outs/raw_peak_bc_matrix/* ${INPUT_DIR}/aligned_reads/${SAMPLE}
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/features.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/matrix.mtx.gz
	# Copy just cell ids for filtered peak matrix
	cp ${SAMPLE}/outs/filtered_peak_bc_matrix/barcodes.tsv.gz ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	# Copy ATAC peak bed file
	cp ${SAMPLE}/outs/peaks.bed ${INPUT_DIR}/aligned_reads/${SAMPLE}/peaks.bed
	#Copy ATAC fragment files
	cp ${SAMPLE}/outs/fragments.tsv.gz ${INPUT_DIR}/aligned_reads/${SAMPLE}/atac_fragments.tsv.gz
	cp ${SAMPLE}/outs/fragments.tsv.gz.tbi ${INPUT_DIR}/aligned_reads/${SAMPLE}/atac_fragments.tsv.gz.tbi
	# Copy cellranger clustering assignments
	cp ${SAMPLE}/outs/analysis/clustering/graphclust/clusters.csv ${INPUT_DIR}/aligned_reads/${SAMPLE}/cellranger_clusters.csv
elif [[ "$TECH" = "RNA_ATAC" ]]; then
	#Copy matrix files for raw count matrix
	cp ${SAMPLE}/outs/raw_feature_bc_matrix/* ${INPUT_DIR}/aligned_reads/${SAMPLE}
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/features.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/matrix.mtx.gz
	# Copy just cell ids for filtered count matrix
	cp ${SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	gunzip ${INPUT_DIR}/aligned_reads/${SAMPLE}/filtered_barcodes.tsv.gz
	# Copy cellranger clustering assignments
	cp ${SAMPLE}/outs/analysis/clustering/gex/graphclust/clusters.csv ${INPUT_DIR}/aligned_reads/${SAMPLE}/cellranger_clusters.csv
	# Copy ATAC peak bed file
	cp ${SAMPLE}/outs/atac_peaks.bed ${INPUT_DIR}/aligned_reads/${SAMPLE}/peaks.bed
	#Copy ATAC fragment files
	cp ${SAMPLE}/outs/atac_fragments.tsv.gz ${INPUT_DIR}/aligned_reads/${SAMPLE}/atac_fragments.tsv.gz
	cp ${SAMPLE}/outs/atac_fragments.tsv.gz.tbi ${INPUT_DIR}/aligned_reads/${SAMPLE}/atac_fragments.tsv.gz.tbi
fi

#Copy qc report
cp ${SAMPLE}/outs/web_summary.html ${INPUT_DIR}/quality_checks/${SAMPLE}_web_summary.html



