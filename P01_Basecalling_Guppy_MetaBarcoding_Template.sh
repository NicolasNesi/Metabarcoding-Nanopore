#!/bin/bash
#SBATCH -J P01_GuppyGPU                           # Job name
#SBATCH -o P01_GuppyGPU.%a.%A.out                 # File to which stdout will be written
#SBATCH -e P01_GuppyGPU.%a.%A.err                 # File to which stderr will be written
#SBATCH -p gpu_all                                # Partition
#SBATCH -N 1                                      # Number of node - 1 Milan node has 32-core, 500 Go RAM and GPU Ampere A100
#SBATCH --gres gpu:8                              # Number GPU
#SBATCH --ntasks-per-node 1                       # Number of cores/cpus
#SBATCH --cpus-per-task 28                        # Number of cpus per task
#SBATCH -t 00-20:00                               # Runtime in DD-HH:MM
#SBATCH --mem 100G                                # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                           # BEGIN,END,ALL
#SBATCH --mail-user=mymail                       # Email address

# ---------------------------------
# Version v2.3
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 05/02/2025
# ---------------------------------

# ---------------------------------
# Modules
# ---------------------------------
module purge
module load bio/guppy/6.5.7_gpu
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
CONFIGURATION="${4}"
BARCODE="${5}"
QSCORE="${6}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="/dlocal/home/2019013/Data/$GROUP/$PROJECT/Scripts"
PATH2="/dlocal/home/2019013/Data/$GROUP/$PROJECT/RawReads_Nanopore/$RUN"

RUNDATE="$(date '+%Y-%m-%d')"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp -R "${FOLDER}"/Nanopore_outputRun/"${RUN}"/no_sample_id/*/fast5 "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/

echo Analysis run the: "${RUNDATE}"
echo Working directory: "${PWD}"
echo Submitted with variables:: Group "${GROUP}" Project "${PROJECT}" Run "${RUN}" Guppy configuration "${CONFIGURATION}" NanoporeKit "${BARCODE}" Qscore "${QSCORE}"
# ---------------------------------

# ---------------------------------
# Make subdirectories
# ---------------------------------
for dir in \
   "${FOLDER}/Guppy_Basecalling/${RUN}" \
   "${FOLDER}/RawReads_Nanopore/${RUN}" \
   "${FOLDER}/CleanedReads_Nanopore" \
   "${FOLDER}/Microbiome_R/Phyloseq" \
   "${FOLDER}/Microbiome_R/Metadata" \
   "${FOLDER}/Microbiome_R/BIOM" \
   "${FOLDER}/QC_Stats/${RUN}" \
   "${FOLDER}/Log/${RUN}_run${RUNDATE}" \
   "${FOLDER}/Scripts/Sbatch_log"
do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Created directory: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done
# ---------------------------------

# ---------------------------------
# Guppy GPU - Trimming (barcodes, adapters) - demultiplexing
# ---------------------------------
guppy_basecaller \
-i "${LOCAL_WORK_DIR}"/fast5 \
-s "${LOCAL_WORK_DIR}"/demultiplex \
-c "${CONFIGURATION}" \
-x "cuda:0,1" \
--num_callers 28 \
--gpu_runners_per_device 8 \
--chunks_per_runner 512 \
--chunk_size 3000 \
--trim_strategy dna \
--trim_adapters \
--detect_mid_strand_barcodes \
--enable_trim_barcodes \
--barcode_kits "${BARCODE}" \
--min_qscore "${QSCORE}"
# ---------------------------------

# ---------------------------------
# Concatenate pass fastq
# ---------------------------------
if [[ -d demultiplex/pass/unclassified ]]; then
  mv demultiplex/pass/unclassified demultiplex/pass/barcode00
fi

for barcode_dir in demultiplex/pass/barcode* ; do
  if [[ -d "${barcode_dir}" ]]; then
   cat "${barcode_dir}"/*.fastq > "${barcode_dir}".fastq
   rm -r "${barcode_dir}"
   mv "${barcode_dir}".fastq demultiplex/
  fi
done
# ---------------------------------

# ---------------------------------
# Move demultiplex outputs
# ---------------------------------
mv demultiplex/ "${FOLDER}"/Guppy_Basecalling/"${RUN}"
cd "${FOLDER}"/Guppy_Basecalling/"${RUN}"/demultiplex/
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate NN_env
# include:
# csvtk=0.31.0
# ---------------------------------

# ---------------------------------
# Rename samples using Sample_MetaBarcodingGuppy_${RUN}.xlsx
# ---------------------------------
# Convert XLSX to CSV and remove the first line (header)
csvtk xlsx2csv "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.xlsx" | tail -n +2 > "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv.temp"

# Remove rows where the second column is empty or contains only whitespace
awk -F',' '!(($2 == "" || $2 ~ /^[[:space:]]*$/))' "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv.temp" > "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv"
rm "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv.temp"

# Remove the first column from Rename_RUN.txt and save to List_Samples_RUN.txt
awk -F',' 'NF > 1 && $2 !~ /^[[:space:]]*$/ { print $2 }' "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv" > "${PATH1}/List_samples_${RUN}.txt"

# Copy Rename_RUN.txt to the Demultiplex directory
cp "${PATH1}/Sample_MetaBarcodingGuppy_${RUN}.csv" "${FOLDER}"/Guppy_Basecalling/"${RUN}"/demultiplex/

# Move to the Demultiplex directory
cd "${FOLDER}"/Guppy_Basecalling/"${RUN}"/demultiplex/

# Prepare OldNewNames.txt with old and new filenames
awk -F',' '{ print $1".fastq", $2".fastq" }' "Sample_MetaBarcodingGuppy_${RUN}.csv" > OldNewNames.txt

# Rename files based on OldNewNames.txt
while read -r old new; do
   mv "$old" "$new"
done < OldNewNames.txt
# ---------------------------------

# ---------------------------------
# Move Fastq to RawReads
# ---------------------------------
mv "${FOLDER}"/Guppy_Basecalling/"${RUN}"/demultiplex/*.fastq "${PATH2}"/
rm "${PATH2}"/barcode*.fastq
pigz -f "${PATH2}"/*.fastq
# ---------------------------------

# ---------------------------------
# Quality control with PycoQC
# ---------------------------------
conda deactivate
conda activate QC-ONT_env
# include:
# PycoQC=2.5.2

pycoQC \
--summary_file "${FOLDER}"/Guppy_Basecalling/"${RUN}"/demultiplex/sequencing_summary.txt \
--report_title "${PROJECT}"_"${RUN}"_"${CONFIGURATION}"_minqscore"${QSCORE}"_report \
--html_outfile "${FOLDER}"/QC_Stats/"${RUN}"/"${PROJECT}"_"${RUN}"_run"${RUNDATE}".html
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P01_GuppyGPU.* "${FOLDER}"/Scripts/Sbatch_log/
# ---------------------------------
