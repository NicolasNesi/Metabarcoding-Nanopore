#!/bin/bash
#SBATCH -J P02_CleanStats                             # Job name
#SBATCH -o P02_CleanStats.%a.%A.out                   # File to which stdout will be written
#SBATCH -e P02_CleanStats.%a.%A.err                   # File to which stderr will be written
#SBATCH -p court                                      # Partition
#SBATCH -N 1                                          # Number of node - 1 GENOA node has 96-core and 750 Go RAM
#SBATCH -n 1                                          # Number of cores/cpus
#SBATCH -c 24                                         # Number of cpus ask per task
#SBATCH -t 00-12:00                                   # Runtime in DD-HH:MM
#SBATCH --mem 120G                                    # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                               # BEGIN,END,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr                            # Email address

# ---------------------------------
# Version v2.3
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date: 03/07/2025
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
MINLENGTH="${4}"
MAXLENGTH="${5}"
QSCORE="${6}"
FWD_PRIMER="${7}"
REV_PRIMER="${8}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="$FOLDER/RawReads_Nanopore/$RUN"
PATH2="$FOLDER/CleanedReads_Nanopore/$RUN"
PATH3="$FOLDER/QC_Stats/$RUN"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" ${FOLDER}/Scripts/List_samples_$RUN.txt)"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/"${SAMPLE}".fastq.gz "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/
gunzip *.fastq.gz

echo -ne Analysis run the: "${RUNDATE}"
echo -ne Working directory: "${PWD}"
echo -ne Submitted with Variables: Group "${GROUP}" Project "${PROJECT}" Run "${RUN}" Minimum read length "${MINLENGTH}" Maximum read length "${MAXLENGTH}" Qscore "${QSCORE}"
echo -ne ""
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate clean-ONT_env
# include:
# cutadapt=5.0
# fastp=0.23.4 
# nanoq=0.10.0 
# porechop_ABI=0.5.0 
# minimap2=2.26-r1175
# seqkit=2.6.1
# ---------------------------------

# ---------------------------------
# Adapter trimming
# ---------------------------------
echo "Starting Nanopore adaptor trimming for sample: $SAMPLE"

porechop_abi -abi \
-i "${SAMPLE}".fastq \
-t "${THREADS}" \
-v 0 \
-o "${SAMPLE}"_TrimAdapter.fastq.gz 2>&1 | tee -a "${SAMPLE}"_porechop-abi_run"${RUNDATE}".log
# ---------------------------------

# ---------------------------------
# Primer trimming
# ---------------------------------
echo "Starting primeres triming of Nanopore reads for sample: $SAMPLE"

# Producing the reverse complement of the primers provided by the user
RC_FWD_PRIMER=$(echo -ne ">fwd_primer\n${FWD_PRIMER}\n" | seqkit seq -r -p -t DNA 2>/dev/null | grep -v ">")
RC_REV_PRIMER=$(echo -ne ">rev_primer\n${REV_PRIMER}\n" | seqkit seq -r -p -t DNA 2>/dev/null | grep -v ">")

cutadapt \
"${SAMPLE}"_TrimAdapter.fastq.gz \
-g "${FWD_PRIMER}" \
-a "${RC_REV_PRIMER}" \
-e 0 \
-n 2 \
-j "${THREADS}" \
-o "${SAMPLE}"_cutadapt_output_tmp.fastq > /dev/null 2>/dev/null \
&& cutadapt \
"${SAMPLE}"_cutadapt_output_tmp.fastq \
-g "${REV_PRIMER}" \
-a "${RC_REV_PRIMER}" \
-e 0 \
-n 2 \
-j "${THREADS}" \
-o "${SAMPLE}"_cutadapt_output.fastq >/dev/null 2>/dev/null \
&& rm "${SAMPLE}"_cutadapt_output_tmp.fastq

mv "${SAMPLE}"_cutadapt_output.fastq "${SAMPLE}"_TrimmedReads.fastq
# ---------------------------------

# ---------------------------------
# Size and quality trimming
# ---------------------------------
echo "Starting quality and size trimming of Nanopore reads for sample: $SAMPLE"

fastp \
--in1 "${SAMPLE}"_TrimmedReads.fastq \
--disable_adapter_trimming \
--length_required "${MINLENGTH}" \
--length_limit "${MAXLENGTH}" \
--qualified_quality_phred "${QSCORE}" \
--low_complexity_filter \
--thread 16 \
--html "${SAMPLE}"_ReportFastp.html \
--report_title "$SAMPLE fastq report" \
--out1 "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats.fastq 2>&1 | tee -a "${SAMPLE}"_fastp_run"${RUNDATE}".log
#--dedup
# ---------------------------------

# ---------------------------------
# Remove simple repetitive sequences (20 stretches bases)
# ---------------------------------
echo "Remove sequences with simple repetitive motifs of 20 or more nucleotide on sample: $SAMPLE"

seqkit grep \
-vsif /home/2019013/PARTAGE/Various_Scripts_for_various_purposes/StretchPatterns_Nanopore_20bp.txt \
"${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats.fastq > "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq

pigz -f "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq
pigz -f "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats.fastq
# ---------------------------------

# ---------------------------------
# Statistic Nanopore reads
# ---------------------------------
echo "Starting statistic of Nanopore raw fastq for sample: $SAMPLE"

nanoq \
-i "${SAMPLE}".fastq \
--stats \
--header \
-vvv > "${SAMPLE}"_RawReads_Report.txt

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "${SAMPLE}"_RawReads_Report.txt
sed -i '2,5d' "${SAMPLE}"_RawReads_Report.txt
sed -i 's/_RawReads_Report//' "${SAMPLE}"_RawReads_Report.txt
sed -i 's/ \+ /\t/g' "${SAMPLE}"_RawReads_Report.txt

echo "Starting statistic of cleaned fastq with repeat for sample: $SAMPLE"

nanoq \
-i "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats.fastq.gz \
--stats \
--header \
-vvv > "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt
sed -i '2,5d' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt
sed -i "s/^$SAMPLE.*/$SAMPLE/" "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt
sed -i 's/ \+ /\t/g' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt

echo "Starting statistic of cleaned fastq no repeat for sample: $SAMPLE"

nanoq \
-i "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq.gz \
--stats \
--header \
-vvv > "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt

gawk -i inplace 'BEGINFILE {split(FILENAME,f,".")} NR==1{print f[1];print;next} {print}' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt
sed -i '2,5d' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt
sed -i "s/^$SAMPLE.*/$SAMPLE/" "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt
sed -i 's/ \+ /\t/g' "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt
# ---------------------------------

# ---------------------------------
# Create output directories and move outputs
# ---------------------------------
# Cleaned Fastq

if [[ ! -d "${PATH2}/withRepeats" ]]; then
   mkdir -p "${PATH2}/withRepeats"
fi

mv "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats.fastq.gz "${PATH2}/withRepeats/"

if [[ ! -d "${PATH2}/NoRepeats" ]]; then
   mkdir -p "${PATH2}/NoRepeats"
fi

mv "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq.gz "${PATH2}/NoRepeats/"

# Statistic
if [[ ! -d "${PATH3}/RawReads" ]]; then
  mkdir -p "${PATH3}/RawReads"
fi

mv "${SAMPLE}"_RawReads_Report.txt "${PATH3}"/RawReads/

if [[ ! -d "${PATH3}/CleanedReads/withRepeats" ]]; then
  mkdir -p "${PATH3}/CleanedReads/withRepeats"
fi

mv "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_withRepeats_Report.txt "${PATH3}/CleanedReads/withRepeats/"
mv "${SAMPLE}"_ReportFastp.html "${PATH3}/CleanedReads/"

if [[ ! -d "${PATH3}/CleanedReads/NoRepeats" ]]; then
  mkdir -p "${PATH3}/CleanedReads/NoRepeats"
fi

mv "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}"_noRepeats_Report.txt "${PATH3}/CleanedReads/NoRepeats/"

if [[ ! -d "${FOLDER}/Log/${RUN}_run${RUNDATE}" ]]; then
   mkdir -p "${FOLDER}/Log/${RUN}_run${RUNDATE}"
fi

# Log files
mv *.log "${FOLDER}/Log/${RUN}_run${RUNDATE}"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P02_CleanStats."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}"/Scripts/Sbatch_log
# ---------------------------------
