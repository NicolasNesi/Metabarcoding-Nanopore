#!/usr/bin/env bash

# -----------------------------------------------
# Control file to config and submit the three scripts
# of the MetaBarcoding Pipeline using Nanopore Input data
# -----------------------------------------------

# -----------------------------------------------
# Version v1.4
# -----------------------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:03/07/2025
# -----------------------------------------------

# -----------------------------------------------#
# -----------------------------------------------#
# EDIT THE FOLLOWING LINES OF THE CONFIG SCRIPTS #
# -----------------------------------------------#
# -----------------------------------------------#

# -----------------------------------------------
# READ INSTRUCTION BEFORE SUBMISSION
# -----------------------------------------------
# Configuration file for the submission of the pipeline for MetaBarconding analyses using Nanopore amplicon with chimestry V14.

# Start point is a folder /dlocal/home/2019013/Data/$GROUP/$PROJECT with:

#   - subfolder Nanopore_outputRun/ where you will copy your Nanopore output folder of the $RUN

#   - subfolder Scripts/ with the MetaBarcoding Pipeline scripts, List_samples_$RUN.txt and Sample_sheet_$RUN.csv
# 
# Sample name format: Start with [a-zA-Z], NO space, NO underscore! (e.g. S23-MTX-COV-NOSE)
#
# .
# └── GROUP
#     └── PROJECT
#         ├── Nanopore_outputRun
#         │   └── RUN
#         │       └── no_sample_id
#         │           └── something
#         │               └── pod5
#         │                   ├── FAV*.pod5
#         │                   ├── FAV*.pod5
#         │                   ├── FAV*.pod5
#         │                   └── FAV*.pod5
#         └── Scripts
#             ├── Sample_MetaBarcodingDorado_$RUN.xlsx
#             ├── P01_Basecalling_Dorado_MetaBarcoding_Template.sh
#             ├── P02_Cleaning_Statistics_MetaBarcoding_Template.sh
#             ├── P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh
#             └── Config_Pipeline_MetaBarcoding_Nanopore_V14.sh
#
# -----------------------------------------------

# -----------------------------------------------
# Submission of the control script from Scripts/ folder
# -----------------------------------------------
# For help type:
# bash Config_Pipeline_MetaBarcoding_Nanopore_V14.sh -h
# -----------------------------------------------

# -----------------------------------------------
# Function to display help message
# -----------------------------------------------
usage() {
    echo "
 __  __        _           ____                                 _  _                  ___   _   _  _____ 
|  \/  |  ___ | |_   __ _ | __ )   __ _  _ __   ___   ___    __| |(_) _ __    __ _   / _ \ | \ | ||_   _|
| |\/| | / _ \| __| / _  ||  _ \  / _  || '__| / __| / _ \  / _  || || '_ \  / _  | | | | ||  \| |  | |  
| |  | ||  __/| |_ | (_| || |_) || (_| || |   | (__ | (_) || (_| || || | | || (_| | | |_| || |\  |  | |  
|_|  |_| \___| \__| \__,_||____/  \__,_||_|    \___| \___/  \__,_||_||_| |_| \__, |  \___/ |_| \_|  |_|  
                                                                             |___/                       
                                                                                                         

    "
    echo "Usage: bash $0 -p PROJECT -r RUN -ar ARRAY ...and options you want to change from default values"
    echo ""
    echo -e "   -p, --project     indicate the name of the project folder (e.g. CCProject) \e[31m(required)\e[0m"
    echo -e "   -r, --run         indicate the name of the folder with input reads to analyses (e.g. 20241217_CCProject_bar1-8_7ng) \e[31m(required)\e[0m"
    echo -e "  -ar, --array       indicate the number of samples in the input folder to analyses (e.g. 1-24) \e[31m(required)\e[0m"
    echo "Cleaning:"
	echo "  -qs, --qscore      indicate minimum quality Phred score for Nanopore reads [default 15]"
    echo "  -mi, --minlength   indicate minimum length of cleaned Nanopore amplicon [default 1000]"
	echo "  -ma, --maxlength   indicate maximum length of cleaned Nanopore amplicon [default 2000]"
    echo "  -br, --barcodekit  indicate name of the barcode kit used for sequencing [default SQK-NBD114-24]"
	echo "Taxonomic Profile:"
	echo "  -em, --emu         indicate if you want to perform taxonomic profile using EMU. yes [default] / no"
	echo "  -kb, --kraken      indicate if you want to perform taxonomic profile using Kraken2/Bracken. yes / no [default]"
	echo "  -bn, --blastn      indicate if you want to perform taxonomic profile using BLASTN. yes / no [default]"
	echo "Part of the pipeline to submit:"
	echo "  -pi, --pipeline    indicate part of the pipeline you want to do. Options available: all [default] - Dorado - Cleaning - Taxo - CleaningTaxo"
    echo ""
	echo "   -h, --help        show this help message"
    echo ""
    echo -e "\e[31m!! DONT FORGET TO ADD IN SCRIPTS/ YOUR Sample_MetaBarcodingONT_RUN.xlsx FILE !!\e[0m"
    echo ""
    exit 1
}
# -----------------------------------------------


# -----------------------------------------------
# SLURM EMAIL
# -----------------------------------------------
# $EMAIL: Indicate your email address for the Slurm SBATCH line --mail-user= (e.g. nicolas.nesi@unicaen.fr)
EMAIL="nicolas.nesi@unicaen.fr"
# -----------------------------------------------

# -----------------------------------------------
# SLURM ARRAY JOBS
# -----------------------------------------------
# $ARRAY: numbers to be given to sbatch argument --array= corresponding to the number of lines in your Sample_MetaBarcodingONT_$RUN.xlsx file
    # A full 96-plate: 1-96
    # Only run some samples: 1,12,35-40
ARRAY=""
# -----------------------------------------------

# -----------------------------------------------
# DIRECTORIES and INPUT FILES PATH
# -----------------------------------------------
# $GROUP: path to the folder in /dlocal/home/201913/Data/ - Options available:
    # BacterioCaen
    # Colab
    # Rouen
    # ViroCaen
    # Labeo
GROUP="BacterioCaen"
# $PROJECT: name of the Project main directory in /dlocal/home/201913/Data/$GROUP (e.g. DYNAMICH, Microbiosthme_alice, DynamiK)
PROJECT=""
# $RUN: name of Run, corresponding to the MinKNOW folder created (e.g. 20231208_DynamiKTest_bar1-8_60ng)
RUN=""
# -----------------------------------------------

# -----------------------------------------------
# DORADO BASECALLING
# -----------------------------------------------
# $MODEL: Dorado model file for basecalling. AnalyteType_PoreType_ChemistryType_TranslocationSpeed_ModelType@ModelVersionNumber	
        # dna_r10.4.1_e8.2_400bps_fast@v5.0.0 (Fast basecalling for pore type R10.4.1 FLO-MIN114)
        # dna_r10.4.1_e8.2_400bps_hac@v5.0.0 (High accuracy basecalling for pore type R10.4.1 FLO-MIN114)
        # dna_r10.4.1_e8.2_400bps_sup@v5.0.0 (Super accurate basecalling for pore type R10.4.1 FLO-MIN114) *Choose this one for best reads*
MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# $BARCODE: name of the barcode kit used (inside the double quotes) - Options available:
        # SQK-NBD114-24
        # SQK-NBD114-96
        # EXP-NBD114-24
        # EXP-NBD196
BARCODE="SQK-NBD114-24"
# $QSCOREDORADO: value for the --min_qscore parameters of Dorado basecaller (e.g. 15)
QSCOREDORADO="15"
# -----------------------------------------------

# -----------------------------------------------
# CLEANING
# -----------------------------------------------
# 16S: 1000-2000  ITS: 500-1500
# $MINLENGTH: minimum length wanted for the Nanopore cleaned reads (e.g. 1000)
MINLENGTH="1000"
# $MAXLENGTH: maximum length wanted for the Nanopore cleaned reads (e.g. 2000)
MAXLENGTH="2000"
# $QSCORE: value for the --qualified_quality_phred argument (similar to Dorado) (e.g. 10)
QSCORE="15"
# $FWD_PRIMER: Sequence of the forward primer used for amplification (e.g. 27F_16S_forward: AGAGTTTGATCCTGGCTCAG)
FWD_PRIMER="AGAGTTTGATCATGGCTCA"
# $REV_PRIMER: Sequence of the reverse primer used for amplification (e.g. 1492R_16S_reverse TACGGYTACCTTGTTACGACTT)
REV_PRIMER="TACGGTTACCTTGTTACGACTT"
# -----------------------------------------------

# -----------------------------------------------
# TAXONOMIC PROFILE
# -----------------------------------------------
# $EMU: indicate if you want to run EMU on your sample.
          # yes
          # no 
EMU="yes"
# $EMUDATABASE: name of the database to use with EMU
   # emu-01022023 - 49301 sequences representing 17555 unique bacterial and archaeal species, collected from NCBI 16S RefSeq (17/09/2020) and rrnDB v5.6
   # unitefungi-01022023 - UNITE general fasta v8.3 fungi has been pre-built for Emu v3.0+
   # uniteall-01022023 - UNITE general fasta v8.3 all eukaryotes has been pre-built for Emu v3.0+
   # silva-01022023 - SILVA v138.1 has been pre-built for Emu v3.0+ from the DADA2 SILVA species-level database
   # rdp-01022023 - RDP v11.5 with NCBI taxonomy has been pre-built for Emu v3.0+
EMUDATABASE="emu-01022023"
# $KRAKEN: indicate if you want to run Kraken2/Bracken on your sample.
          # yes
          # no 
KRAKEN="no"
# $KRAKENDATABASE: name of the database to use with Kraken2/Bracken
   # Kraken2_16SEMU - Conversion of EMU 16S database for use with Kraken2/Bracken. Include Bracken for 100,150,200,250 and 300-mers
   # Kraken2_ITS - Database builded from Fungi ITS ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/ (29/10/2024). Include Bracken for 100,150,200,250 and 300-mers
KRAKENDATABASE="Kraken2_16SEMU"
# $CONFIDENCE: confidence score value for Kraken2, threshold score [0,1] interval [(use as default 0.0(ONT) or 0.51(Illumina)] - sequence called as unclassified if score below the threshold (https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring)
CONFIDENCE="0"
# $THRESHOLDKRACKEN: number of minimum hit groups (overlaping kmers sharing a common minimizer) before declaring a sequence classified (default = 2)
THRESHOLDKRACKEN="2"
# $THRESHOLDBRACKEN: number of reads required PRIOR to abundance estimation to perform reestimation with Bracken (default = 10)
THRESHOLDBRACKEN="10"
# $READLEN: Bracken kmers length parameter (check for the presence of kmer_distrib files with the database) (default 100)
KMERLENGTH="100"
# $BLASTN: indicate if you want to run blastn on your sample. Carreful can get a very long time if a lot of data.
          # yes
          # no 
BLASTN="no"
# $BLASTDATABASE: name of the database to use for blast analyses
   # BLAST_16S_Custom - Conversion of EMU 16S database for use with BLAST+
   # 16S_ribosomal_RNA - 16S ribosomal RNA (Bacteria and Archaea type strains) 0.0179 2023-10-10T00:00:00
   # ITS_RefSeq_Fungi - Internal transcribed spacer region (ITS) from Fungi type and reference material 0.0067 2022-10-28T00:00:00
BLASTDATABASE="BLAST_16S_Custom"
#$EVALUE: Expect value (E) for saving hits (e.g. 1e-5)
EVALUE="1e-5"
# -----------------------------------------------

# -----------------------------------------------
# SELECT WHAT YOU WANT TO RUN
# -----------------------------------------------
# $PIPELINE: indidicate if run all the pipeline
   # Dorado
   # Cleaning
   # Taxo
   # CleaningTaxo
PIPELINE="all"
# -----------------------------------------------

# -----------------------------------------------
# Parse command-line arguments
# -----------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project)  PROJECT="$2"; shift 2 ;;
        -r|--run)  RUN="$2"; shift 2 ;;
        -ar|--array)  ARRAY="$2"; shift 2 ;;
        -pi|--pipeline)  PIPELINE="$2"; shift 2 ;;
        -qs|--qscore)  QSCOREONT="$2"; shift 2 ;;
        -mi|--minlength)  MINLENGTH="$2"; shift 2 ;;
        -ma|--maxlength)  MAXLENGTH="$2"; shift 2 ;;
        -br|--barcodekit)  BARCODE="$2"; shift 2 ;;
        -em|--emu)  EMU="$2"; shift 2 ;;
		-kb|--kraken)  KRAKEN="$2"; shift 2 ;;
		-bl|--blastn)  BLASTN="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) 
            echo "Invalid option: $1" >&2
            usage
			exit 1
            ;;
    esac
done
# -----------------------------------------------

# -----------------------------------------------
# Submission of the scripts
# -----------------------------------------------

# $PATH-SCRIPTS: path to scripts and Sample_HybridPipeline_$RUN.xlsx files
PATHSCRIPTS="/dlocal/home/2019013/Data/${GROUP}/${PROJECT}/Scripts"

# Ensure required variables are set
MISSING_VARS=()

for VAR in PROJECT EMAIL ARRAY GROUP RUN PATHSCRIPTS MODEL BARCODE QSCORE \
           MINLENGTH MAXLENGTH FWD_PRIMER REV_PRIMER EMU EMUDATABASE KRAKEN KRAKENDATABASE CONFIDENCE \
           KMERLENGTH THRESHOLDKRACKEN THRESHOLDBRACKEN BLASTN BLASTDATABASE EVALUE PIPELINE; do
   if [ -z "${!VAR}" ]; then
       echo -e "\e[31mError: $VAR is not set!\e[0m"
       MISSING_VARS+=("$VAR")
   fi
done

# If there are missing variables, call usage and exit
if [[ ${#MISSING_VARS[@]} -gt 0 ]]; then
    usage
    exit 1
fi

sleep 1

RUNDATE="$(date '+%Y-%m-%d %H:%M:%S')"

# Show variables indicated
echo ""
echo "Description of the analysis submitted the $RUNDATE:"
echo "Project folder: $PROJECT"
echo "Input folder: $RUN"
echo "Size array job: $ARRAY"
echo "Part of the pipeline submitted: $PIPELINE"
echo "Minimum Qscore Nanopore reads: $QSCORE"
echo "Minimum length of the Nanopore reads: $MINLENGTH"
echo "Maximum length of the Nanopore reads: $MAXLENGTH"
echo "Barcode kit used: $BARCODE"
if [[ "$EMU" == "yes" ]]; then
  echo "Taxoprofile using EMU and database $EMUDATABASE"
fi
if [[ "$KRAKEN" == "yes" ]]; then
  echo "Taxoprofile using Kraken2/Bracken and database $KRAKENDATABASE"
fi
if [[ "$BLASTN" == "yes" ]]; then
  echo "Taxoprofile using Blastn and database $BLASTDATABASE with a evalue of $EVALUE"
fi
echo ""

# Create List samples files in Scripts folder
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate NN_env

csvtk xlsx2csv "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".xlsx > "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv.temp
awk -F',' '!(($5 == "" || $5 == " ") && ($6 == "" || $6 == " ") && ($7 == "" || $7 == " "))' "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv.temp > "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv
rm "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv.temp
csvtk cut -f 1-5 "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv > "${PATHSCRIPTS}"/Sample_sheet_"$RUN".csv # For P01 Dorado demultiplexing
csvtk cut -f 5 "${PATHSCRIPTS}"/Sample_MetaBarcodingDorado_"$RUN".csv | csvtk del-header > "${PATHSCRIPTS}"/List_samples_"$RUN".txt # For P02 and P03 Array job

conda deactivate
conda deactivate

# indicate your email address 
if [[ -n "$EMAIL" ]]; then
    # Loop through all matching files
    for file in "${PATHSCRIPTS}"/P*_Template.sh; do
        if [[ -f "$file" ]]; then
            if ! grep -q "#SBATCH --mail-user=$EMAIL" "$file"; then
                if grep -q "#SBATCH --mail-user=mymail" "$file"; then
                    sed -i "s|#SBATCH --mail-user=mymail|#SBATCH --mail-user=$EMAIL|" "$file"
                    echo "Replaced 'mymail' with $EMAIL in $file"
                else
                    echo "Warning: 'mymail' not found in $file. Skipping."
                fi
            else
                echo "Email $EMAIL is already present in $file"
            fi
        else
            echo "Warning: $file does not exist or is not accessible."
        fi
    done
else
    echo "Warning: \$EMAIL is not set. Please provide an email address."
fi

if [[ "$PIPELINE" == "all" ]]; then
   # Step 1: Submission of the script P01_Basecalling_Dorado_MetaBarcoding_Template.sh
   echo "STEP 1 - Basecalling of pod5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Basecalling_Dorado_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MODEL $BARCODE $QSCORE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

   # Step 2: Submission of the script P02_Cleaning_Statistics_MetaBarcoding_Template.sh
   echo "STEP 2 - Cleaning and trimming reads with qcat and QC reads with NanoQ"
   jid_step2=$(sbatch --parsable --dependency=afterok:${jid_step1} --array=$ARRAY $PATHSCRIPTS/P02_Cleaning_Statistics_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $QSCORE $FWD_PRIMER $REV_PRIMER)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"

   # Step 3: Submission of the script P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_Template.sh
   echo "STEP 3 - Taxonomic profile analyses using EMU, Kraken2/Bracken and BLAST, extraction of reads of a specific OTU group"
   jid_step3=$(sbatch --parsable --dependency=afterany:${jid_step2} --array=$ARRAY $PATHSCRIPTS/P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMU $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

elif [[ "$PIPELINE" == "Dorado" ]]; then
   # Step 1: Submission of the script P01_Basecalling_Dorado_MetaBarcoding_Template.sh
   echo "STEP 1 - Basecalling of pod5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Basecalling_Dorado_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MODEL $BARCODE $QSCORE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

elif [[ "$PIPELINE" == "Cleaning" ]]; then
   # Step 2: Submission of the script P02_Cleaning_Statistics_MetaBarcoding_Template.sh
   echo "STEP 2 - Cleaning and trimming reads with qcat and QC reads with NanoQ"
   jid_step2=$(sbatch --parsable --array=$ARRAY $PATHSCRIPTS/P02_Cleaning_Statistics_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $QSCORE $FWD_PRIMER $REV_PRIMER)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"

elif [[ "$PIPELINE" == "Taxo" ]]; then
   # Step 3: Submission of the script P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_Template.sh
   echo "STEP 3 - Taxonomic profile analyses using EMU, Kraken2/Bracken and BLAST, extraction of reads of a specific OTU group"
   jid_step3=$(sbatch --parsable --array=$ARRAY $PATHSCRIPTS/P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMU $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

elif [[ "$PIPELINE" == "CleaningTaxo" ]]; then
   # Step 2: Submission of the script P02_Cleaning_Statistics_MetaBarcoding_Template.sh
   echo "STEP 2 - Cleaning and trimming reads with qcat and QC reads with NanoQ"
   jid_step2=$(sbatch --parsable --array=$ARRAY $PATHSCRIPTS/P02_Cleaning_Statistics_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $QSCORE $FWD_PRIMER $REV_PRIMER)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"

   # Step 3: Submission of the script P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_Template.sh
   echo "STEP 3 - Taxonomic profile analyses using EMU, Kraken2/Bracken and BLAST, extraction of reads of a specific OTU group"
   jid_step3=$(sbatch --parsable --dependency=afterany:${jid_step2} --array=$ARRAY $PATHSCRIPTS/P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMU $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)


   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

else
  echo "Invalid value for ANALYSES variable. Please set it to 'all', 'Dorado', 'Cleaning','Taxo' or 'CleaningTaxo'."
  exit 1
fi
# -----------------------------------------------
