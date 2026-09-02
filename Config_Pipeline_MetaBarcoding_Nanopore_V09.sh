#!/usr/bin/env bash

# -----------------------------------------------
# Control file to config and submit the three scripts
# of the MetaBarcoding Pipeline using Nanopore Input data
# -----------------------------------------------

# -----------------------------------------------
# Version v1.5
# -----------------------------------------------
# authors: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:20/01/2025
# -----------------------------------------------

# -----------------------------------------------#
# -----------------------------------------------#
# EDIT THE FOLLOWING LINES OF THE CONFIG SCRIPTS #
# -----------------------------------------------#
# -----------------------------------------------#

# -----------------------------------------------
# READ INSTRUCTION BEFORE SUBMISSION
# -----------------------------------------------
# Configuration file for the submission of the pipeline for MetaBarconding analyses using Nanopore amplicon with chimestry V09.

# Start point is a folder /dlocal/home/2019013/Data/$GROUP/$PROJECT with:

#   - subfolder Nanopore_outputRun/ where you will copy your Nanopore output folder of the $RUN

#   - subfolder Scripts/ with the MetaBarcoding Pipeline scripts and Sample_MetaBarcodingGuppy_$RUN.xlsx
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
#         │               └── fast5
#         │                   ├── FAV*.fast5
#         │                   ├── FAV*.fast5
#         │                   ├── FAV*.fast5
#         │                   └── FAV*.fast5
#         └── Scripts
#             ├── Sample_MetaBarcodingGuppy_$RUN.xlsx
#             ├── P01_Basecalling_Guppy_MetaBarcoding_Template.sh
#             ├── P02_Cleaning_Statistics_MetaBarcoding_Template.sh
#             ├── P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh
#             └── Config_Pipeline_MetaBarcoding_Nanopore_V09.sh
#
# -----------------------------------------------

# -----------------------------------------------
# Submission of the control script
# -----------------------------------------------
# from Scripts/ run:
# bash Config_Pipeline_MetaBarcoding_Nanopore_V09.sh
# -----------------------------------------------

# -----------------------------------------------
# SLURM EMAIL
# -----------------------------------------------
# $EMAIL: Indicate your email address for the Slurm SBATCH line --mail-user= (e.g. nicolas.nesi@unicaen.fr)
EMAIL=""
# -----------------------------------------------

# -----------------------------------------------
# SLURM ARRAY JOBS
# -----------------------------------------------
# $ARRAY: numbers to be given to sbatch argument --array= corresponding to the number of lines in your Sample_MetaBarcodingGuppy_$RUN.xlsx file !!Skipe 1 if not want no barcode!!
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
GROUP=""
# $PROJECT: name of the Project main directory in /dlocal/home/201913/Data/$GROUP (e.g. DYNAMICH, Microbiosthme_alice, DynamiK)
PROJECT=""
# $RUN: name of Run, corresponding to the MinKNOW folder created (e.g. 20231208_DynamiKTest_bar1-8_60ng)
RUN=""
# $PATH-SCRIPTS: path to scripts and Sample_MetaBarcoding_$RUN.xlsx
PATHSCRIPTS="/dlocal/home/2019013/Data/${GROUP}/${PROJECT}/Scripts"
# -----------------------------------------------

# -----------------------------------------------
# GUPPY BASECALLING
# -----------------------------------------------
# $CONFIGURATION: Guppy configuration file for basecalling.
    # dna_r9.4.1_450bps_fast.cfg (Fast basecalling for pore type R9.4.1 FLO-MIN106)
    # dna_r9.4.1_450bps_hac.cfg (High accuracy basecalling for pore type R9.4.1 FLO-MIN106)
    # dna_r9.4.1_450bps_sup.cfg (Super-accurate basecalling for pore type R9.4.1 FLO-MIN106) !!!Choose this one!!!
CONFIGURATION="dna_r9.4.1_450bps_hac.cfg"
# $BARCODE: name of the barcode kit used (inside the double quotes) - Options available:
    # EXP-NBD104
    # EXP-NBD114
    # EXP-NBD104 EXP-NBD114
    # EXP-NBD196
BARCODE="EXP-NBD196"
# $QSCORE: value for the --min_qscore parameters of guppy_basecaller (e.g. 7)
QSCORE="7"
# -----------------------------------------------

# -----------------------------------------------
# CLEANING
# -----------------------------------------------
# 16S: 1000-2000  ITS: 500-1500
# $MINLENGTH: minimum length wanted for the Nanopore cleaned reads (e.g. 1000)
MINLENGTH="1000"
# $MAXLENGTH: maximum length wanted for the Nanopore cleaned reads (e.g. 2000)
MAXLENGTH="2000"
# $FWD_PRIMER: Sequence of the forward primer used for amplification (e.g. 27F_16S_forward: AGAGTTTGATCCTGGCTCAG)
FWD_PRIMER="AGAGTTTGATCCTGGCTCAG"
# $REV_PRIMER: Sequence of the reverse primer used for amplification (e.g. 1492R_16S_reverse TACGGYTACCTTGTTACGACTT)
REV_PRIMER="TACGGYTACCTTGTTACGACTT"
# -----------------------------------------------

# -----------------------------------------------
# TAXONOMIC PROFILE
# -----------------------------------------------
# $EMUDATABASE: name of the database to use with EMU
   # emu-01022023 - 49243 sequences representing 17555 unique bacterial and archaeal species, collected from NCBI 16S RefSeq (17/09/2020) and rrnDB v5.6
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
   # 16SEMU - Conversion of EMU 16S database for use with Kraken2/Bracken. Include Bracken for 100,150,200,250 and 300-mers
   # ITS - Database builded from Fungi ITS ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/ (29/10/2024). Include Bracken for 100,150,200,250 and 300-mers
KRAKENDATABASE="16SEMU"
# $CONFIDENCE: confidence score value for Kraken2, threshold score [0,1] interval [(use as default 0.0(ONT) or 0.51(Illumina)] - sequence called as unclassified if score below the threshold (https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring)
CONFIDENCE="1"
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
# $ANALYSES: indidicate if run all the pipeline
   # all
   # guppy
   # cleaning
   # taxo
   # cleaningtaxo
ANALYSES=""
# -----------------------------------------------

# -----------------------------------------------
# Submission of the scripts
# -----------------------------------------------

# Ensure required variables are set
if [ -z "$EMAIL" ] || [ -z "$ARRAY" ] || [ -z "$GROUP" ] || [ -z "$PROJECT" ] || [ -z "$RUN" ] || [ -z "$PATHSCRIPTS" ] || \
   [ -z "$CONFIGURATION" ] || [ -z "$BARCODE" ] || [ -z "$QSCORE" ] || \
   [ -z "$MINLENGTH" ] || [ -z "$MAXLENGTH" ] || [ -z "FWD_PRIMER" ] || [ -z "$REV_PRIMER" ] \
   [ -z "$EMUDATABASE" ] || \
   [ -z "$KRAKEN" ] || [ -z "$KRAKENDATABASE" ] [ -z "$CONFIDENCE" ] || [ -z "$THRESHOLDKRACKEN" ] || [ -z "$THRESHOLDBRACKEN" ] || [ -z "$KMERLENGTH" ] || \
   [ -z "$BLASTN" ] || [ -z "$BLASTDATABASE" ] || [ -z "$EVALUE" ] || \
   [ -z "$ANALYSES" ]; then
  echo -ne "\n\033[31m\033[1mError: one or more required variables are not set. \033[0m\n\n"
  exit 1
fi

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

bash /home/2019013/nnesi02/Test/BL/BL.sh

# sbatch submissions
if [[ "$ANALYSES" == "all" ]]; then
   # Step 1: Submission of the script P01_Basecalling_Guppy_MetaBarcoding_Template.sh
   echo "STEP 1 - Basecalling of fast5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Basecalling_Guppy_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $CONFIGURATION $BARCODE $QSCORE)

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
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

elif [[ "$ANALYSES" == "guppy" ]]; then
   # Step 1: Submission of the script P01_Basecalling_Guppy_MetaBarcoding_Template.sh
   echo "STEP 1 - Basecalling of fast5 and rename of fastq"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_Basecalling_Guppy_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $CONFIGURATION $BARCODE $QSCORE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

elif [[ "$ANALYSES" == "cleaning" ]]; then
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

elif [[ "$ANALYSES" == "taxo" ]]; then
   # Step 3: Submission of the script P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_Template.sh
   echo "STEP 3 - Taxonomic profile analyses using EMU, Kraken2/Bracken and BLAST, extraction of reads of a specific OTU group"
   jid_step3=$(sbatch --parsable --array=$ARRAY $PATHSCRIPTS/P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh \
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

elif [[ "$ANALYSES" == "cleaningtaxo" ]]; then
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
   $GROUP $PROJECT $RUN $MINLENGTH $MAXLENGTH $EMUDATABASE $KRAKEN $KRAKENDATABASE $CONFIDENCE $THRESHOLDKRACKEN $THRESHOLDBRACKEN $KMERLENGTH $BLASTN $BLASTDATABASE $EVALUE)


   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

else
  echo "Invalid value for ANALYSES variable. Please set it to 'all', 'guppy', 'cleaning','taxo' or 'cleaningtaxo'."
  exit 1
fi
# -----------------------------------------------
