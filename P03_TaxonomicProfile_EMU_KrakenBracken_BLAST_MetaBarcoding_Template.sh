#!/bin/bash
#SBATCH -J P03_TaxoProfile                          # Job name
#SBATCH -o P03_TaxoProfile.%a.%A.out                # File to which stdout will be written
#SBATCH -e P03_TaxoProfile.%a.%A.err                # File to which stderr will be written
#SBATCH -p court                                    # Partition
#SBATCH -N 1                                        # Number of node - 1 GENOA node has 96-core and 750 Go RAM
#SBATCH -n 1                                        # Number of cores/cpus
#SBATCH -c 96                                       # Number of cpus ask per task
#SBATCH -t 02-00:00                                 # Runtime in DD-HH:MM
#SBATCH --mem 600G                                  # Memory for all cores in Gbytes
#SBATCH --mail-type=ALL                             # BEGIN,END,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr                          # Email address

# ---------------------------------
# Version v1.7
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, DYNAMICURE INSERM UMR 1311
# Date:03/07/2025
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
EMU="${6}"
EMUDATABASE="${7}"
KRAKEN="${8}"
KRAKENDATABASE="${9}"
CONFIDENCE="${10}"
THRESHOLDKRACKEN="${11}"
THRESHOLDBRACKEN="${12}"
KMERLENGTH="${13}"
BLASTN="${14}"
BLASTDATABASE="${15}"
EVALUE="${16}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="$FOLDER/CleanedReads_Nanopore/$RUN/NoRepeats"
PATH2="$FOLDER/TaxonomicProfile"
SCRIPTS="/home/2019013/PARTAGE/Various_Scripts_for_various_purposes"

# Databases
PATHDBEMU=/dlocal/home/2019013/Databases/EMU/"${EMUDATABASE}"
PATHDBK=/dlocal/home/2019013/Databases/"${KRAKENDATABASE}"
export BLASTDB=/dlocal/home/2019013/Databases/BLAST_NCBI/"${BLASTDATABASE}"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" $FOLDER/Scripts/List_samples_$RUN.txt)"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/"${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq.gz "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/

gunzip *.fastq.gz

echo Working directory: "${PWD}"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate emu_env
# include:
# emu=3.5.1
# ---------------------------------

# ---------------------------------
# Classify Reads using Emu
# ---------------------------------
if [[ "$EMU" == "yes" ]]; then
   echo "Running EMU with database $EMUDATABASE on sample: $SAMPLE"
   
   emu abundance \
   "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq \
   --type map-ont \
   --db "${PATHDBEMU}" \
   --keep-files \
   --keep-counts \
   --keep-read-assignments \
   --output-unclassified \
   --output-dir ./results \
   --N 50 \
   --min-abundance 0.0001 \
   --output-basename "${SAMPLE}" \
   --threads "${THREADS}" 2>&1 | tee -a "${SAMPLE}_EMU_db${EMUDATABASE}_run${RUNDATE}.log"

elif [[ "$EMU" == "no" ]]; then
    echo "Not running EMU"
else
    echo "invalid value for EMU variable. Please set it to 'yes' or 'no'"
fi
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module purge
module load py_env/anaconda3
eval "$(conda shell.bash hook)"
conda activate taxoprofile_env
# include:
# BLAST+=2.16.0
# Bracken=2.9.0
# KrakenTools=1.2
# Kraken2=2.1.3 
# KronaTools=2.8.1
# ---------------------------------

# ---------------------------------
# Classify Reads using Kraken2
# ---------------------------------
if [[ "$KRAKEN" == "yes" ]]; then
   if [[ ! -e "$PATHDBK"/hash.k2d ]]; then
     echo "Kraken2 database not found at $PATHDBK"
     exit 1
   fi

   if [[ ! -e "$PATHDBK"/database"${KMERLENGTH}"mers.kmer_distrib ]]; then
    echo "Bracken database not found at $PATHDBK"
    exit 1
   fi

   echo "Running Kraken2/Bracken with database $KRAKENDATABASE, a confidence value of $CONFIDENCE, a minimum hit group of $THRESHOLDKRACKEN and a number of read prior abundance reestimation of $THRESHOLDBRACKEN for sample: $SAMPLE"

    kraken2 \
    --db "${PATHDBK}" \
    --memory-mapping \
    --confidence "${CONFIDENCE}" \
    --threads "${THREADS}" \
    --minimum-hit-groups "${THRESHOLDKRACKEN}" \
    --output "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_Kraken2.tsv" \
    --report "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" \
    --unclassified-out "${SAMPLE}_Kraken2_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_NoHit.fastq" \
    "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq 2>&1 | tee -a "${SAMPLE}_Kraken2_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# Estimate Abundance using Bracken (number of reads originating for each OTU)
# ---------------------------------
    echo "Bracken on sample: $SAMPLE"
    
    # Species level
    bracken \
    -d "${PATHDBK}" \
    -i "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species_Bracken.tsv" \
    -w "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species.breport" \
    -r "${KMERLENGTH}" \
    -l S \
    -t "${THRESHOLDBRACKEN}" 2>&1 | tee -a "${SAMPLE}_Bracken_Species_run${RUNDATE}.log"
    
    # Genus level
    bracken \
    -d "${PATHDBK}" \
    -i "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus_Bracken.tsv" \
    -w "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus.breport" \
    -r "${KMERLENGTH}" \
    -l G \
    -t "${THRESHOLDBRACKEN}" 2>&1 | tee -a "${SAMPLE}_Bracken_Genus_run${RUNDATE}.log"
    
    # Phylum level
    bracken \
    -d "${PATHDBK}" \
    -i "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum_Bracken.tsv" \
    -w "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum.breport" \
    -r "${KMERLENGTH}" \
    -l P \
    -t "${THRESHOLDBRACKEN}" 2>&1 | tee -a "${SAMPLE}_Bracken_Phylum_run${RUNDATE}.log"
	
elif [[ "$KRAKEN" == "no" ]]; then
    echo "Not running Kraken2/Bracken"
else
    echo "invalid value for KRAKEN variable. Please set it to 'yes' or 'no'"
fi
# ---------------------------------

# ---------------------------------
# Classify Reads using BLASTN
# ---------------------------------
if [[ "$BLASTN" == "yes" ]]; then
   echo "Running Blastn against $BLASTDATABASE with evalue $EVALUE for sample: $SAMPLE"
   
   # fastq to fasta
   seqkit fq2fa "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq -o "${SAMPLE}".fasta
   
   blastn \
   -query "${SAMPLE}.fasta" \
   -db BLAST_"${BLASTDATABASE}" \
   -num_threads "${THREADS}" \
   -evalue "${EVALUE}" \
   -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms' \
   -out "${SAMPLE}_vs_${BLASTDATABASE}_blastn.tsv" 2>&1 | tee -a "${SAMPLE}_vs_${BLASTDATABASE}_blastn.log"
   
   # Get top bitscore per query
   cat "${SAMPLE}_vs_${BLASTDATABASE}_blastn.tsv" | sort -u -k1,1 --merge > "${SAMPLE}_vs_${BLASTDATABASE}_blastn_filter_topBitscore.tsv"
   awk 'BEGIN{FS=OFS="\t"; print "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tsscinames\tsskingdoms"} {print}' \
   "${SAMPLE}_vs_${BLASTDATABASE}_blastn_filter_topBitscore.tsv" > tmp_file && mv tmp_file "${SAMPLE}_vs_${BLASTDATABASE}_blastn_filter_topBitscore.tsv"
   
   # Create List of Hit per Taxa
   sed '1d' "${SAMPLE}_vs_${BLASTDATABASE}_blastn_filter_topBitscore.tsv" | cut -f 13,14,15 -d $'\t' | sort | uniq -c | sort -nrk1 > "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv"
   sed -i 's/^[ ]*//' "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv"
   sed -i 's/ /\t/' "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv"
   awk 'BEGINFILE {split(FILENAME,f,"_")} {$0=f[1]"\t"$0}1' "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv" > tmp_file && mv tmp_file "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv"
   awk 'BEGIN{FS=OFS="\t"; print "SampleID\tNbHit\tNCBITaxID\tScientificNames\tKingdoms"} {print}' "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv" > tmp_file && mv tmp_file "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv"

elif [[ "$BLASTN" == "no" ]]; then
  echo "Not running Blastn"
else
  echo "invalid value for BLAST variable. Please set it to 'yes' or 'no'"
fi
# ---------------------------------

# ---------------------------------
# Generate Krona outputs
# ---------------------------------
# EMU Krona
seqkit stats -T "${SAMPLE}"_TrimmedCleanedReads_min"${MINLENGTH}"_max"${MAXLENGTH}".fastq | csvtk cut -t -f "num_seqs" | csvtk del-header > "${SAMPLE}_readcount.txt"
python "${SCRIPTS}"/create_tsv_for_krona.py ./results/"${SAMPLE}_rel-abundance-threshold-0.0001.tsv" "${SAMPLE}_readcount.txt" "${SAMPLE}_EMU_Krona_input_0.0001.txt"
python "${SCRIPTS}"/create_tsv_for_krona.py ./results/"${SAMPLE}_rel-abundance.tsv" "${SAMPLE}_readcount.txt" "${SAMPLE}_EMU_Krona_input.txt"

ktImportText "${SAMPLE}_EMU_Krona_input_0.0001.txt" -o "${SAMPLE}_db${EMUDATABASE}_EMU_0.0001_Krona.html"
ktImportText "${SAMPLE}_EMU_Krona_input.txt" -o "${SAMPLE}_db${EMUDATABASE}_EMU_Krona.html"

if [[ "$KRAKEN" == "yes" ]]; then
    echo "Generate Krona from Kraken2, Bracken EMU and BLAST for sample: $SAMPLE"
    
    # Kraken2 Krona
    kreport2krona.py \
    -r "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_Tr${THRESHOLDKRACKEN}_Kraken2_Krona.txt"
    
    ktImportText "${SAMPLE}_db${KRAKENDATABASE}_Tr${THRESHOLDKRACKEN}_Kraken2_Krona.txt" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_Tr${THRESHOLDKRACKEN}_Kraken2_Krona.html"
    
    # Bracken Krona
    kreport2krona.py \
    -r "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species.breport" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species_Bracken_Krona.txt"
    
    ktImportText "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species_Bracken_Krona.txt" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species_Bracken_Krona.html"
    
    kreport2krona.py \
    -r "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus.breport" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus_Bracken_Krona.txt"
    
    ktImportText "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus_Bracken_Krona.txt" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus_Bracken_Krona.html"
    
    kreport2krona.py \
    -r "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum.breport" \
    -o "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum_Bracken_Krona.txt"
    
    ktImportText "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum_Bracken_Krona.txt" \
    -o "${SAMPLE}_db${TKRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum_Bracken_Krona.html"
fi

# BLAST Krona
if [[ "$BLASTN" == "yes" ]]; then
   awk -F "\t" '{print $2 "\t" $3}' "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv" > "${SAMPLE}_nbHit_per_Taxa_forKrona.tsv"
   cat "${SAMPLE}_nbHit_per_Taxa_forKrona.tsv" | taxonkit lineage -i 2 | taxonkit reformat -P -i 3 > "${SAMPLE}_nbHit_per_Taxa_TaxID_forKrona.tsv"
   cat "${SAMPLE}_nbHit_per_Taxa_TaxID_forKrona.tsv" | csvtk cut -t -f 1,4 | sed 's/;/\t/g' > "${SAMPLE}_nbHit_per_Taxa_TaxID-OK_forKrona.tsv"
   
   ktImportText "${SAMPLE}_nbHit_per_Taxa_TaxID-OK_forKrona.tsv" -o "${SAMPLE}_vs_${BLASTDATABASE}_blastn_best-hit_Krona.html"
fi
# ---------------------------------

# ---------------------------------
# Create output directories and move outputs
# ---------------------------------
# Outputs Emu
if [[ "$EMU" == "yes" ]]; then
     mkdir -p \
         "${PATH2}/EMU_${RUN}_database${TYPEDB}/Krona" \
         "${PATH2}/EMU_${RUN}_database${TYPEDB}/sam" \
         "${PATH2}/EMU_${RUN}_database${TYPEDB}/tsv" \
         "${PATH2}/EMU_${RUN}_database${TYPEDB}/NoHit"
     
     mv "${LOCAL_WORK_DIR}/results/${SAMPLE}_emu_alignments.sam" "${PATH2}/EMU_${RUN}_database${TYPEDB}/sam"
     mv "${LOCAL_WORK_DIR}/results/${SAMPLE}_unclassified_mapped.fastq" "${PATH2}/EMU_${RUN}_database${TYPEDB}/NoHit"
     mv "${LOCAL_WORK_DIR}/results/${SAMPLE}_mapped.fastq" "${PATH2}/EMU_${RUN}_database${TYPEDB}/NoHit"
     mv "${LOCAL_WORK_DIR}"/results/"${SAMPLE}"_rel-abundance* "${PATH2}/EMU_${RUN}_database${TYPEDB}/tsv"
     mv "${LOCAL_WORK_DIR}"/results/"${SAMPLE}_read-assignment-distributions.tsv" "${PATH2}/EMU_${RUN}_database${TYPEDB}/tsv"
     mv "${LOCAL_WORK_DIR}/${SAMPLE}_readcount.txt" "${PATH2}/EMU_${RUN}_database${TYPEDB}"
     mv "${SAMPLE}_EMU_Krona_input.txt" "${PATH2}/EMU_${RUN}_database${TYPEDB}/Krona"
     mv "${SAMPLE}_db${EMUDATABASE}_EMU_Krona.html" "${PATH2}/EMU_${RUN}_database${TYPEDB}/Krona"
fi

if [[ "$KRAKEN" == "yes" ]]; then
    base_dir="${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}"

    mkdir -p \
        "${base_dir}/Reports/Kraken2" \
        "${base_dir}/Reports/Bracken/Species" \
        "${base_dir}/Reports/Bracken/Genus" \
        "${base_dir}/Reports/Bracken/Phylum" \
        "${base_dir}/NoHit" \
        "${base_dir}/Krona/Kraken2" \
        "${base_dir}/Krona/Bracken"
	
   # Outputs Kraken2/Bracken
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_Kraken2.tsv" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Kraken2"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}.k2report" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Kraken2"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species_Bracken.tsv" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Species"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Species.breport" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Species"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus_Bracken.tsv" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Genus"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Genus.breport" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Genus"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum_Bracken.tsv" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Phylum"
   mv "${SAMPLE}_db${KRAKENDATABASE}_tr${THRESHOLDBRACKEN}_Phylum.breport" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Reports/Bracken/Phylum"
   mv "${SAMPLE}_Kraken2_db${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_NoHit.fastq" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/NoHit"
   mv "${SAMPLE}_db${KRAKENDATABASE}_Tr${THRESHOLDKRACKEN}_Kraken2_Krona.html" "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Krona/Kraken2"
   mv "${SAMPLE}"_db"${KRAKENDATABASE}"_tr"${THRESHOLDBRACKEN}"_*_Bracken_Krona.html "${PATH2}/Kraken2Bracken_${RUN}_database${KRAKENDATABASE}_tr${THRESHOLDKRACKEN}_tr${THRESHOLDBRACKEN}/Krona/Bracken"
fi

if [[ "$BLASTN" == "yes" ]]; then
    mkdir -p "${PATH2}/BLAST_${RUN}_database${BLASTDATABASE}/Krona"
		
	# Outputs BLAST
    mv "${SAMPLE}_vs_${BLASTDATABASE}_blastn.tsv" "${PATH2}/BLAST_${RUN}_database${BLASTDATABASE}"
    mv "${SAMPLE}_vs_${BLASTDATABASE}_blastn_filter_topBitscore.tsv" "${PATH2}/BLAST_${RUN}_database${BLASTDATABASE}"
    mv "${SAMPLE}_vs_${BLASTDATABASE}_List_nbHit_per_Taxa.tsv" "${PATH2}/BLAST_${RUN}_database${BLASTDATABASE}"
    mv "${SAMPLE}_vs_${BLASTDATABASE}_blastn_best-hit_Krona.html" "${PATH2}/BLAST_${RUN}_database${BLASTDATABASE}/Krona"
fi

if [[ ! -d "${FOLDER}/Log/${RUN}_run${RUNDATE}" ]]; then
   mkdir -p "${FOLDER}/Log/${RUN}_run${RUNDATE}"
fi
mv *.log "${FOLDER}/Log/${RUN}_run${RUNDATE}"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P03_TaxoProfile."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}/Scripts/Sbatch_log"
rm -f "${FOLDER}"/Scripts/P01_DoradoGPU.*.out
rm -f "${FOLDER}"/Scripts/P02_CleanStats.*.out
# ---------------------------------
