# MetaBarcoding Nanopore Data

Pipeline for MetaBarcoding analyses of 16S or ITS using Nanopore long reads amplicon data in HPC Austral (group 2019013)

## Configuration file ##
The configuration file is needed for the MetaBarcoding pipeline to function properly, where the following settings can be modified and customized to run the entire pipeline in a single command. Preferences of analyses can be modified to your projects specifications.
Select the Configuration file based on the Nanopore chemistry used to generate your sequencing data:
- ``Config_Pipeline_MetaBarcoding_Nanopore_Template_v09.sh`` for chemistry 109 
- ``Config_Pipeline_MetaBarcoding_Nanopore_Template_v14.sh`` for chemistry 114.

## Create Sample_MetaBarcodingGuppy_$RUN.xlsx for V09 ##
Save a ```Sample_MetaBarcodingGuppy_$RUN.xlsx``` Excel file with name of your sample for the corresponding barcode used in the sequencing $RUN. Keep the first line ```barcode00 Nobarcode```! Name must start with [a-zA-Z] have NO space and NO underscore:
```
S23-MTX-COV-NOSE
S59-MTX-HRV-LUNG
```

## Create SampleMetaBarcodingDorado_$RUN.xlsx for V14 ##
Save a ```SampleMetaBarcodingDorado_$RUN.xlsx``` Excel file with name of your sample for the corresponding barcode used in the sequencing $RUN. For barcode not used indicate NoSample1, NoSample2 etc... Name must start with [a-zA-Z] have NO space and NO underscore:
```
S23-MTX-COV-NOSE
S59-MTX-HRV-LUNG
```

## Directory structure ##
The beginning directory structure for chemistry 109 (FLOMIN-106D) should be:
```
 .
 └── GROUP/
     └── PROJECT/
         ├── Nanopore_outputRun/
         │   └── RUN/
         │       └── no_sample_id/
         │           └── something/
         │               └── fast5/
         │                   ├── FAV*.fast5
         │                   ├── FAV*.fast5
         │                   ├── FAV*.fast5
         │                   └── FAV*.fast5
         └── Scripts/
             ├── Sample_MetaBarcodingGuppy_$RUN.xlsx
             ├── P01_Basecalling_Guppy_MetaBarcoding_Template.sh
             ├── P02_Cleaning_Statistics_MetaBarcoding_Template.sh
             ├── P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh
             └── Config_Pipeline_MetaBarcoding_Nanopore_V09.sh

*/ denotes directories
```
The beginning directory structure for chemistry 114 (FLOMIN-114) should be:
```
 .
 └── GROUP/
     └── PROJECT/
         ├── Nanopore_outputRun/
         │   └── RUN/
         │       └── no_sample_id/
         │           └── something/
         │               └── pod5_skip/
         │                   ├── FAV*.pod5
         │                   ├── FAV*.pod5
         │                   ├── FAV*.pod5
         │                   └── FAV*.pod5
         └── Scripts/
             ├── Sample_MetaBarcodingDorado_$RUN.xlsx
             ├── P01_Basecalling_Dorado_MetaBarcoding_Template.sh
             ├── P02_Cleaning_Statistics_MetaBarcoding_Template.sh
             ├── P03_TaxonomicProfile_EMU_KrakenBracken_BLAST_MetaBarcoding_Template.sh
             └── Config_Pipeline_MetaBarcoding_Nanopore_V14.sh

*/ denotes directories
```

## Taxonomic Profile performed with
- [EMU] (https://github.com/treangenlab/emu)
- [Kraken2/Bracken] (https://github.com/DerrickWood/kraken2) (https://github.com/jenniferlu717/Bracken)
- [BLAST] (https://blast.ncbi.nlm.nih.gov/Blast.cgi)

## Pipeline workflow ##

![Figure_Pipeline_MetaBarcoding_Nanopore](https://github.com/user-attachments/assets/fedc78d7-0f0b-41dd-be07-238e781890e3)


## To do ##

- Generate high-accuracy consensus sequences and remove chimera before taxo profile [VSEARCH / MEDAKA] see Ohta et al. 2023
- Add in config file the option to run P2 and P3 scripts with Illumina inputs amplicon reads

## Author ##

Nicolas Nesi — UMR 1311 Inserm DYNAMICURE, Université de Caen Normandie
