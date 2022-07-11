# CEd_TEd_analysis

Analysis script on data analysis of gene-editing experiments in CEd and TEd manuscript.

# High-throughput Sequencing Data Analysis 
The indel analysis is basically based on CRISPREsso2: https://github.com/pinellolab/CRISPResso2

# LM-PCR and Analysis
**Input**： Fastq files 
1. cutadapt.sh: Script for preprocessing of LM-PCR raw fastq files, including: trim the reads with the first 200bp with Cutadapt, transform the fastq to fasta with Bioawk, and align the sequences against the human genome with blastn.
