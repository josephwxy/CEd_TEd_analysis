# CEd_TEd_analysis

Analysis script on data analysis of gene-editing experiments in CEd and TEd manuscript.

# High-throughput Sequencing Data Analysis 
The indel analysis is basically based on CRISPResso2: https://github.com/pinellolab/CRISPResso2

CRISPResso2-HTS.sh： Script for indel analysis with CRISPResso2 on sample：lamin-IsceI-CED5/lamin-IsceI-TED6/mcherry-TED/G40/G41/G50
 
# LM-PCR and Analysis
**Input**： Fastq files

**software**： cutadapt, bioawk, R, seqkit
1. cutadapt.sh: Script for preprocessing of LM-PCR raw fastq files, including: trim the reads with the first 200bp with Cutadapt, transform the fastq to fasta with Bioawk, and align the sequences against the human genome with blastn.
2. read_annotation.R: R script on annotating the reads against the human genome. The next section is to index the header line of the most frequent genes.
3. seqkit.sh: The tiny script to extract the corresponding sequences from the header line.
