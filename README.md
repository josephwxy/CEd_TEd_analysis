# CEd_TEd_analysis

Analysis script on manuscript: Transcription-coupled donor DNA expression increases homologous recombination for efficient genome editing. 
Access: https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac676/6656067

# High-throughput Sequencing Data Analysis 
The indel analysis is basically based on CRISPResso2: https://github.com/pinellolab/CRISPResso2

CRISPResso2-HTS.sh: Script for indel analysis with CRISPResso2 on samples: lamin-IsceI-CED5/lamin-IsceI-TED6/mcherry-TED/G40/G41/G50
 
# LM-PCR and Analysis
**Input**: Fastq files

**Software**: cutadapt, bioawk, R, seqkit
1. cutadapt.sh: Script for preprocessing of LM-PCR raw fastq files, including: trim the reads with the first 200bp with Cutadapt, transform the fastq to fasta with Bioawk, and align the sequences against the human genome with blastn.
2. read_annotation.R: R script on annotating the reads against the human genome. Then to index the header line of the most frequent genes.
3. seqkit.sh: The tiny script to extract the corresponding sequences from the header line.
