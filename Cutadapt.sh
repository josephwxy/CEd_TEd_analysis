# Remain the first 200bp reads
cutadapt -j 30 -u 200 -o CE.fastq barcode1_CE.fastq
cutadapt -j 30 -u 200 -o FBL-CE.fastq barcode1_FBL-CE.fastq
cutadapt -j 30 -u 200 -o TE.fastq barcode1_TE.fastq
cutadapt -j 30 -u 200 -o FBL-TE.fastq barcode1_FBL-TE.fastq

# fastq to fasta
bioawk -c fastx '{print ">"$name; print $seq}' CE.fastq > CE.fasta
bioawk -c fastx '{print ">"$name; print $seq}' TE.fastq > TE.fasta
bioawk -c fastx '{print ">"$name; print $seq}' FBL-CE.fastq > FBL-CE.fasta
bioawk -c fastx '{print ">"$name; print $seq}' FBL-TE.fastq > FBL-TE.fasta

# blastn - 
# output format 6
#-max_hsps 1 -> remain the result with the lowest evalue）

makeblastdb -in hg19.fa -dbtype nucl -parse_seqids -out hg19.blastdb

blastn -query CE.fasta -max_hsps 1 -num_threads 8 -outfmt 6 \
-db hg19.blastdb -out CE_blastn.xls -task blastn
blastn -query TE.fasta -max_hsps 1 -num_threads 8 -outfmt 6 \
-db hg19.blastdb -out TE_blastn.xls -task blastn
blastn -query FBL-CE.fasta -max_hsps 1 -evalue 1e-6 -num_threads 8 -outfmt 6 \
-db hg19.blastdb -out FBL-CE_blastn.xls -task blastn
blastn -query FBL-TE.fasta -max_hsps 1 -evalue 1e-6 -num_threads 8 -outfmt 6 \
-db hg19.blastdb -out FBL-TE_blastn.xls -task blastn