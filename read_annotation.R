library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(AnnotationDbi)
library(ChIPpeakAnno)
library(dplyr)
library(plyr)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotate the reads with their corresponding genes, and statistics
# CE
CE <- read.table("CE_blastn.xls")
head(CE)
CE2 <- CE[CE$V4>90,]
CE2 <- CE2[CE2$V3 == 100,]
CE2 <- CE2[nchar(CE2$V2)<6,]
CE_grange <- GRanges(seqnames = CE2$V2, ranges = IRanges(start = CE2$V9,end = CE2$V9+1))
CE_anno <- annotatePeakInBatch(CE_grange,genes(txdb))
mcols(CE_anno)$gene <- AnnotationDbi::select(org.Hs.eg.db, keys = CE_anno$feature,
                                             keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))$SYMBOL
CE_anno <- as.data.frame(CE_anno)
CE_anno <- CE_anno[,c(1,2,15)]
CE_count <- plyr::count(CE_anno,names(CE_anno))
CE_list <- split(CE_count, CE_count$gene)
CE_list <- lapply(CE_list, function(x) arrange(x,-freq)[1,])
CE_table <- NULL
for (i in 1:length(CE_list)){
  CE_table <- rbind(CE_table,CE_list[[i]])
}
CE_table <- arrange(CE_table,-freq)
#write.csv(CE_table,"CE_genetable1.csv",row.names = F)
write.csv(arrange(CE_gene,-Freq),"CE_gene_frequency.csv") 

# TE
TE <- read.table("TE_blastn.xls")
head(TE)
TE2 <- TE[TE$V4>90,]
TE2 <- TE2[TE2$V3 == 100,]
TE2 <- TE2[nchar(TE2$V2)<6,]
TE2 <- GRanges(seqnames = TE2$V2, ranges = IRanges(start = TE2$V9,end = TE2$V9+1))
TE_anno <- annotatePeakInBatch(TE2,genes(txdb))
mcols(TE_anno)$gene <- AnnotationDbi::select(org.Hs.eg.db, keys = TE_anno$feature,
                                             keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))$SYMBOL
TE_anno <- as.data.frame(TE_anno)
TE_anno <- TE_anno[,c(1,2,15)]
TE_count <- plyr::count(TE_anno,names(TE_anno))
TE_list <- split(TE_count, TE_count$gene)
TE_list <- lapply(TE_list, function(x) arrange(x,-freq)[1,])
TE_table <- NULL
for (i in 1:length(TE_list)){
  TE_table <- rbind(TE_table,TE_list[[i]])
}
TE_table <- arrange(TE_table,-freq)
#write.csv(TE_table,"TE_genetable1.csv",row.names = F)
write.csv(arrange(TE_gene,-Freq),"TE_gene_frequency.csv") 

# FBL-CE
FBL_CE <- read.table("FBL-CE_blastn.xls")
head(FBL_CE)
FBL_CE2 <- FBL_CE[FBL_CE$V4>90,]
FBL_CE2 <- FBL_CE2[FBL_CE2$V3 == 100,]
FBL_CE2 <- FBL_CE2[nchar(FBL_CE2$V2)<6,]
FBL_CE2 <- GRanges(seqnames = FBL_CE2$V2, ranges = IRanges(start = FBL_CE2$V9,end = FBL_CE2$V9+1))
FBL_CE_anno <- annotatePeakInBatch(FBL_CE2,genes(txdb))
mcols(FBL_CE_anno)$gene <- AnnotationDbi::select(org.Hs.eg.db, keys = FBL_CE_anno$feature,
                                                 keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))$SYMBOL
FBL_CE_anno <- as.data.frame(FBL_CE_anno)
FBL_CE_anno <- FBL_CE_anno[,c(1,2,15)]
FBL_CE_count <- plyr::count(FBL_CE_anno,names(FBL_CE_anno))
FBL_CE_list <- split(FBL_CE_count, FBL_CE_count$gene)
FBL_CE_list <- lapply(FBL_CE_list, function(x) arrange(x,-freq)[1,])
FBL_CE_table <- NULL
for (i in 1:length(FBL_CE_list)){
  FBL_CE_table <- rbind(FBL_CE_table,FBL_CE_list[[i]])
}
FBL_CE_table <- arrange(FBL_CE_table,-freq)
#write.csv(FBL_CE_table,"FBL_CE_genetable1.csv",row.names = F)
write.csv(arrange(FBL_CE_gene,-Freq),"FBL_CE_gene_frequency.csv") 

# FBL-TE
FBL_TE <- read.table("FBL-TE_blastn.xls")
head(FBL_TE)
FBL_TE2 <- FBL_TE[FBL_TE$V4>90,]
FBL_TE2 <- FBL_TE2[FBL_TE2$V3 == 100,]
FBL_TE2 <- FBL_TE2[nchar(FBL_TE2$V2)<6,]
FBL_TE2 <- GRanges(seqnames = FBL_TE2$V2, ranges = IRanges(start = FBL_TE2$V9,end = FBL_TE2$V9+1))
FBL_TE_anno <- annotatePeakInBatch(FBL_TE2,genes(txdb))
mcols(FBL_TE_anno)$gene <- AnnotationDbi::select(org.Hs.eg.db, keys = FBL_TE_anno$feature,
                                                 keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))$SYMBOL
FBL_TE_anno <- as.data.frame(FBL_TE_anno)
FBL_TE_anno <- FBL_TE_anno[,c(1,2,15)]
FBL_TE_count <- plyr::count(FBL_TE_anno,names(FBL_TE_anno))
FBL_TE_list <- split(FBL_TE_count, FBL_TE_count$gene)
FBL_TE_list <- lapply(FBL_TE_list, function(x) arrange(x,-freq)[1,])
FBL_TE_table <- NULL
for (i in 1:length(FBL_TE_list)){
  FBL_TE_table <- rbind(FBL_TE_table,FBL_TE_list[[i]])
}
FBL_TE_table <- arrange(FBL_TE_table,-freq)
#write.csv(FBL_TE_table,"FBL_TE_genetable1.csv",row.names = F)
write.csv(arrange(FBL_TE_gene,-Freq),"FBL_TE_gene_frequency.csv") 

# Index the header line of the most frequent genes for furthering sequence extraction
CE_start <- CE_table[1:13,]$start
CE_ID <- NULL
for (i in CE_start){
  CE_ID <- rbind(CE_ID, CE2[which(CE2$V9 == i),]$V1[1])
}
write.table(CE_ID, 'CE_ID.txt',col.names = F,row.names = F,quote = F)

TE_start <- TE_table[1:14,]$start
TE_ID <- NULL
for (i in TE_start){
  TE_ID <- rbind(TE_ID, TE2[which(TE2$V9 == i),]$V1[1])
}
write.table(TE_ID, 'TE_ID.txt',col.names = F,row.names = F,quote = F)

FBL_CE_start <- FBL_CE_table[1:15,]$start
FBL_CE_ID <- NULL
for (i in FBL_CE_start){
  FBL_CE_ID <- rbind(FBL_CE_ID, FBL_CE2[which(FBL_CE2$V9 == i),]$V1[1])
}
write.table(FBL_CE_ID, 'FBL_CE_ID.txt',col.names = F,row.names = F,quote = F)

FBL_TE_start <- FBL_TE_table[1:15,]$start
FBL_TE_ID <- NULL
for (i in FBL_TE_start){
  FBL_TE_ID <- rbind(FBL_TE_ID, FBL_TE2[which(FBL_TE2$V9 == i),]$V1[1])
}
write.table(FBL_TE_ID, 'FBL_TE_ID.txt',col.names = F,row.names = F,quote = F)
