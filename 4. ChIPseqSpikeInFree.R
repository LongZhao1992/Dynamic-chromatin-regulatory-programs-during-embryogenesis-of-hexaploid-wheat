# Prerequisites

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicAlignments")

# Installation
# Install this package from GitHub
install.packages("devtools")
library(devtools)
install_github("stjude/ChIPseqSpikeInFree")

# Usage

library("ChIPseqSpikeInFree")
library("tidyverse")
metaFile <- "/path/sample_meta.txt"
metaData <- read_tsv("/path/sample_meta.txt")
bams <- metaData$ID
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "/path/cs.genome", metaFile = metaFile, prefix = "data1")

# Usage for differential analysis with DESeq2

library(DESeq2)
library(plyranges)
library(ChIPseeker)
library(GenomicFeatures)
txdb <- loadDb("/path/wheat.CS.txDb.sqlite")

data <- read_tsv("/path/_SF.txt")
SF <- dat$SF
countData <- read_tsv("/path/readCount.txt")
countData_mt <- countData[,-1] %>% as.matrix()
condition <- factor(c("DPA0","DPA0","DPA2","DPA2"))
dds <- DESeqDataSetFromMatrix(countData_mt, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)
coldata<- colData(dds)
coldata$sizeFactor <- coldata$sizeFactor *SF
colData(dds) <- coldata
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
result_table <- results(dds) %>% as.data.frame()

der_up <- result_table %>% filter(log2FoldChange >= 1 & padj <= 0.05)
der_down <- result_table %>% filter(log2FoldChange <= -1 & padj <= 0.05)

der_up <- der_up %>% rownames_to_column("genes")
der_down <- der_down %>% rownames_to_column("genes")

der_up_gr <- str_split(der_up$genes,"_", simplify = T) %>% as_tibble() %>% 
  mutate(start = as.integer(V2),end=as.integer(V3)) %>% bind_cols(der_up) %>% 
  dplyr::rename(seqnames="V1",start = "V2",end = "V3") %>% as_granges()

der_up_ann <- annotatePeak(der_up_gr,tssRegion=c(-3000,1000),TxDb=txdb)

der_down_gr <- str_split(der_down$genes,"_", simplify = T) %>% as_tibble() %>% 
  mutate(start = as.integer(V2),end=as.integer(V3)) %>% bind_cols(der_down) %>% 
  dplyr::rename(seqnames="V1",start = "V2",end = "V3") %>% as_granges()

der_down_ann <- annotatePeak(der_down_gr,tssRegion=c(-3000,1000),TxDb=txdb)

