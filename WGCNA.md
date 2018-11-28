# Building co-expression network using WGCNA
BioProject PRJEB25639 contained 69 samples with three biological replicates and one sample with two biological replicates.
We use these data to constructed co-expression network.
Firstly, we got gene level counts and then got normolized counts using DESeq2.

```R
# working dir is in PRJEB25639
library("DESeq2")
library(tximportData)
library(tximport)
library(readr)
samples <- read.table(file.path("./salmon_quant", "sample_information.txt"), header = TRUE, sep = '\t')
samples
files <- file.path("./salmon_quant", samples$File, "quant.sf")
files
names(files) <- paste0(samples$sample)
samples$condition <- factor(samples$Order)
all(file.exists(files))
t2g <- read.table(file.path("/data2/user_data/rna_seq", "gene_transcript_list.txt"), header = TRUE)
t2g <- dplyr::rename(t2g,target_id = transcript_id, gene_id = gene_id)
head(t2g)
txi <- tximport(files, type = "salmon", tx2gene = t2g)
dds <- DESeqDataSetFromTximport(txi, colData = samples,design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
nor_counts <- counts(dds, normalized=TRUE)
write.table(as.data.frame(nor_counts), file="./salmon_quant/BCS_normalized_count.tsv")
```