library("ggplot2")
library("RColorBrewer")
library("factoextra")
library("DESeq2")
library("dplyr")
library("edgeR")
library("pheatmap")
library("fgsea")
library("GSA")
library("dplyr")
library("GO.db")
library("tidyr")
library("msigdbr")
library("Hmisc")
library("psych")
library("corrplot")
library("ggcorrplot")

# load iNeuron pool counttables 3 flow cells combined mapped to hg38p12_rn6 annotation file #

metadata <- read.table(".../0.Data/RNA_seq/counts_MCS/190820_EnsemblID_genes_GRCh38p12_ensembl94_metadata_biomart.txt", header=TRUE, sep="\t")
colnames(metadata) <- c("Gene_stable_ID", "Gene_stable_ID_version", "Gene_name", "Chromosome", "Gene_type")

pool1 <- read.table(".../0.Data/RNA_seq/counts_MCS/combined_tn5-KdVs-run1-KdVs-iNeurons-pool1-16967.mapped_unique.htseq_counts.nonunique_all.txt", header=TRUE, row.names=1)
colnames(pool1) <- c("run1_pool1_empty", "run1_pool1_CRISPR_2309_B3", "run1_pool1_409B2_2311_C1", "run1_pool1_70_2311_B4", "run1_pool1_71_2309_B5", "run1_pool1_WTC_2311_A2", "run1_pool1_76_2311_D4", "run1_pool1_77_2309_C1", 
                     "run2_pool1_empty", "run2_pool1_CRISPR_2309_B3", "run2_pool1_409B2_2311_C1", "run2_pool1_70_2311_B4", "run2_pool1_71_2309_B5", "run2_pool1_WTC_2311_A2", "run2_pool1_76_2311_D4", "run2_pool1_77_2309_C1", 
                     "run3_pool1_empty", "run3_pool1_CRISPR_2309_B3", "run3_pool1_409B2_2311_C1", "run3_pool1_70_2311_B4", "run3_pool1_71_2309_B5", "run3_pool1_WTC_2311_A2", "run3_pool1_76_2311_D4", "run3_pool1_77_2309_C1",
                     "run4_pool1_empty", "run4_pool1_CRISPR_2309_B3", "run4_pool1_409B2_2311_C1", "run4_pool1_70_2311_B4", "run4_pool1_71_2309_B5", "run4_pool1_WTC_2311_A2", "run4_pool1_76_2311_D4", "run4_pool1_77_2309_C1")

pool3 <- read.table(".../0.Data/RNA_seq/counts_MCS/combined_tn5-KdVs-run1-KdVs-iNeurons-pool3-16965.mapped_unique.htseq_counts.nonunique_all.txt", header=TRUE, row.names=1)
colnames(pool3) <- c("run1_pool3_71_2309_B6", "run1_pool3_409B2_2311_D2", "run1_pool3_77_2309_D1", "run1_pool3_WTC_2311_B1", "run1_pool3_empty", "run1_pool3_76_2311_D5", "run1_pool3_CRISPR_2309_A4", "run1_pool3_76_2311_C4",  
                     "run2_pool3_71_2309_B6", "run2_pool3_409B2_2311_D2", "run2_pool3_77_2309_D1", "run2_pool3_WTC_2311_B1", "run2_pool3_empty", "run2_pool3_76_2311_D5", "run2_pool3_CRISPR_2309_A4", "run2_pool3_76_2311_C4", 
                     "run3_pool3_71_2309_B6", "run3_pool3_409B2_2311_D2", "run3_pool3_77_2309_D1", "run3_pool3_WTC_2311_B1", "run3_pool3_empty", "run3_pool3_76_2311_D5", "run3_pool3_CRISPR_2309_A4", "run3_pool3_76_2311_C4",
                     "run4_pool3_71_2309_B6", "run4_pool3_409B2_2311_D2", "run4_pool3_77_2309_D1", "run4_pool3_WTC_2311_B1", "run4_pool3_empty", "run4_pool3_76_2311_D5", "run4_pool3_CRISPR_2309_A4", "run4_pool3_76_2311_C4")

pool4 <- read.table(".../0.Data/RNA_seq/counts_MCS/combined_tn5-KdVs-run1-KdVs-iNeurons-pool4-16964.mapped_unique.htseq_counts.nonunique_all.txt", header=TRUE, row.names=1)
colnames(pool4) <- c("run1_pool4_empty", "run1_pool4_70_2311_A4", "run1_pool4_WTC_2311_B2", "run1_pool4_70_2311_A5", "run1_pool4_71_2309_A5", "run1_pool4_77_2309_C2", "run1_pool4_CRISPR_2309_B4", "run1_pool4_409B2_2311_D1",  
                     "run2_pool4_empty", "run2_pool4_70_2311_A4", "run2_pool4_WTC_2311_B2", "run2_pool4_70_2311_A5", "run2_pool4_71_2309_A5", "run2_pool4_77_2309_C2", "run2_pool4_CRISPR_2309_B4", "run2_pool4_409B2_2311_D1", 
                     "run3_pool4_empty", "run3_pool4_70_2311_A4", "run3_pool3_WTC_2311_B2", "run3_pool4_70_2311_A5", "run3_pool4_71_2309_A5", "run3_pool4_77_2309_C2", "run3_pool4_CRISPR_2309_B4", "run3_pool4_409B2_2311_D1",
                     "run4_pool4_empty", "run4_pool4_70_2311_A4", "run4_pool3_WTC_2311_B2", "run4_pool4_70_2311_A5", "run4_pool4_71_2309_A5", "run4_pool4_77_2309_C2", "run4_pool4_CRISPR_2309_B4", "run4_pool4_409B2_2311_D1")


# combine replicates measured on multiple flow cells #

pool1_comb <- cbind.data.frame(rowSums(cbind(pool1[,1], pool1[,9], pool1[,17], pool1[,25])), rowSums(cbind(pool1[,2], pool1[,10], pool1[,18], pool1[,26])), rowSums(cbind(pool1[,3], pool1[,11], pool1[,19], pool1[,27])), rowSums(cbind(pool1[,4], pool1[,12], pool1[,20], pool1[,28])), rowSums(cbind(pool1[,5], pool1[,13], pool1[,21], pool1[,29])), rowSums(cbind(pool1[,6], pool1[,14], pool1[,22], pool1[,30])), rowSums(cbind(pool1[,7], pool1[,15], pool1[,23], pool1[,31])), rowSums(cbind(pool1[,8], pool1[,16], pool1[,24], pool1[,32])))
colnames(pool1_comb) <- c("pool1_empty", "pool1_CRISPR_2309_B3", "pool1_409B2_2311_C1", "pool1_70_2311_B4", "pool1_71_2309_B5", "pool1_WTC_2311_A2", "pool1_76_2311_D4", "pool1_77_2309_C1")

pool3_comb <- cbind.data.frame(rowSums(cbind(pool3[,1], pool3[,9], pool3[,17], pool3[,25])), rowSums(cbind(pool3[,2], pool3[,10], pool3[,18], pool3[,26])), rowSums(cbind(pool3[,3], pool3[,11], pool3[,19], pool3[,27])), rowSums(cbind(pool3[,4], pool3[,12], pool3[,20], pool3[,28])), rowSums(cbind(pool3[,5], pool3[,13], pool3[,21], pool3[,29])), rowSums(cbind(pool3[,6], pool3[,14], pool3[,22], pool3[,30])), rowSums(cbind(pool3[,7], pool3[,15], pool3[,23], pool3[,31])), rowSums(cbind(pool3[,8], pool3[,16], pool3[,24], pool3[,32])))
colnames(pool3_comb) <- c("pool3_71_2309_B6", "pool3_409B2_2311_D2", "pool3_77_2309_D1", "pool3_WTC_2311_B1", "pool3_empty", "pool3_76_2311_D5", "pool3_CRISPR_2309_A4", "pool3_76_2311_C4")

pool4_comb <- cbind.data.frame(rowSums(cbind(pool4[,1], pool4[,9], pool4[,17], pool4[,25])), rowSums(cbind(pool4[,2], pool4[,10], pool4[,18], pool4[,26])), rowSums(cbind(pool4[,3], pool4[,11], pool4[,19], pool4[,27])), rowSums(cbind(pool4[,4], pool4[,12], pool4[,20], pool4[,28])), rowSums(cbind(pool4[,5], pool4[,13], pool4[,21], pool4[,29])), rowSums(cbind(pool4[,6], pool4[,14], pool4[,22], pool4[,30])), rowSums(cbind(pool4[,7], pool4[,15], pool4[,23], pool4[,31])), rowSums(cbind(pool4[,8], pool4[,16], pool4[,24], pool4[,32])))
colnames(pool4_comb) <- c("pool4_empty", "pool4_70_2311_A4", "pool4_WTC_2311_B2", "pool4_70_2311_A5", "pool4_71_2309_A5", "pool4_77_2309_C2", "pool4_CRISPR_2309_B4", "pool4_409B2_2311_D1")


# select genes of interest

all <- cbind.data.frame(pool1_comb, pool3_comb, pool4_comb)
rownames(all) <- rownames(pool1)

all_hgenes <- all[grepl("ENSG", rownames(all)),] #  58721 genes

iNeuron_hgenes <- cbind.data.frame(all_hgenes[,2:12], all_hgenes[,14:16], all_hgenes[,18:24]) 

iNeuron_hgenes2 <- cbind.data.frame(Gene_stable_ID_version = rownames(iNeuron_hgenes), iNeuron_hgenes)
iNeuron_hgenes2 <- merge(metadata, iNeuron_hgenes2, by="Gene_stable_ID_version")

iNeuron_hgenes3 <- cbind.data.frame(iNeuron_hgenes2[,1:5], iNeuron_hgenes2[,grepl("WTC", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("409B2", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("77", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("CRISPR", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("71", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("76", colnames(iNeuron_hgenes2))], iNeuron_hgenes2[,grepl("70", colnames(iNeuron_hgenes2))])

# write.table(iNeuron_hgenes3, "200109_rawcounts_allgenes.txt")


#### cpm data transformation ####

cpm <- cpm(iNeuron_hgenes)

cpm2 <- cbind.data.frame(Gene_stable_ID_version = rownames(cpm), cpm)
cpm2 <- merge(metadata, cpm2, by="Gene_stable_ID_version")
cpm3 <- cbind.data.frame(cpm2[,1:5], cpm2[,grepl("WTC", colnames(cpm2))], cpm2[,grepl("409B2", colnames(cpm2))], cpm2[,grepl("77", colnames(cpm2))], cpm2[,grepl("CRISPR", colnames(cpm2))], cpm2[,grepl("71", colnames(cpm2))], cpm2[,grepl("76", colnames(cpm2))], cpm2[,grepl("70", colnames(cpm2))])


cpm_pc <- subset(cpm2, cpm2$Gene_type=="protein_coding") # 19922


# select rows with cpm > 2 in at least 3 columns #

cpm_thr2_TF <- cpm > 2
cpm_thr2_select <- cpm_thr2_TF[rowSums(cpm_thr2_TF) > 2, ] # 17601

cpm_thr2_select2 <- cbind.data.frame(Gene_stable_ID_version = rownames(cpm_thr2_select), cpm_thr2_select)
cpm_thr2 <- cpm2[cpm2$Gene_stable_ID_version %in% cpm_thr2_select2$Gene_stable_ID_version,]
# 17601 genes left, of which 13721 protein_coding

cpm_thr3 <- cbind.data.frame(cpm_thr2[,1:5], cpm_thr2[,grepl("WTC", colnames(cpm_thr2))], cpm_thr2[,grepl("409B2", colnames(cpm_thr2))], cpm_thr2[,grepl("77", colnames(cpm_thr2))], cpm_thr2[,grepl("CRISPR", colnames(cpm_thr2))], cpm_thr2[,grepl("71", colnames(cpm_thr2))], cpm_thr2[,grepl("76", colnames(cpm_thr2))], cpm_thr2[,grepl("70", colnames(cpm_thr2))])

# write.table(cpm_thr3, "200109_cpm_above2_in_3repl.txt")


### subset counts table for genes with cpm > 2 in 3 repl ###

counts_cpm2 <- iNeuron_hgenes2[iNeuron_hgenes2$Gene_stable_ID_version %in% cpm_thr2$Gene_stable_ID_version,]
counts_final <- cbind.data.frame(counts_cpm2[,1:5], counts_cpm2[,grepl("WTC", colnames(counts_cpm2))], counts_cpm2[,grepl("409B2", colnames(counts_cpm2))], counts_cpm2[,grepl("77", colnames(counts_cpm2))], counts_cpm2[,grepl("CRISPR", colnames(counts_cpm2))], counts_cpm2[,grepl("71", colnames(counts_cpm2))], counts_cpm2[,grepl("76", colnames(counts_cpm2))], counts_cpm2[,grepl("70", colnames(counts_cpm2))])


#### DE analysis ####

# remove Y genes from count table to correct for gender effect #

counts_gencor <- subset(counts_final, counts_final$Chromosome != "Y")  # 17.584 genes (17 Y genes removed)

# make DGE list #

groups <- as.factor(c(paste(rep("CTL",9)), paste(rep("KdVS", 12))))
iNeurons_DGE <- DGEList(counts=counts_gencor[,6:26], genes=counts_gencor[,1:5], group=groups)


# metadata regarding samples #

samplenames <- colnames(iNeurons_DGE$counts)
condition <- as.factor(c(paste(rep("control", 9)), paste(rep("patient", 12))))
cell_line <- as.factor(c(paste(rep("C1", 3)), paste(rep("C2", 3)), paste(rep("C3",3)), paste(rep("C3.CR",3)), paste(rep("P1", 3)), paste(rep("P2",3)), paste(rep("P3",3))))
cell_line <- factor(cell_line, levels=c("C1", "C2","C3", "C3.CR", "P1", "P2", "P3"))
gender <- c(paste(rep("male", 3)), paste(rep("female", 18)))
MEA <-  c(paste(rep("2311", 6)), paste(rep("2309", 9)), paste(rep("2311", 6)))
purdate <-  c("set1", "set3", "set2", "set2", "set2", "set1", "set1", "set1", "set2", "set3", "set2", "set1", "set3", "set3", "set2", "set1", "set3", "set3", "set2", "set3", "set1")
pool <- c(paste(rep(c("pool1", "pool3", "pool4"), 5)), "pool1", "pool3", "pool3" , "pool1", "pool4", "pool4")
age <- c(paste(rep("30", 3)), paste(rep("36", 3)), paste(rep("41", 6)), paste(rep("7",3)), paste(rep("10",3)), paste(rep("21",3)))
age <- factor(age, levels=c("7", "10", "21", "30", "36", "41"))

iNeurons_DGE$samples$samplenames <- samplenames
iNeurons_DGE$samples$condition <- condition
iNeurons_DGE$samples$cell_line <- cell_line
iNeurons_DGE$samples$gender <- gender
iNeurons_DGE$samples$MEA <- MEA
iNeurons_DGE$samples$purdate <- purdate
iNeurons_DGE$samples$pool <- pool
iNeurons_DGE$samples$age <- age

iNeurons_DGE <- calcNormFactors(iNeurons_DGE, method="TMM")


### VOOM transform data ###

design <- model.matrix(~0+condition+MEA)
voom <- voom(iNeurons_DGE, design, save.plot=TRUE)


# heatmap(cor(voom$E)) # not corrected for batch


### Differential expression analysis - all lines ###

blocks <- c(paste(rep("1", 3)), paste(rep("2",3)), paste(rep("3",3)), paste(rep("4",3)), paste(rep("5", 3)), paste(rep("6",3)), paste(rep("7",3)))
corfit <- duplicateCorrelation(voom, design, block=blocks) # corfit$consensus = 0.28
fit <- lmFit(voom, design, block=blocks, correlation = corfit$consensus)

contr <- makeContrasts(conditionpatient - conditioncontrol, levels=colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contr)
fit3 <- eBayes(fit2)

DE_res <- topTable(fit3, sort.by="P", adjust="BH", number=Inf) # 20 DE genes, patient vs control
DE_genes <- subset(DE_res, DE_res$adj.P.Val < 0.05)
DE_genes <- DE_genes[order(DE_genes$adj.P.Val),]

# correct voom data for batch for data visualization #

voom_cor <- removeBatchEffect(voom, batch=iNeurons_DGE$samples$MEA)
voom_cor2 <- cbind.data.frame(voom$genes, voom_cor)

# pheatmap(cor(voom_cor))

DE_voom <- merge(DE_res, voom_cor2)
DE_genes_voom <- merge(DE_genes, voom_cor2)
DE_genes_voom <- DE_genes_voom[order(DE_genes_voom$adj.P.Val),]
rownames(DE_genes_voom) <- DE_genes_voom$Gene_name


# scale per gene #

DE_genes_voom2 <- as.data.frame(t(DE_genes_voom[,12:32]))

gene_scale <- scale(DE_genes_voom2[,1])
colnames(gene_scale) <- colnames(DE_genes_voom2)[1]
rownames(gene_scale) <- rownames(DE_genes_voom2)
DE_genes_voom_scgen <- data.frame(gene_scale)
for (i in 2:20){
  gene_scale <- scale(DE_genes_voom2[,i])
  gene_scale2 <- data.frame(gene_scale)
  colnames(gene_scale2) <- colnames(DE_genes_voom2)[i]
  DE_genes_voom_scgen <- cbind(DE_genes_voom_scgen, gene_scale2)
}
DE_genes_voom_scgen <- t(DE_genes_voom_scgen)

DE_genes_voom_scgen2 <- DE_genes_voom_scgen
colnames(DE_genes_voom_scgen2) <- cell_line

breaksList = seq(-2.5, 2, by = 0.5)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
pheatmap(DE_genes_voom_scgen2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames=TRUE, annotation_col=NA, col=colors, 
         breaks=breaksList, annotation_names_col=FALSE, main="DE genes in iNeurons", fontsize=18)


DE_voom2 <- DE_voom[order(DE_voom$adj.P.Val),]
DE_voom2 <- as.data.frame(t(DE_voom2[,12:32]))

gene_scale <- scale(DE_voom2[,1])
colnames(gene_scale) <- colnames(DE_voom2)[1]
rownames(gene_scale) <- rownames(DE_voom2)
DE_voom_scgen <- data.frame(gene_scale)
for (i in 2:200){
  gene_scale <- scale(DE_voom2[,i])
  gene_scale2 <- data.frame(gene_scale)
  colnames(gene_scale2) <- colnames(DE_voom2)[i]
  DE_voom_scgen <- cbind(DE_voom_scgen, gene_scale2)
}
DE_voom_scgen <- t(DE_voom_scgen)

DE_voom_scgen2 <- DE_voom_scgen
colnames(DE_voom_scgen2) <- cell_line

breaksList = seq(range(DE_voom_scgen2)[1], range(DE_voom_scgen2)[2], by = 0.5)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
pheatmap(DE_voom_scgen2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames=TRUE, annotation_col=NA, col=colors, 
         breaks=breaksList, annotation_names_col=FALSE, main="Top 100 DE genes in iNeurons", fontsize=18)


# plot boxplots per gene #

perDEgene <- list()
for(i in 1:17584){
  norm_counts <- t(DE_voom[i,12:32])
  perline <- cbind.data.frame(cell_line = cell_line, norm_counts)
  colnames(perline)[2] <- "voom_cor"
  perDEgene[[i]] <- perline
}
names(perDEgene) <- DE_voom$Gene_name

boxplot(voom_cor~cell_line, perDEgene[[1]], main=DE_voom$Gene_name[[1]], xlab="", ylab="Normalized counts (voom)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)


#### perform GSEA using fgsea function on t-statistic ####

metadata97 <- read.table(".../0.Data/RNA_seq/counts_MCS/190123_EnsemblID_genes_GRCh38p12_ensembl97_metadata_biomart.txt", header=TRUE, sep="\t")
colnames(metadata97) <- c("Gene_stable_ID", "Gene_stable_ID_version", "Gene_name", "Chromosome", "Gene_type")

DE_res_tmp <- cbind.data.frame(Gene_stable_ID=DE_res$Gene_stable_ID, Gene_name=DE_res$Gene_name)
metadata97_tmp <- cbind.data.frame(Gene_stable_ID=metadata97$Gene_stable_ID, Gene_name=metadata97$Gene_name)

same <- intersect(DE_res_tmp, metadata97_tmp) # 17278 out of 17584
diff <- setdiff(DE_res_tmp, metadata97_tmp) # 306 out of 17584

same_tmp <- merge(same, metadata97_tmp, by="Gene_stable_ID")
colnames(same_tmp) <- c("Gene_stable_ID", "Gene_name", "Gene_name97")

diff_gs <- merge(diff, metadata97_tmp, by="Gene_stable_ID") # 273 have same ID, different name
colnames(diff_gs) <- c("Gene_stable_ID", "Gene_name", "Gene_name97")

same_final <- rbind.data.frame(same_tmp, diff_gs) # 17551 with unique Gene stable ID (same as rows), 17521 unique gene symbols
# 33 genes removed with different ID or no match at all 

DE_res97 <- merge(DE_res, same_final) # 17551 unique Gene stable ID, 17521 unique Gene name

ranks_symbol <- DE_res97$t
names(ranks_symbol) <- DE_res97$Gene_name97


# load gene set info from msigdbr #

human <- msigdbr(species = "Homo sapiens")
human_Reac <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
human_TFT <- msigdbr(species="Homo sapiens", category="C3", subcategory="TFT")
human_GO <- msigdbr(species="Homo sapiens", category="C5")


# Reactome #

Reactome <-  human_Reac %>% split(x = .$gene_symbol, f = .$gs_name)

Reactome_overlap <- c()
for (i in 1:1499){
  Reactome_overlap[[i]] <- Reactome[[i]][which(Reactome[[i]] %in% DE_res97$Gene_name97)]
}
names(Reactome_overlap) <- names(Reactome)

Reactome_tmp <- sapply(Reactome_overlap, function(i) length(i) > 5 & length(i) < 500)
Reactome_final <- Reactome_overlap[Reactome_tmp]


# TFT #

TFT <-  human_TFT %>% split(x = .$gene_symbol, f = .$gs_name)

TFT_overlap <- c()
for (i in 1:610){
  TFT_overlap[[i]] <- TFT[[i]][which(TFT[[i]] %in% DE_res97$Gene_name97)]
}
names(TFT_overlap) <- names(TFT)

TFT_tmp <- sapply(TFT_overlap, function(i) length(i) > 5 & length(i) < 500)
TFT_final <- TFT_overlap[TFT_tmp]


# GO #

GO <-  human_GO %>% split(x = .$gene_symbol, f = .$gs_name)

GO_overlap <- c()
for (i in 1:9996){
  GO_overlap[[i]] <- GO[[i]][which(GO[[i]] %in% DE_res97$Gene_name97)]
}
names(GO_overlap) <- names(GO)

GO_tmp <- sapply(GO_overlap, function(i) length(i) > 5 & length(i) < 500)
GO_final <- GO_overlap[GO_tmp]


# combine all gene sets on which the final analysis is performed #

genesets <- do.call(c, list(Reactome_final, TFT_final, GO_final))

GSEA_res <- fgsea(genesets, ranks_symbol, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig <- subset(GSEA_res, GSEA_res$padj < 0.05)
GSEA_sig <- GSEA_sig[order(GSEA_sig$padj), ] # 159 significant

GSEA_LEgenes <- c()
for (i in 1:159){
  GSEA_LEgenes[[i]] <- unnest(GSEA_sig[i,8])$leadingEdge
}
names(GSEA_LEgenes) <- GSEA_sig$pathway


# load MEA data #
# exclude 76 2311 D4 well #

MEA <- read.table(".../0.Data/MCS_MEA_data/200407_MEAdata_RNAseq_samples_outlExcl.txt", header=TRUE, row.names=1)
MEA2 <- t(MEA[,c(9:14,16:19)])

MEA_data <- MEA[,c(9:14,16:19)]
MEA_data2 <- cbind.data.frame(Status = c(paste(rep("Control", 9)), paste(rep("KdVS", 11))), Line=cell_line[c(1:16,18:21)], Line2 = cell_line[c(1:9,7:9,13:16,18:21)], MEA_data)

par(mfrow=c(2,2))

ggplot(MEA_data2, aes(y=mean_percentage_of_random_spikes, x=Status)) + geom_boxplot() + geom_point(aes(colour=Status, shape=Line2), position=position_jitterdodge(), size=5) + theme_minimal() + scale_color_manual(values=c("darkgrey", "darkred")) + ylab("Inter-burst interval (s)")
ggplot(MEA_data2, aes(y=mean_inter_burst_interval_.ms./1000, x=Status)) + geom_boxplot() + geom_point(aes(colour=Status, shape=Line2), position=position_jitterdodge(), size=5) + theme_minimal() + scale_color_manual(values=c("darkgrey", "darkred")) + ylab("Inter-burst interval (s)")
ggplot(MEA_data2, aes(y=mean_NB_rate_.nb.min., x=Status)) + geom_boxplot() + geom_point(aes(colour=Status, shape=Line2), position=position_jitterdodge(), size=5) + theme_minimal() + scale_color_manual(values=c("darkgrey", "darkred")) + ylab("Network burst rate (NB/min)")
ggplot(MEA_data2, aes(y=mean_NB_inter_burst_interval_.ms./1000, x=Status)) + geom_boxplot() + geom_point(aes(colour=Status, shape=Line2), position=position_jitterdodge(), size=5) + theme_minimal() + scale_color_manual(values=c("darkgrey", "darkred")) + ylab("Inter-network burst interval (s)")
ggplot(MEA_data2, aes(y=Coefficient_of_variation, x=Status)) + geom_boxplot() + geom_point(aes(colour=Status, shape=Line2), position=position_jitterdodge(), size=5) + theme_minimal() + scale_color_manual(values=c("darkgrey", "darkred")) + ylab("Coefficient of variation (CVNIBI)")



# use not-batch corrected counts #

rownames(voom$E) <- voom$genes$Gene_stable_ID_version
MEAseq <- rbind.data.frame(voom$E[,c(1:15,17:21)], MEA2)

MEAseq_cor <- corr.test(t(MEA2), t(voom$E[,c(1:15,17:21)]), method="spearman", adjust="BH", ci=FALSE)
MEAseq_cor_r <- t(MEAseq_cor$r)
MEAseq_cor_r <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_r), MEAseq_cor_r)
MEAseq_cor_r <- merge(metadata, MEAseq_cor_r, by="Gene_stable_ID_version")
MEAseq_cor_p <- t(MEAseq_cor$p)
MEAseq_cor_p <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_p), MEAseq_cor_p)
MEAseq_cor_p <- merge(metadata, MEAseq_cor_p, by="Gene_stable_ID_version")

NBR <- cbind.data.frame(MEAseq_cor_r[,c(1:5,12)], MEAseq_cor_p[,c(12)])
NBR2 <- NBR[order(NBR$mean_NB_rate_.nb.min.),]


# select significant per parameter #

param_cor <- vector("list", length=10)
for(i in 1:10){
    temp <- MEAseq_cor_p[MEAseq_cor_p[,5+i] < 0.05,]
    temp2 <- MEAseq_cor_r[MEAseq_cor_r$Gene_stable_ID_version %in% temp$Gene_stable_ID_version,]
    temp3 <- temp2[order(temp2[,5+i], decreasing=TRUE),]
    param_cor[[i]] <- cbind.data.frame(temp3[,1:5], cor = temp3[,5+i], padj = temp[,5+i])
}
names(param_cor) <- c("MFR", "PRS", "BR", "BD", "BSR", "IBI", "NBR", "NBD", "NIBI", "CV")

param_cor_sig_final <- do.call("rbind", param_cor) # 666
param_cor_sig_final <- cbind.data.frame(param = rownames(param_cor_sig_final), param_cor_sig_final)
param_cor_sig_final$param <- gsub("\\..*", "", param_cor_sig_final$param)

allgenes <- rbind(param_cor$PRS[,1:3], param_cor$BR[,1:3], param_cor$BD[,1:3], param_cor$IBI[,1:3], param_cor$NBR[,1:3], param_cor$NBD[,1:3], param_cor$NIBI[,1:3], param_cor$CV[,1:3]) # 461 unique
allgenes <- allgenes[!duplicated(allgenes$Gene_stable_ID_version),] # 435 unique genes

allgenes_cor_r <- MEAseq_cor_r[which(MEAseq_cor_r$Gene_stable_ID_version %in% allgenes$Gene_stable_ID_version),]
allgenes_cor_p <- MEAseq_cor_p[which(MEAseq_cor_p$Gene_stable_ID_version %in% allgenes$Gene_stable_ID_version),]



# plot correlations # 

rownames(MEA2) <- names(param_cor)

voom2 <- voom$E[,c(1:15,17:21)]
voom3 <- cbind.data.frame(voom$genes, voom2)
rownames(voom3) <- voom3$Gene_stable_ID_version

par(mfrow=c(3,3))

for(j in 4){
    for(i in 81:89){
      gene_tmp <- voom3[which(voom3$Gene_stable_ID_version %in% param_cor[[j]]$Gene_stable_ID_version[i]),]
      gene <- as.numeric(voom3[which(voom3$Gene_stable_ID_version %in% param_cor[[j]]$Gene_stable_ID_version[i]),c(6:25)])
      PRS <- as.numeric(MEA2[j,])
      # plot(PRS, gene, main=gene_tmp$Gene_name, xlab=rownames(MEA2)[j])
      scatter.smooth(PRS, gene, main=gene_tmp$Gene_name, xlab=rownames(MEA2)[j])
  }
}


# correlation plots #

corrplot_all <- ggcorrplot(t(allgenes_cor_r[,6:15]), method="circle", p.mat = t(allgenes_cor_p[,6:15]), insig="blank", ggtheme=theme_bw)


#########################################
############ posthoc analysis ###########
#########################################


### confirm findings in controls only ###

# per parameter #

param_cor2 <- param_cor[lapply(param_cor,nrow)>0]

param_ctl_cor <- vector("list", length=10)
for(i in 1:10){
    voom3_param <- voom3[which(voom3$Gene_stable_ID_version %in% param_cor2[[i]]$Gene_stable_ID_version),]
    voom3_param_ctl <- voom3_param[,1:14]
    MEA2_ctl <- as.data.frame(t(MEA2[,1:9]))
    MEAseq_cor_param_ctl <- corr.test(MEA2_ctl[,i], t(voom3_param_ctl[,c(6:14)]), method = "spearman", adjust="none", ci=FALSE)
    MEAseq_cor_param_ctl_r <- t(MEAseq_cor_param_ctl$r)
    MEAseq_cor_param_ctl_r <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_ctl_r), MEAseq_cor_param_ctl_r)
    MEAseq_cor_param_ctl_r <- merge(metadata, MEAseq_cor_param_ctl_r, by="Gene_stable_ID_version")
    MEAseq_cor_param_ctl_p <- t(MEAseq_cor_param_ctl$p)
    MEAseq_cor_param_ctl_p <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_ctl_p), MEAseq_cor_param_ctl_p)
    MEAseq_cor_param_ctl_p <- merge(metadata, MEAseq_cor_param_ctl_p, by="Gene_stable_ID_version")
    param_ctl_cor[[i]] <- cbind.data.frame(MEAseq_cor_param_ctl_r[,1:5], cor = MEAseq_cor_param_ctl_r$MEAseq_cor_param_ctl_r, padj = MEAseq_cor_param_ctl_p$MEAseq_cor_param_ctl_p)
}
names(param_ctl_cor) <- names(param_cor2)

param_ctl_cor_final <- do.call("rbind", param_ctl_cor) # 327
param_ctl_cor_final <- cbind.data.frame(param = rownames(param_ctl_cor_final), param_ctl_cor_final)
param_ctl_cor_final$param <- gsub("\\..*", "", param_ctl_cor_final$param)
param_ctl_cor_final2 <- merge(param_cor_sig_final, param_ctl_cor_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_ctl_cor_final2)[7:10] <- c("cor_all", "padj_all", "cor_ctl", "padj_ctl")

param_ctl_cor_sig <- vector("list", length=10)
for(i in 1:10){
    param_ctl_cor_sig[[i]] <- subset(param_ctl_cor[[i]], param_ctl_cor[[i]]$padj < 0.05)
}
names(param_ctl_cor_sig) <- names(param_ctl_cor)


param_ctl_cor_sig_final <- do.call("rbind", param_ctl_cor_sig) # 417
param_ctl_cor_sig_pos <- subset(param_ctl_cor_sig_final, param_ctl_cor_sig_final$cor > 0)
param_ctl_cor_sig_neg <- subset(param_ctl_cor_sig_final, param_ctl_cor_sig_final$cor < 0)

param_ctl_cor_sig_final <- cbind.data.frame(param = rownames(param_ctl_cor_sig_final), param_ctl_cor_sig_final)
param_ctl_cor_sig_final$param <- gsub("\\..*", "", param_ctl_cor_sig_final$param)

param_ctl_cor_sig_final2 <- merge(param_cor_sig_final, param_ctl_cor_sig_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_ctl_cor_sig_final2)[7:10] <- c("cor_all", "padj_all", "cor_ctl", "pval_ctl")




### IN ADDITION: confirm findings in patients ###

param_pat_cor <- vector("list", length=10)
for(i in 1:10){
  voom3_param <- voom3[which(voom3$Gene_stable_ID_version %in% param_cor2[[i]]$Gene_stable_ID_version),]
  voom3_param_pat <- voom3_param[,c(1:5,15:25)]
  MEA2_pat <- as.data.frame(t(MEA2[,10:20]))
  MEAseq_cor_param_pat <- corr.test(MEA2_pat[,i], t(voom3_param_pat[,c(6:16)]), method = "spearman", adjust="none", ci=FALSE)
  MEAseq_cor_param_pat_r <- t(MEAseq_cor_param_pat$r)
  MEAseq_cor_param_pat_r <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_pat_r), MEAseq_cor_param_pat_r)
  MEAseq_cor_param_pat_r <- merge(metadata, MEAseq_cor_param_pat_r, by="Gene_stable_ID_version")
  MEAseq_cor_param_pat_p <- t(MEAseq_cor_param_pat$p)
  MEAseq_cor_param_pat_p <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_pat_p), MEAseq_cor_param_pat_p)
  MEAseq_cor_param_pat_p <- merge(metadata, MEAseq_cor_param_pat_p, by="Gene_stable_ID_version")
  param_pat_cor[[i]] <- cbind.data.frame(MEAseq_cor_param_pat_r[,1:5], cor = MEAseq_cor_param_pat_r$MEAseq_cor_param_pat_r, pval = MEAseq_cor_param_pat_p$MEAseq_cor_param_pat_p)
}
names(param_pat_cor) <- names(param_cor2)

param_pat_cor_final <- do.call("rbind", param_pat_cor) # 327
param_pat_cor_final <- cbind.data.frame(param = rownames(param_pat_cor_final), param_pat_cor_final)
param_pat_cor_final$param <- gsub("\\..*", "", param_pat_cor_final$param)
param_pat_cor_final2 <- merge(param_ctl_cor_final2, param_pat_cor_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_pat_cor_final2)[11:12] <- c("cor_pat", "padj_pat")


param_pat_cor_sig <- vector("list", length=10)
for(i in 1:10){
  param_pat_cor_sig[[i]] <- subset(param_pat_cor[[i]], param_pat_cor[[i]]$pval < 0.05)
}
names(param_pat_cor_sig) <- names(param_pat_cor)


param_pat_cor_sig_final <- do.call("rbind", param_pat_cor_sig) # 327
param_pat_cor_sig_pos <- subset(param_pat_cor_sig_final, param_pat_cor_sig_final$cor > 0)
param_pat_cor_sig_neg <- subset(param_pat_cor_sig_final, param_pat_cor_sig_final$cor < 0)

param_pat_cor_sig_final <- cbind.data.frame(param = rownames(param_pat_cor_sig_final), param_pat_cor_sig_final)
param_pat_cor_sig_final$param <- gsub("\\..*", "", param_pat_cor_sig_final$param)

param_pat_cor_sig_final2 <- merge(param_ctl_cor_sig_final2, param_pat_cor_sig_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_pat_cor_sig_final2)[11:12] <- c("cor_pat", "pval_pat")




### confirm findings in female lines only ###

# per parameter #

param_fem_cor <- vector("list", length=10)
for(i in 1:10){
  voom3_param <- voom3[which(voom3$Gene_stable_ID_version %in% param_cor2[[i]]$Gene_stable_ID_version),]
  voom3_param_fem <- voom3_param[,c(1:5,9:25)]
  MEA2_fem <- as.data.frame(t(MEA2[,4:20]))
  MEAseq_cor_param_fem <- corr.test(MEA2_fem[,i], t(voom3_param_fem[,c(6:22)]), method = "spearman", adjust="none", ci=FALSE)
  MEAseq_cor_param_fem_r <- t(MEAseq_cor_param_fem$r)
  MEAseq_cor_param_fem_r <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_fem_r), MEAseq_cor_param_fem_r)
  MEAseq_cor_param_fem_r <- merge(metadata, MEAseq_cor_param_fem_r, by="Gene_stable_ID_version")
  MEAseq_cor_param_fem_p <- t(MEAseq_cor_param_fem$p)
  MEAseq_cor_param_fem_p <- cbind.data.frame(Gene_stable_ID_version = rownames(MEAseq_cor_param_fem_p), MEAseq_cor_param_fem_p)
  MEAseq_cor_param_fem_p <- merge(metadata, MEAseq_cor_param_fem_p, by="Gene_stable_ID_version")
  param_fem_cor[[i]] <- cbind.data.frame(MEAseq_cor_param_fem_r[,1:5], cor = MEAseq_cor_param_fem_r$MEAseq_cor_param_fem_r, pval = MEAseq_cor_param_fem_p$MEAseq_cor_param_fem_p)
}
names(param_fem_cor) <- names(param_cor2)

param_fem_cor_final <- do.call("rbind", param_fem_cor) # 327
param_fem_cor_final <- cbind.data.frame(param = rownames(param_fem_cor_final), param_fem_cor_final)
param_fem_cor_final$param <- gsub("\\..*", "", param_fem_cor_final$param)
param_fem_cor_final2 <- merge(param_pat_cor_final2, param_fem_cor_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_fem_cor_final2)[13:14] <- c("cor_fem", "padj_fem")

param_fem_cor_sig <- vector("list", length=10)
for(i in 1:10){
  param_fem_cor_sig[[i]] <- subset(param_fem_cor[[i]], param_fem_cor[[i]]$pval < 0.05)
}
names(param_fem_cor_sig) <- names(param_fem_cor)


param_fem_cor_sig_final <- do.call("rbind", param_fem_cor_sig) # 327
param_fem_cor_sig_pos <- subset(param_fem_cor_sig_final, param_fem_cor_sig_final$cor > 0)
param_fem_cor_sig_neg <- subset(param_fem_cor_sig_final, param_fem_cor_sig_final$cor < 0)

param_fem_cor_sig_final <- cbind.data.frame(param = rownames(param_fem_cor_sig_final), param_fem_cor_sig_final)
param_fem_cor_sig_final$param <- gsub("\\..*", "", param_fem_cor_sig_final$param)

param_fem_cor_sig_final2 <- merge(param_pat_cor_sig_final2, param_fem_cor_sig_final, by=c("param", "Gene_stable_ID_version", "Gene_stable_ID", "Gene_name", "Chromosome", "Gene_type"))
colnames(param_fem_cor_sig_final2)[13:14] <- c("cor_fem", "pval_fem")


## final posthoc results (control only, patient only, female only) ##

param_cor_posthoc_control <- subset(param_fem_cor_final2, param_fem_cor_final2$padj_all < 0.05 & param_fem_cor_final2$padj_ctl < 0.05) # 327
param_cor_posthoc_control2 <- subset(param_cor_posthoc_control, param_cor_posthoc_control$cor_all > 0 & param_cor_posthoc_control$cor_ctl > 0 | param_cor_posthoc_control$cor_all < 0 & param_cor_posthoc_control$cor_ctl < 0 )
param_cor_posthoc_patient <- subset(param_fem_cor_final2, param_fem_cor_final2$padj_all < 0.05 & param_fem_cor_final2$padj_pat < 0.05) # 562
param_cor_posthoc_patient2 <- subset(param_cor_posthoc_patient, param_cor_posthoc_patient$cor_all > 0 & param_cor_posthoc_patient$cor_pat > 0 | param_cor_posthoc_patient$cor_all < 0 & param_cor_posthoc_patient$cor_pat < 0 )
param_cor_posthoc_female <- subset(param_fem_cor_final2, param_fem_cor_final2$padj_all < 0.05 & param_fem_cor_final2$padj_fem < 0.05) # 666
param_cor_posthoc_female2 <- subset(param_cor_posthoc_female, param_cor_posthoc_female$cor_all > 0 & param_cor_posthoc_female$cor_fem > 0 | param_cor_posthoc_female$cor_all < 0 & param_cor_posthoc_female$cor_fem < 0 )
# all correlations significant are in the same direction = checked #

param_cor_posthoc_final <- subset(param_fem_cor_final2, param_fem_cor_final2$padj_all < 0.05 & param_fem_cor_final2$padj_ctl < 0.05 & param_fem_cor_final2$padj_pat < 0.05 & param_fem_cor_final2$padj_fem < 0.05)

param_cor_posthoc_final_param <- subset(param_cor_posthoc_final, param_cor_posthoc_final$param == "PRS" | param_cor_posthoc_final$param == "IBI" | param_cor_posthoc_final$param =="NBR"  | param_cor_posthoc_final$param == "NIBI" | param_cor_posthoc_final$param == "CV")
DEres_cor_posthoc_final_param <- DE_res[which(DE_res$Gene_stable_ID_version %in% param_cor_posthoc_final_param$Gene_stable_ID_version),]

# write.table(param_cor_posthoc_final, "C:/Users/Anouk/Google Drive/Documenten/Experiments/Sequencing_results/KdVs_RNAseq_iPSCs_iNeurons_final/Results/iNeurons/Tables/MEAseq/200629_MEAseq_correlations_final_with_posthoc.txt")

# correlation plots #

IBINBRNIBICV <- param_cor_posthoc_final[which(param_cor_posthoc_final$param == "PRS" | param_cor_posthoc_final$param == "IBI" | param_cor_posthoc_final$param == "NBR" | param_cor_posthoc_final$param == "NIBI" |  param_cor_posthoc_final$param == "CV"),]
MEAseq_cor_r_NNC <- MEAseq_cor_r[which(MEAseq_cor_r$Gene_stable_ID_version %in% IBINBRNIBICV$Gene_stable_ID_version),]
MEAseq_cor_r_NNC <- MEAseq_cor_r_NNC[!grepl("ENSG00000284779.1", MEAseq_cor_r_NNC$Gene_stable_ID_version),]
MEAseq_cor_r_NNC2 <- MEAseq_cor_r_NNC[,c(1:5,7,11:12,14:15)]
colnames(MEAseq_cor_r_NNC2)[6:10] <- c("PRS", "IBI", "NBR", "NIBI", "CV")
rownames(MEAseq_cor_r_NNC2) <- MEAseq_cor_r_NNC2$Gene_name

MEAseq_cor_p_NNC <- MEAseq_cor_p[which(MEAseq_cor_p$Gene_stable_ID_version %in% IBINBRNIBICV$Gene_stable_ID_version),]
MEAseq_cor_p_NNC <- MEAseq_cor_p_NNC[!grepl("ENSG00000284779.1", MEAseq_cor_p_NNC$Gene_stable_ID_version),]
MEAseq_cor_p_NNC2 <- MEAseq_cor_p_NNC[,c(1:5,7,11:12,14:15)]
colnames(MEAseq_cor_p_NNC2)[6:10] <- c("PRS", "IBI", "NBR", "NIBI", "CV")
rownames(MEAseq_cor_p_NNC2) <- MEAseq_cor_p_NNC2$Gene_name

NBR_not <- c("STS", "TLE4", "CTSF")
NIBI_not <- c("RBBP7", "UROD", "PRKAB2", "TBX2-AS1", "AP000350.6", "AL031658.2")

MEAseq_cor_p_NNC2.1 <- within(MEAseq_cor_p_NNC2, {
  f <- Gene_name == 'STS' | Gene_name == "TLE4" | Gene_name== "CTSF"
  NBR[f] <- 1
}) 

MEAseq_cor_p_NNC2.2 <- within(MEAseq_cor_p_NNC2.1, {
  f <- Gene_name == "RBBP7" | Gene_name == "UROD" | Gene_name== "TBX2-AS1" | Gene_name ==  "AP000350.6" | Gene_name == "AL031658.2"
  NIBI[f] <- 1
}) 

corrplot_all <- ggcorrplot(as.matrix(MEAseq_cor_r_NNC2[,6:10]), insig="blank", p.mat = as.matrix(MEAseq_cor_p_NNC2.2[,6:10]), method="circle", ggtheme=theme_bw)

DEres_IBINBRNIBICV <- DE_res[which(DE_res$Gene_stable_ID_version %in% IBINBRNIBICV$Gene_stable_ID_version),]
DEres_IBINBRNIBICV_trend <- subset(DEres_IBINBRNIBICV, DEres_IBINBRNIBICV$P.Value < 0.05) 
MEAseq_cor_r_NNC3 <- MEAseq_cor_r_NNC2[which(MEAseq_cor_r_NNC2$Gene_stable_ID_version %in% DEres_IBINBRNIBICV_trend$Gene_stable_ID_version),]
MEAseq_cor_p_NNC3 <- MEAseq_cor_p_NNC2.2[which(MEAseq_cor_p_NNC2.2$Gene_stable_ID_version %in% DEres_IBINBRNIBICV_trend$Gene_stable_ID_version),]

corrplot_DE <- ggcorrplot(as.matrix(MEAseq_cor_r_NNC3[,6:10]), insig="blank", p.mat = as.matrix(MEAseq_cor_p_NNC3[,6:10]), method="circle", ggtheme=theme_bw)

# heatmap DE trend genes linked to MEA parameters 4 #

DEres_IBINBRNIBICV_trendvoom <- DE_voom[which(DE_voom$Gene_stable_ID_version %in% DEres_IBINBRNIBICV_trend$Gene_stable_ID_version),]
rownames(DEres_IBINBRNIBICV_trendvoom) <- DEres_IBINBRNIBICV_trendvoom$Gene_name
DEres_IBINBRNIBICV_trendvoom2 <- DEres_IBINBRNIBICV_trendvoom[order(DEres_IBINBRNIBICV_trendvoom$adj.P.Val),]
DEres_IBINBRNIBICV_trendvoom2 <- as.data.frame(t(DEres_IBINBRNIBICV_trendvoom2[,12:32]))

gene_scale <- scale(DEres_IBINBRNIBICV_trendvoom2[,1])
colnames(gene_scale) <- colnames(DEres_IBINBRNIBICV_trendvoom2)[1]
rownames(gene_scale) <- rownames(DEres_IBINBRNIBICV_trendvoom2)
DEres_IBINBRNIBICV_trendvoom_scgen <- data.frame(gene_scale)
for (i in 2:11){
  gene_scale <- scale(DEres_IBINBRNIBICV_trendvoom2[,i])
  gene_scale2 <- data.frame(gene_scale)
  colnames(gene_scale2) <- colnames(DEres_IBINBRNIBICV_trendvoom2)[i]
  DEres_IBINBRNIBICV_trendvoom_scgen <- cbind(DEres_IBINBRNIBICV_trendvoom_scgen, gene_scale2)
}
DEres_IBINBRNIBICV_trendvoom_scgen <- t(DEres_IBINBRNIBICV_trendvoom_scgen)

DEres_IBINBRNIBICV_trendvoom_scgen2 <- DEres_IBINBRNIBICV_trendvoom_scgen
colnames(DEres_IBINBRNIBICV_trendvoom_scgen2) <- cell_line

breaksList = seq(range(DEres_IBINBRNIBICV_trendvoom_scgen2)[1], range(DEres_IBINBRNIBICV_trendvoom_scgen2)[2], by = 0.5)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
pheatmap(DEres_IBINBRNIBICV_trendvoom_scgen2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, show_colnames=TRUE, annotation_col=NA, col=colors, 
         breaks=breaksList, annotation_names_col=FALSE, main="Trend DE genes correlated to 5 MEA param", fontsize=14)




# perform gene set enrichment analysis on correlated genes - NBR #

NBR_cor <- MEAseq_cor_r[,c(1:5,12)]
colnames(NBR_cor)[6] <- "NBR"

NBR_cor_p <- MEAseq_cor_p[,c(1:5,12)]
colnames(NBR_cor_p)[6] <- "NBR_padj"

NBR_cor2 <- merge(NBR_cor, DE_res97[,c(1:5,12)])
NBR_cor_p2 <- merge(NBR_cor2, NBR_cor_p)

### formula not correct anymore, check what you did !!! ### 
ranks_symbol_NBR <- NBR_cor_p2$NBR #* -log10
names(ranks_symbol_NBR) <- NBR_cor_p2$Gene_name97
ranks_symbol_NBR2 <- sort(ranks_symbol_NBR)


# load gene set info from msigdbr #

human <- msigdbr(species = "Homo sapiens")
human_Reac <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
human_GO <- msigdbr(species="Homo sapiens", category="C5")


# Reactome #

Reactome <-  human_Reac %>% split(x = .$gene_symbol, f = .$gs_name)

Reactome_overlap <- c()
for (i in 1:length(Reactome)){
  Reactome_overlap[[i]] <- Reactome[[i]][which(Reactome[[i]] %in% DE_res97$Gene_name97)]
}
names(Reactome_overlap) <- names(Reactome)

Reactome_tmp <- sapply(Reactome_overlap, function(i) length(i) > 5 & length(i) < 500)
Reactome_final <- Reactome_overlap[Reactome_tmp]

# GO #

GO <-  human_GO %>% split(x = .$gene_symbol, f = .$gs_name)

GO_overlap <- c()
for (i in 1:length(GO)){
  GO_overlap[[i]] <- GO[[i]][which(GO[[i]] %in% DE_res97$Gene_name97)]
}
names(GO_overlap) <- names(GO)

GO_tmp <- sapply(GO_overlap, function(i) length(i) > 5 & length(i) < 500)
GO_final <- GO_overlap[GO_tmp]


# combine all gene sets on which the final analysis is performed #

genesets <- do.call(c, list(Reactome_final, GO_final))

GSEA_res_NBR <- fgsea(genesets, ranks_symbol_NBR2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_NBR <- subset(GSEA_res_NBR, GSEA_res_NBR$padj < 0.05)
GSEA_sig_NBR <- GSEA_sig_NBR[order(GSEA_sig_NBR$NES), ] # 159 significant

GSEA_LEgenes_NBR <- c()
for (i in 1:nrow(GSEA_LEgenes_NBR)){
 GSEA_LEgenes_NBR[[i]] <- unnest(GSEA_sig_NBR[i,8])$leadingEdge
}
names(GSEA_LEgenes_NBR) <- GSEA_sig_NBR$pathway




# perform gene set enrichment analysis on correlated genes - IBI #

IBI_cor <- MEAseq_cor_r[,c(1:5,11)]
colnames(IBI_cor)[6] <- "IBI"

IBI_cor_p <- MEAseq_cor_p[,c(1:5,11)]
colnames(IBI_cor_p)[6] <- "IBI_padj"

IBI_cor2 <- merge(IBI_cor, DE_res97[,c(1:5,12)])
IBI_cor_p2 <- merge(IBI_cor2, IBI_cor_p)

ranks_symbol_IBI <- IBI_cor_p2$IBI #* -log10
names(ranks_symbol_IBI) <- IBI_cor_p2$Gene_name97
ranks_symbol_IBI2 <- sort(ranks_symbol_IBI)

GSEA_res_IBI <- fgsea(genesets, ranks_symbol_IBI2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_IBI <- subset(GSEA_res_IBI, GSEA_res_IBI$padj < 0.05)
GSEA_sig_IBI <- GSEA_sig_IBI[order(GSEA_sig_IBI$NES), ] # 159 significant

GSEA_LEgenes_IBI <- c()
for (i in 1:nrow(GSEA_LEgenes_IBI)){
  GSEA_LEgenes_IBI[[i]] <- unnest(GSEA_sig_IBI[i,8])$leadingEdge
}
names(GSEA_LEgenes_IBI) <- GSEA_sig_IBI$pathway



# perform gene set enrichment analysis on correlated genes - PRS #

PRS_cor <- MEAseq_cor_r[,c(1:5,7)]
colnames(PRS_cor)[6] <- "PRS"

PRS_cor_p <- MEAseq_cor_p[,c(1:5,7)]
colnames(PRS_cor_p)[6] <- "PRS_padj"

PRS_cor2 <- merge(PRS_cor, DE_res97[,c(1:5,12)])
PRS_cor_p2 <- merge(PRS_cor2, PRS_cor_p)

ranks_symbol_PRS <- PRS_cor_p2$PRS #* -log10
names(ranks_symbol_PRS) <- PRS_cor_p2$Gene_name97
ranks_symbol_PRS2 <- sort(ranks_symbol_PRS)

GSEA_res_PRS <- fgsea(genesets, ranks_symbol_PRS2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_PRS <- subset(GSEA_res_PRS, GSEA_res_PRS$padj < 0.05)
GSEA_sig_PRS <- GSEA_sig_PRS[order(GSEA_sig_PRS$NES), ] # 159 significant

GSEA_LEgenes_PRS <- c()
for (i in 1:nrow(GSEA_LEgenes_PRS)){
  GSEA_LEgenes_PRS[[i]] <- unnest(GSEA_sig_PRS[i,8])$leadingEdge
}
names(GSEA_LEgenes_PRS) <- GSEA_sig_PRS$pathway



# perform gene set enrichment analysis on correlated genes - NIBI #

NIBI_cor <- MEAseq_cor_r[,c(1:5,14)]
colnames(NIBI_cor)[6] <- "NIBI"

NIBI_cor_p <- MEAseq_cor_p[,c(1:5,14)]
colnames(NIBI_cor_p)[6] <- "NIBI_padj"

NIBI_cor2 <- merge(NIBI_cor, DE_res97[,c(1:5,12)])
NIBI_cor_p2 <- merge(NIBI_cor2, NIBI_cor_p)

ranks_symbol_NIBI <- NIBI_cor_p2$NIBI #* -log10
names(ranks_symbol_NIBI) <- NIBI_cor_p2$Gene_name97
ranks_symbol_NIBI2 <- sort(ranks_symbol_NIBI)

GSEA_res_NIBI <- fgsea(genesets, ranks_symbol_NIBI2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_NIBI <- subset(GSEA_res_NIBI, GSEA_res_NIBI$padj < 0.05)
GSEA_sig_NIBI <- GSEA_sig_NIBI[order(GSEA_sig_NIBI$NES), ] # 159 significant

GSEA_LEgenes_NIBI <- c()
for (i in 1:nrow(GSEA_LEgenes_NIBI)){
  GSEA_LEgenes_NIBI[[i]] <- unnest(GSEA_sig_NIBI[i,8])$leadingEdge
}
names(GSEA_LEgenes_NIBI) <- GSEA_sig_NIBI$pathway



# perform gene set enrichment analysis on correlated genes - CV #

CV_cor <- MEAseq_cor_r[,c(1:5,15)]
colnames(CV_cor)[6] <- "CV"

CV_cor_p <- MEAseq_cor_p[,c(1:5,15)]
colnames(CV_cor_p)[6] <- "CV_padj"

CV_cor2 <- merge(CV_cor, DE_res97[,c(1:5,12)])
CV_cor_p2 <- merge(CV_cor2, CV_cor_p)

ranks_symbol_CV <- CV_cor_p2$CV #* -log10
names(ranks_symbol_CV) <- CV_cor_p2$Gene_name97
ranks_symbol_CV2 <- sort(ranks_symbol_CV)

GSEA_res_CV <- fgsea(genesets, ranks_symbol_CV2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_CV <- subset(GSEA_res_CV, GSEA_res_CV$padj < 0.05)
GSEA_sig_CV <- GSEA_sig_CV[order(GSEA_sig_CV$NES), ] # 159 significant

GSEA_LEgenes_CV <- c()
for (i in 1:nrow(GSEA_LEgenes_CV)){
  GSEA_LEgenes_CV[[i]] <- unnest(GSEA_sig_CV[i,8])$leadingEdge
}
names(GSEA_LEgenes_CV) <- GSEA_sig_CV$pathway



# perform gene set enrichment analysis on correlated genes - MFR #

MFR_cor <- MEAseq_cor_r[,c(1:6)]
colnames(MFR_cor)[6] <- "MFR"

MFR_cor_p <- MEAseq_cor_p[,c(1:6)]
colnames(MFR_cor_p)[6] <- "MFR_padj"

MFR_cor2 <- merge(MFR_cor, DE_res97[,c(1:5,12)])
MFR_cor_p2 <- merge(MFR_cor2, MFR_cor_p)

ranks_symbol_MFR <- MFR_cor_p2$MFR #* -log10
names(ranks_symbol_MFR) <- MFR_cor_p2$Gene_name97
ranks_symbol_MFR2 <- sort(ranks_symbol_MFR)

GSEA_res_MFR <- fgsea(genesets, ranks_symbol_MFR2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_MFR <- subset(GSEA_res_MFR, GSEA_res_MFR$padj < 0.05)
GSEA_sig_MFR <- GSEA_sig_MFR[order(GSEA_sig_MFR$NES), ] # 159 significant

GSEA_LEgenes_MFR <- c()
for (i in 1:nrow(GSEA_LEgenes_MFR)){
  GSEA_LEgenes_MFR[[i]] <- unnest(GSEA_sig_MFR[i,8])$leadingEdge
}
names(GSEA_LEgenes_MFR) <- GSEA_sig_MFR$pathway



# perform gene set enrichment analysis on correlated genes - BR #

BR_cor <- MEAseq_cor_r[,c(1:5,8)]
colnames(BR_cor)[6] <- "BR"

BR_cor_p <- MEAseq_cor_p[,c(1:5,8)]
colnames(BR_cor_p)[6] <- "BR_padj"

BR_cor2 <- merge(BR_cor, DE_res97[,c(1:5,12)])
BR_cor_p2 <- merge(BR_cor2, BR_cor_p)

ranks_symbol_BR <- BR_cor_p2$BR #* -log10
names(ranks_symbol_BR) <- BR_cor_p2$Gene_name97
ranks_symbol_BR2 <- sort(ranks_symbol_BR)

GSEA_res_BR <- fgsea(genesets, ranks_symbol_BR2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_BR <- subset(GSEA_res_BR, GSEA_res_BR$padj < 0.05)
GSEA_sig_BR <- GSEA_sig_BR[order(GSEA_sig_BR$NES), ] # 159 significant

GSEA_LEgenes_BR <- c()
for (i in 1:nrow(GSEA_LEgenes_BR)){
  GSEA_LEgenes_BR[[i]] <- unnest(GSEA_sig_BR[i,8])$leadingEdge
}
names(GSEA_LEgenes_BR) <- GSEA_sig_BR$pathway


# perform gene set enrichment analysis on correlated genes - BD #

BD_cor <- MEAseq_cor_r[,c(1:5,9)]
colnames(BD_cor)[6] <- "BD"

BD_cor_p <- MEAseq_cor_p[,c(1:5,9)]
colnames(BD_cor_p)[6] <- "BD_padj"

BD_cor2 <- merge(BD_cor, DE_res97[,c(1:5,12)])
BD_cor_p2 <- merge(BD_cor2, BD_cor_p)

ranks_symbol_BD <- BD_cor_p2$BD #* -log10
names(ranks_symbol_BD) <- BD_cor_p2$Gene_name97
ranks_symbol_BD2 <- sort(ranks_symbol_BD)

GSEA_res_BD <- fgsea(genesets, ranks_symbol_BD2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_BD <- subset(GSEA_res_BD, GSEA_res_BD$padj < 0.05)
GSEA_sig_BD <- GSEA_sig_BD[order(GSEA_sig_BD$NES), ] # 159 significant

GSEA_LEgenes_BD <- c()
for (i in 1:nrow(GSEA_LEgenes_BD)){
  GSEA_LEgenes_BD[[i]] <- unnest(GSEA_sig_BD[i,8])$leadingEdge
}
names(GSEA_LEgenes_BD) <- GSEA_sig_BD$pathway




# perform gene set enrichment analysis on correlated genes - BSR #

BSR_cor <- MEAseq_cor_r[,c(1:5,10)]
colnames(BSR_cor)[6] <- "BSR"

BSR_cor_p <- MEAseq_cor_p[,c(1:5,10)]
colnames(BSR_cor_p)[6] <- "BSR_padj"

BSR_cor2 <- merge(BSR_cor, DE_res97[,c(1:5,12)])
BSR_cor_p2 <- merge(BSR_cor2, BSR_cor_p)

ranks_symbol_BSR <- BSR_cor_p2$BSR #* -log10
names(ranks_symbol_BSR) <- BSR_cor_p2$Gene_name97
ranks_symbol_BSR2 <- sort(ranks_symbol_BSR)

GSEA_res_BSR <- fgsea(genesets, ranks_symbol_BSR2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_BSR <- subset(GSEA_res_BSR, GSEA_res_BSR$padj < 0.05)
GSEA_sig_BSR <- GSEA_sig_BSR[order(GSEA_sig_BSR$NES), ] # 159 significant

GSEA_LEgenes_BSR <- c()
for (i in 1:nrow(GSEA_LEgenes_BSR)){
  GSEA_LEgenes_BSR[[i]] <- unnest(GSEA_sig_BSR[i,8])$leadingEdge
}
names(GSEA_LEgenes_BSR) <- GSEA_sig_BSR$pathway



# perform gene set enrichment analysis on correlated genes - NBD #

NBD_cor <- MEAseq_cor_r[,c(1:5,13)]
colnames(NBD_cor)[6] <- "NBD"

NBD_cor_p <- MEAseq_cor_p[,c(1:5,13)]
colnames(NBD_cor_p)[6] <- "NBD_padj"

NBD_cor2 <- merge(NBD_cor, DE_res97[,c(1:5,12)])
NBD_cor_p2 <- merge(NBD_cor2, NBD_cor_p)

ranks_symbol_NBD <- NBD_cor_p2$NBD #* -log10
names(ranks_symbol_NBD) <- NBD_cor_p2$Gene_name97
ranks_symbol_NBD2 <- sort(ranks_symbol_NBD)

GSEA_res_NBD <- fgsea(genesets, ranks_symbol_NBD2, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig_NBD <- subset(GSEA_res_NBD, GSEA_res_NBD$padj < 0.05)
GSEA_sig_NBD <- GSEA_sig_NBD[order(GSEA_sig_NBD$NES), ] # 159 significant

GSEA_LEgenes_NBD <- c()
for (i in 1:nrow(GSEA_LEgenes_NBD)){
  GSEA_LEgenes_NBD[[i]] <- unnest(GSEA_sig_NBD[i,8])$leadingEdge
}
names(GSEA_LEgenes_NBD) <- GSEA_sig_NBD$pathway




### PLOTS ###


## make bubble plot with ggplot2 ##

library(dplyr)
library(ggplot2)
library(forcats)

# subset reactome pathways and BP GO terms #

human_GOBP <- msigdbr(species="Homo sapiens", category="C5", subcategory="BP")

geneset_finalBP <- GSEA_sig_NBR[which(GSEA_sig_NBR$pathway %in% human_GOBP$gs_name),] # 54
geneset_finalRe <- GSEA_sig_NBR[grepl("^REACTOME_", GSEA_sig_NBR$pathway),] # 39
geneset_final_tmp <- rbind.data.frame(geneset_finalBP, geneset_finalRe) # 93

geneset_final <- geneset_final_tmp[order(geneset_final_tmp$NES, decreasing=TRUE),]
geneset_final_top10 <- geneset_final[1:20,]

generatio_top <- c()
for (i in 1:20){
  temp <- GSEA_LEgenes[grepl(paste(geneset_final_top10$pathway[i],"\\.", sep=""), rownames(GSEA_LEgenes_NBR)),]
  generatio_top[[i]] <- length(temp)
} 

geneset_final2_top10 <- cbind.data.frame(geneset_final_top10, GeneRatio = generatio_top / geneset_final_top10$size, logp = -log10(geneset_final_top10$padj))

ggplot(geneset_final2_top10[1:20,], aes(x=NES, y=pathway, size=size, fill=GeneRatio)) + 
  geom_point(alpha=1, shape=21, color="black") + aes(y = fct_reorder(pathway, NES, .desc=F)) + 
  scale_fill_gradient(low="blue", high="red", name="Gene ratio", limits=c(0.2,1.0), breaks=c(0.2,0.4,0.6,0.8,1.0)) + 
  scale_size(range=c(4,10),  name="Size", limits=c(0,450), breaks=c(0,150,300,450)) + 
  xlab("Normalized Enrichment Score (NES)") + ylab("") + theme(axis.title.x=element_text(size=12, face="bold"), axis.text=element_text(size=12), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12)) +
  guides(size=guide_legend(order=1)) + scale_x_continuous(breaks=c(1.8,1.9, 2.0, 2.1, 2.2), limits=c(1.8,2.2)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size=15, face="bold")) + ggtitle("")

geneset_final_bot <- geneset_final_tmp[order(geneset_final_tmp$NES),]
geneset_final_bot10 <- geneset_final_bot[1:20,]

generatio_bot <- c()
for (i in 1:20){
  temp <- GSEA_LEgenes[grepl(paste(geneset_final_bot10$pathway[i],"\\.", sep=""), rownames(GSEA_LEgenes_NBR)),]
  generatio_bot[[i]] <- length(temp)
} 

geneset_final2_bot10 <- cbind.data.frame(geneset_final_bot10, GeneRatio = generatio_bot / geneset_final_bot10$size, logp = -log10(geneset_final_bot10$padj))

ggplot(geneset_final2_bot10[1:20,], aes(x=NES, y=pathway, size=size, fill=GeneRatio)) + 
  geom_point(alpha=1, shape=21, color="black") + aes(y = fct_reorder(pathway, NES, .desc=F)) + 
  scale_fill_gradient(low="blue", high="red", name="Gene ratio", limits=c(0.2,1.0), breaks=c(0.2,0.4,0.6,0.8,1.0)) + 
  scale_size(range=c(4,10),  name="Size", limits=c(0,450), breaks=c(0,150,300,450)) + 
  xlab("Normalized Enrichment Score (NES)") + ylab("") + theme(axis.title.x=element_text(size=12, face="bold"), axis.text=element_text(size=12), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12)) +
  guides(size=guide_legend(order=1)) + scale_x_continuous(breaks=c(-2.2, -2.1, -2.0, -1.9), limits=c(-2.2,-1.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size=15, face="bold")) + ggtitle("")






