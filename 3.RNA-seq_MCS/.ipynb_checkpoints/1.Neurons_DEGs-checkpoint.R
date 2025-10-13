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
library("userfriendlyscience")
library("GOplot")


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

write.table(iNeuron_hgenes3, "../0.Data/RNA_seq/counts_MCS/rawcounts_allgenes.txt")


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



### subset counts table for genes with cpm > 2 in 3 repl ###

counts_cpm2 <- iNeuron_hgenes2[iNeuron_hgenes2$Gene_stable_ID_version %in% cpm_thr2$Gene_stable_ID_version,]
counts_final <- cbind.data.frame(counts_cpm2[,1:5], counts_cpm2[,grepl("WTC", colnames(counts_cpm2))], counts_cpm2[,grepl("409B2", colnames(counts_cpm2))], counts_cpm2[,grepl("77", colnames(counts_cpm2))], counts_cpm2[,grepl("CRISPR", colnames(counts_cpm2))], counts_cpm2[,grepl("71", colnames(counts_cpm2))], counts_cpm2[,grepl("76", colnames(counts_cpm2))], counts_cpm2[,grepl("70", colnames(counts_cpm2))])
setwd("C:/Users/Anouk/Google Drive/Documenten/Experiments/Sequencing_results/KdVs_RNAseq_iPSCs_iNeurons_final/Results/iNeurons/Tables/")
# write.table(counts_final, "200109_rawcounts_genes_cpm_above2_in_3repl.txt")


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

breaksList = seq(range(DE_genes_voom_scgen2)[1], range(DE_genes_voom_scgen2)[2], by = 0.5)
# breaksList = seq(-2.5, 2, by = 0.5)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
pheatmap(DE_genes_voom_scgen2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames=TRUE, annotation_col=NA, col=colors, 
         breaks=breaksList, annotation_names_col=FALSE, main="DE genes in iNeurons", fontsize=18)

DE_genes_voom2 <- as.data.frame(t(DE_genes_voom[,12:17]))

gene_scale <- scale(DE_genes_voom2[,1])
colnames(gene_scale) <- colnames(DE_genes_voom2)[1]
rownames(gene_scale) <- rownames(DE_genes_voom2)
DE_genes_voom_scgen <- data.frame(gene_scale)
for (i in 2:442){
  gene_scale <- scale(DE_genes_voom2[,i])
  gene_scale2 <- data.frame(gene_scale)
  colnames(gene_scale2) <- colnames(DE_genes_voom2)[i]
  DE_genes_voom_scgen <- cbind(DE_genes_voom_scgen, gene_scale2)
}
DE_genes_voom_scgen <- t(DE_genes_voom_scgen)

DE_genes_voom_scgen2 <- DE_genes_voom_scgen
colnames(DE_genes_voom_scgen2) <- cell_line

breaksList = seq(-2, 2, by = 0.5)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList))
pheatmap(DE_genes_voom_scgen2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, show_colnames=TRUE, annotation_col=NA, col=colors, 
         breaks=breaksList, annotation_names_col=FALSE, main="DE genes in iNeurons", fontsize=18)


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



# volcano plot DE results #

plot(DE_res$logFC, -log10(DE_res$P.Value), pch=20, ylab="-log10(p-value)", xlab="LogFC", cex=1, xlim=c(-7,7))
ind <- DE_res$adj.P.Val < 0.05 & DE_res$logFC >0
ind2 <- DE_res$adj.P.Val < 0.05 & DE_res$logFC < 0
points(DE_res$logFC[ind], -log10(DE_res$P.Value)[ind], pch = 20, col = "#A50026", cex=1.5)
points(DE_res$logFC[ind2], -log10(DE_res$P.Value)[ind2], pch = 20, col = "#313695", cex=1.5)
sign.genes=which(DE_res$adj.P.Val<0.05)[1:20]
text(x=DE_res$logFC[sign.genes]+0.2, y=-log10(DE_res$P.Value[sign.genes])+0.2, label=DE_res$Gene_name[sign.genes], cex=0.5)

#### perform GSEA using fgsea function on t-statistic ####

metadata97 <- read.table("../0.Data/RNA_seq/counts_MCS/190123_EnsemblID_genes_GRCh38p12_ensembl97_metadata_biomart.txt", header=TRUE, sep="\t")
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
human_TFT_GTRD <- msigdbr(species="Homo sapiens", category="C3", subcategory="TFT:GTRD")
human_TFT_leg <- msigdbr(species="Homo sapiens", category="C3", subcategory="TFT:TFT_Legacy")
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


# TFT:GTRD #

TFT_GTRD <-  human_TFT_GTRD %>% split(x = .$gene_symbol, f = .$gs_name)

TFT_GTRD_overlap <- c()
for (i in 1:length(TFT_GTRD)){
  TFT_GTRD_overlap[[i]] <- TFT_GTRD[[i]][which(TFT_GTRD[[i]] %in% DE_res97$Gene_name97)]
}
names(TFT_GTRD_overlap) <- names(TFT_GTRD)

TFT_GTRD_tmp <- sapply(TFT_GTRD_overlap, function(i) length(i) > 5 & length(i) < 500)
TFT_GTRD_final <- TFT_GTRD_overlap[TFT_GTRD_tmp]


# TFT_Legacy #

TFT_leg <-  human_TFT_leg %>% split(x = .$gene_symbol, f = .$gs_name)

TFT_leg_overlap <- c()
for (i in 1:length(TFT_leg)){
  TFT_leg_overlap[[i]] <- TFT_leg[[i]][which(TFT_leg[[i]] %in% DE_res97$Gene_name97)]
}
names(TFT_leg_overlap) <- names(TFT_leg)

TFT_leg_tmp <- sapply(TFT_leg_overlap, function(i) length(i) > 5 & length(i) < 500)
TFT_leg_final <- TFT_leg_overlap[TFT_leg_tmp]


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

genesets <- do.call(c, list(Reactome_final, TFT_GTRD_final, TFT_leg_final, GO_final)) # 9881

GSEA_res <- fgsea(genesets, ranks_symbol, minSize=5, maxSize = 500, nperm=1000000)
GSEA_sig <- subset(GSEA_res, GSEA_res$padj < 0.05)
GSEA_sig <- GSEA_sig[order(GSEA_sig$NES, decreasing=TRUE), ] # 289 significant

GSEA_LEgenes <- c()
for (i in 1:nrow(GSEA_sig)){
   GSEA_LEgenes[[i]] <- unnest(GSEA_sig[i,8])$leadingEdge
}
names(GSEA_LEgenes) <- GSEA_sig$pathway



## make bubble plot with ggplot2 ##

library(dplyr)
library(ggplot2)
library(forcats)

# subset reactome pathways and BP GO terms #

human_GOBP <- msigdbr(species="Homo sapiens", category="C5", subcategory="BP")

geneset_finalBP <- GSEA_sig[which(GSEA_sig$pathway %in% human_GOBP$gs_name),] # 67
geneset_finalRe <- GSEA_sig[grepl("^REACTOME_", GSEA_sig$pathway),] # 43
geneset_final_tmp <- rbind.data.frame(geneset_finalBP, geneset_finalRe) # 93

# top upregulated #

geneset_final <- geneset_final_tmp[order(geneset_final_tmp$NES, decreasing=TRUE),]
geneset_final_top10 <- geneset_final[c(1:20,((nrow(geneset_final)-19):nrow(geneset_final))),]

geneset_final2_top10 <- cbind.data.frame(geneset_final_top10,logp = -log10(geneset_final_top10$padj))

geneset_final2_top10$pathway <- sub("REACTOME_", "", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("GO_", "", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("_", " ", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- tolower(geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub(" nmd", "", geneset_final2_top10$pathway)  
geneset_final2_top10$pathway <- sub(" ejc", "", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("mrna", "mRNA", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("nonsense mediated", "nonsense-mediated", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("nuclear trans", "nuclear-trans", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("catabolic process", "catabolic process,", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("srp ", "SRP-", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("dna", "DNA", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("dna", "DNA", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("tca", "TCA", geneset_final2_top10$pathway)
geneset_final2_top10$pathway <- sub("atp", "ATP", geneset_final2_top10$pathway)

ggplot(geneset_final2_top10, aes(x=pathway, y=NES, fill=logp)) +
  geom_bar(stat="identity") + coord_flip() + aes(x = fct_reorder(pathway, NES, .desc=F)) + 
  scale_fill_gradient(low="blue", high="red", name="-log(p-value)") + 
  ylab("Normalized Enrichment Score (NES)") + xlab("") + theme(axis.title.x=element_text(size=12, face="bold"), axis.text=element_text(size=12), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12)) +
  guides(size=guide_legend(order=1)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size=15, face="bold")) + ggtitle("")



### HPO ###

geneset_finalHPO_tmp <- GSEA_sig[grepl("^HP_", GSEA_sig$pathway),] # 103

# top upregulated #

geneset_finalHPO <- geneset_finalHPO_tmp[order(geneset_finalHPO_tmp$NES, decreasing=TRUE),]
geneset_finalHPO_top10 <- geneset_finalHPO[c(1:20,((nrow(geneset_finalHPO)-19):nrow(geneset_finalHPO))),]

geneset_finalHPO2_top10 <- cbind.data.frame(geneset_finalHPO_top10,logp = -log10(geneset_finalHPO_top10$padj))

geneset_finalHPO2_top10$pathway <- sub("HP_", "", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- sub("_", " ", geneset_finalHPO2_top10$pathway)
geneset_finalHPO2_top10$pathway <- tolower(geneset_finalHPO2_top10$pathway)

ggplot(geneset_finalHPO2_top10, aes(x=pathway, y=NES, fill=logp)) +
  geom_bar(stat="identity") + coord_flip() + aes(x = fct_reorder(pathway, NES, .desc=F)) + 
  scale_fill_gradient(low="blue", high="red", name="-log(p-value)") + 
  ylab("Normalized Enrichment Score (NES)") + xlab("") + theme(axis.title.x=element_text(size=12, face="bold"), axis.text=element_text(size=12), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12)) +
  guides(size=guide_legend(order=1)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size=15, face="bold")) + ggtitle("")



### TFT ###

geneset_finalTFT_tmp <- GSEA_sig[!grepl("^HP_", GSEA_sig$pathway),] # 103
geneset_finalTFT_tmp <- geneset_finalTFT_tmp[!grepl("^GO_", geneset_finalTFT_tmp$pathway),] # 103
geneset_finalTFT_tmp <- geneset_finalTFT_tmp[!grepl("^REACTOME_", geneset_finalTFT_tmp$pathway),] # 103

# top upregulated #

geneset_finalTFT_top10 <- geneset_finalTFT_tmp[order(geneset_finalTFT_tmp$NES, decreasing=TRUE),]

geneset_finalTFT2_top10 <- cbind.data.frame(geneset_finalTFT_top10,logp = -log10(geneset_finalTFT_top10$padj))

geneset_finalTFT2_top10$pathway <- sub("_", "", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- sub("_", " ", geneset_finalTFT2_top10$pathway)
geneset_finalTFT2_top10$pathway <- tolower(geneset_finalTFT2_top10$pathway)

ggplot(geneset_finalTFT2_top10, aes(x=pathway, y=NES, fill=logp)) +
  geom_bar(stat="identity") + coord_flip() + aes(x = fct_reorder(pathway, NES, .desc=F)) + 
  scale_fill_gradient(low="blue", high="red", name="-log(p-value)") + 
  ylab("Normalized Enrichment Score (NES)") + xlab("") + theme(axis.title.x=element_text(size=12, face="bold"), axis.text=element_text(size=12), legend.title=element_text(size=12, face="bold"), legend.text=element_text(size=12)) +
  guides(size=guide_legend(order=1)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5, size=15, face="bold")) + ggtitle("")


# make circle plot #
  
### circle plot TFT genes ###
  
GSEA_sig2 <- GSEA_sig[order(GSEA_sig$NES, decreasing=TRUE),]
  
GSEA_sig3 <- GSEA_sig2[!grepl("^GO_", GSEA_sig2$pathway),]
GSEA_sig3 <- GSEA_sig3[!grepl("^REACTOME_", GSEA_sig3$pathway),]
GSEA_sig3 <- GSEA_sig3[!grepl("^HP_", GSEA_sig3$pathway),]

GSEA_sig3 <- subset(GSEA_sig3, GSEA_sig3$pathway != "AGCYRWTTC_UNKNOWN" & GSEA_sig3$pathway != "EVI1_02" & GSEA_sig3$pathway != "EVI1_05" & GSEA_sig3$pathway != "GATA_C" & GSEA_sig3$pathway != "GATA1_02" & GSEA_sig3$pathway != "GATA6" &
         GSEA_sig3$pathway != "GGCNNMSMYNTTG_UNKNOWN" & GSEA_sig3$pathway != "OCT1_04" & GSEA_sig3$pathway != "SMTTTTGT_UNKNOWN" & GSEA_sig3$pathway != "TGANNYRGCA_TCF11MAFG_01" & GSEA_sig3$pathway != "TTAYRTAA_E4BP4_01" & 
         GSEA_sig3$pathway != "YCATTAA_UNKNOWN")

GSEA_sig_new <- cbind.data.frame(Term = GSEA_sig3$pathway, adj_pval = GSEA_sig3$padj)

GSEA_sig_new$Term <- gsub("^GO_", "",GSEA_sig_new$Term)
GSEA_sig_new$Term <- gsub("^REACTOME_", "",GSEA_sig_new$Term)
GSEA_sig_new$Term <- tolower(GSEA_sig_new$Term)
GSEA_sig_new$Term <- gsub("_", " ", GSEA_sig_new$Term)
  
  
GSEA_LEgenes_TFT <- c()
for (i in 1:23){
  GSEA_LEgenes_TFT_tmp <- unnest(GSEA_sig3[i,8])$leadingEdge
  LFC <- subset(DE_res, DE_res$logFC > 0.5 | DE_res$logFC < -0.5)
  GSEA_LEgenes_TFT[[i]] <- GSEA_LEgenes_TFT_tmp[which(GSEA_LEgenes_TFT_tmp %in% LFC$Gene_name)]
}
names(GSEA_LEgenes_TFT) <- GSEA_sig3$pathway
  
GSEA_LEgenes_TFT2 <- do.call("rbind", lapply(GSEA_LEgenes_TFT, as.data.frame))
  
  
  
GSEA_sig_new2 <- vector()
  
GSEA_sig_new2 <- data.frame(GSEA_sig_new[1,], Genes = paste(GSEA_LEgenes_TFT[[1]], collapse=", "))
for (i in 2:23){
  GSEA_sig_new2_tmp <- data.frame(GSEA_sig_new[i,], Genes = paste(GSEA_LEgenes_TFT[[i]], collapse=", "))
  GSEA_sig_new2 <- rbind.data.frame(GSEA_sig_new2, GSEA_sig_new2_tmp)
}
colnames(GSEA_sig_new2) <- c("Term", "adj_pval", "Genes")  
  
GSEA_sig_new2 <- cbind.data.frame(category = paste(rep("TFT", 23)), GSEA_sig_new2)
  
  
LEgenes <- unique(GSEA_LEgenes_TFT2$`X[[i]]`)
DE_LEgenes <- DE_res[which(DE_res$Gene_name %in% LEgenes),c("Gene_name","logFC")]
colnames(DE_LEgenes)[1] <- c("ID")
  
DE_res2 <- cbind.data.frame(ID = DE_res$Gene_name, logFC = DE_res$logFC)
  

circle <- list()
circle$david <- GSEA_sig_new2$Genes
circle$genelist <- DE_LEgenes
circle$process <- GSEA_sig_new2$Term
circle$genes <- DE_res2
  
circ <- circle_dat(GSEA_sig_new2, DE_LEgenes)
chord <- chord_dat(circ, circle$genes, circle$process)
  
col <- c(brewer.pal(n=11, "Spectral"), brewer.pal(n=9, "Reds"), brewer.pal(n=9, "Blues")[c(5,7,9)])
  
GOCluster(circ, circle$process, clust.by="term", term.width=2, term.col=col)
  
GOCluster(circ, circle$process, clust.by="logFC", term.width=2, term.col = col)
# col <- c("#3b0000", "#4e0000", "#620000", "#760000", "#890000", "#9d0000", "#b10000", "#c40000", "#d80000", "#eb0000", "#ff0000", "#ff1414", "#ff2727", "#ff3b3b", "#ff4e4e", "#ff6262" , "#ff7676", "#ff8989", "#ff9d9d", "#ffb1b1", 
         "#000076", "#0000ff", "#8989ff")

  
  
  
  
  
### select gene sets with genes strong logFC ###
  
DE_LEgenes2 <- subset(DE_LEgenes, DE_LEgenes$logFC > 1 | DE_LEgenes$logFC < -1)


GSEA_LEgenes_GO3 <- cbind.data.frame(geneset = rownames(GSEA_LEgenes_GO2), gene = GSEA_LEgenes_GO2$`X[[i]]`)
GSEA_LEgenes_GO3_lfc <- GSEA_LEgenes_GO3[which(GSEA_LEgenes_GO3$gene %in% DE_LEgenes2$ID),]

GSEA_LEgenes_GO3_lfc$geneset <- gsub("\\..*", "", GSEA_LEgenes_GO3_lfc$geneset)

# genesets of interest #

genesets_lfc <- names(table(GSEA_LEgenes_GO3_lfc$geneset))

  
  
plot(DE_res$logFC, -log10(DE_res$P.Value), pch=20, ylab="-log10(p-value)", xlab="LogFC", cex=0.6, xlim=c(-8,8))
ind <- DE_res$adj.P.Val < 0.05 & DE_res$logFC >0
ind2 <- DE_res$adj.P.Val < 0.05 & DE_res$logFC < 0
points(DE_res$logFC[ind], -log10(DE_res$P.Value)[ind], pch = 20, col = "#A50026", cex=0.8, xlim=c(-8,8))
points(DE_res$logFC[ind2], -log10(DE_res$P.Value)[ind2], pch = 20, col = "#313695", cex=0.8, xlim=c(-8,8))

  
## plot enrichment ##

# plotEnrichment(genesets[["GO_AXON_DEVELOPMENT"]], sort(ranks_symbol)) + labs(title="Axon development") + theme(text = element_text(size=15))
# plotEnrichment(genesets[["GO_SYNAPSE_ORGANIZATION"]], sort(ranks_symbol)) + labs(title="Synapse organization") + theme(text = element_text(size=15))
# plotEnrichment(genesets[["REACTOME_SIGNALING_BY_RHO_GTPASES"]], sort(ranks_symbol)) + labs(title="Signaling by Rho GTPases") + theme(text = element_text(size=15))




# sample distances #

rownames(voom_cor) <- voom$genes$Gene_stable_ID_version

sampleDists <- dist(t(voom_cor[,10:21]))

sampleDistMatrix <- as.matrix(sampleDists )
rownames(sampleDistMatrix) <- cell_line[10:21]
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)





