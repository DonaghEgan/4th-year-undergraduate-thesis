# imported libraries 
library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("pheatmap")
library("RColorBrewer")
library("tm")
library("tidyverse")
library("survminer")
library("survival")
library("dplyr")
library("rpart")
library("devtools")
library("gage")
library("pathview")
library("ggplot2")
library("EBcoexpress")
library("clusterProfiler")
library("rWikiPathways")
library("rlang")
library(rpart)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library("reshape2")
library("edgeR")
library("limma")
library("apeglm")
library("VennDiagram")
library("EnhancedVolcano")
library("RAM")
library("enrichplot")
library("org.Hs.eg.db")
library("fgsea")
library("foreign")
library("dafs")
library("stats")
library("matrixStats")
library("ggsignif")
library("plotrix")

# raw count files and clinical data (both in same file).
load("C:/Users/ThinkPad/Input_Files/ucec_raw_counts.RData") 
new_clincial_data <- ucec_cdr_immune_clin

#finding expression values of crabp2 and appending to clinical data. 
row_number <- which(ucec_raw_counts$external_gene_name=="CRABP2",arr.ind = FALSE)
CRABP2_values <- ucec_raw_counts[row_number,3:dim(ucec_raw_counts)[2]]
expression_levels <- as.numeric(as.vector(CRABP2_values[1,]))
CRABP2_values_remodel <- as.data.frame(t(CRABP2_values))
new_clincial_data <- left_join(new_clincial_data, data.frame("CRABP2_counts" = CRABP2_values_remodel$`8331`, "sample" = rownames(CRABP2_values_remodel)))
new_clincial_data$LOG_CRABP2 <- log2(new_clincial_data$CRABP2_counts)

#removing 0s from os.time (for survival analysis) in clinical data by creating a new data set 
new_clinical_OS <- new_clincial_data[new_clincial_data$OS.time != 0 & !is.na(new_clincial_data$OS.time),]

# tree analysis to determine CRABP2 threshold.
fit_tree <- rpart(OS ~ LOG_CRABP2, data= new_clinical_OS, method="anova") # could lead to over-fitting. Ask Bruce about testing model on another data set.
fancyRpartPlot(fit_tree)  # graph showing how patients are dichotomised 
decision_values <- fit_tree$splits
cut_off <- decision_values[1,4]

# grouping patients into high and low cohorts and appending to data frame.
high_low <- ifelse(log2(CRABP2_values) >= cut_off,"High","Low") 
new_clincial_data <- left_join(new_clincial_data, data.frame("CRABP2_val" = as.vector(high_low), "sample" = colnames(high_low)))

# editing clinical_stage in ucec_cdr_immune_clin so there is less bins
v1 <- c("I", "II", "III", "IV")
new_clincial_data$clinical_stage <- gsub(paste0("[^", paste(v1, collapse = ""), "]+"), "", new_clincial_data$clinical_stage)
staging <- new_clincial_data$clinical_stage
make_numeric <- factor(staging, levels=c("I","II","III","IV"))
make_numeric <- as.numeric(make_numeric)-2
new_clincial_data$grouped_staging <- ifelse(make_numeric == -1 | 0 , "1&2", "3&4")

#grouping grades
grading <- new_clincial_data$histological_grade
assign_numeric <- factor(grading, levels=c("G1","G2","G3","High Grade"))
assign_numeric <- as.numeric(assign_numeric)-2
new_clincial_data$grouped_grades <- ifelse(assign_numeric == -1 |assign_numeric == 0, "G1&G2", "G3")

# creating endometrioid and non-endometrioid 
endometrioid <- new_clincial_data$histological_type
endometrioid_numeric <- factor(endometrioid, levels=c("Endometrioid endometrial adenocarcinoma","Serous endometrial adenocarcinoma","Mixed serous and endometrioid"))
endometrioid_numeric <- as.numeric(endometrioid_numeric)-2
new_clincial_data$Endometrioid <- ifelse(endometrioid_numeric == -1, "Endometrioid", "Non_Endometrioid")

# dealing with NA in age (turning into NAs into medians)
median_age <- median(as.numeric(new_clincial_data$age_at_initial_pathologic_diagnosis), na.rm = TRUE)
new_clincial_data$age_at_initial_pathologic_diagnosis <- replace(new_clincial_data$age_at_initial_pathologic_diagnosis, is.na(new_clincial_data$age_at_initial_pathologic_diagnosis), median_age)

# adding variabls to survival data frame (survival analysis)
new_clinical_OS <- new_clincial_data[new_clincial_data$OS.time != 0 & !is.na(new_clincial_data$OS.time),]

# Survival analysis using the Kaplan-Meier method for CRABP2
surv_object <- Surv(time = new_clinical_OS$OS.time, event = new_clinical_OS$OS)
fit1 <- survfit(surv_object ~ CRABP2_val, data = new_clinical_OS)
survdiff(surv_object ~ CRABP2_val, data = new_clinical_OS)  #log rank test
ggsurvplot(fit1, data = new_clinical_OS, pval = TRUE,
           title = "Survival Curve (CRABP2 High vs CRABP2 Low)",
           legend = "bottom",
           xlab = "Time(Days)",
           ylab = "Survival Probability",
           legend.title = "CRABP2 expression",
           legend.labs = c("High Expression (N=146)", "Low Expression (N=380)"),
           font.x = c(15, "plain", "black"),
           font.y = c(15, "plain", "black"),
           pval.size = 5,
           font.legend = c(10, "plain", "black"))
ggsave("CRABP2_curve.png", path = "C:/Users/ThinkPad/output_results", dpi = 800, width =12, height = 9)

# multivariate cox regression overall survival analysis (including confounding variables)
new_clinical_OS$CRABP2_val = relevel(new_clinical_OS$CRABP2_val, ref = "Low")

multi_cox <- coxph(Surv(OS.time, OS) ~ CRABP2_val + grouped_staging + grouped_grades + age_at_initial_pathologic_diagnosis + Endometrioid,
                   data=new_clinical_OS)

# relationship with clinicopathological features for high and low CRABP2 (Pearson Chi-squared test)
mean_BMI <- mean(new_clincial_data$BMI, na.rm = TRUE)
median_age_CRABP2 <- tapply(new_clincial_data$age_at_initial_pathologic_diagnosis, new_clincial_data$CRABP2_val, median, na.rm=TRUE)
range_age_CRABP2 <- tapply(new_clincial_data$age_at_initial_pathologic_diagnosis, new_clincial_data$CRABP2_val, range, na.rm=TRUE)
BMI_mean <- tapply(new_clincial_data$BMI, new_clincial_data$CRABP2_val, mean, na.rm=TRUE)
BMI_SEM <- tapply(new_clincial_data$BMI, new_clincial_data$CRABP2_val, std.error, na.rm=TRUE)
table_histology <- table(new_clincial_data$CRABP2_val,new_clincial_data$histological_type)
hist_chi_test <- chisq.test(table_histology)
table_stage <- table(new_clincial_data$CRABP2_val,new_clincial_data$clinical_stage)
stage_chi_test <- chisq.test(table_stage)
new_clincial_data$histological_grade <- gsub("High Grade", "G3", new_clincial_data$histological_grade)
table_grade <- table(new_clincial_data$CRABP2_val,new_clincial_data$histological_grade)
grade_chi_test <- chisq.test(table_grade)
table_tcga_subtype <- table(new_clincial_data$CRABP2_val, new_clincial_data$`TCGA Subtype`)
table_tcga_subtype <- table_tcga_subtype[ , (3:6)]
tcga_chi_test <- chisq.test(table_tcga_subtype)

#######################DESeq2#################################################################
new_clincial_data$grouped_grades <- ifelse(assign_numeric == -1 | 0, "G1&G2", "G3")  # grades 2&3 grouped togther (addressed in limitations)

# NOTE: Before carrying out differential expression, check options('contrasts'). Should equal "contr.treatment"
# Creating DESeqDataSet 
count_matrix <- as.data.frame(lapply(ucec_raw_counts[,3:dim(ucec_raw_counts)[2]],as.integer))
rownames(count_matrix) <- ucec_raw_counts[,1] 
dds_crabp2 <- DESeqDataSetFromMatrix(countData = count_matrix,
                                     colData = new_clincial_data,
                                     design = ~ grouped_staging + grouped_grades + CRABP2_val)

#Removing no counts
dds_crabp2 <- dds_crabp2[rowSums(counts(dds_crabp2)) > 1,]

#Variance Stabilization
trans_object <- vst(dds_crabp2, blind = FALSE)

#SampleDistances
sampleDists <- dist(t(assay(trans_object)))

# PCA plot 2.0
pcaData <- plotPCA(trans_object, intgroup = c("CRABP2_val"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = CRABP2_val)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave("PCA_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600, width = 6, height = 6)

# Differential Exp Pipeline for DESeq
dds_crabp2 <- DESeq(dds_crabp2)

#extracting results for CRABP2 high vs low
res <- results(dds_crabp2, alpha = .01)
resup <- res[ which(res$padj < .01 & res$log2FoldChange > 1), ]
resdown <- res[ which(res$padj < .01 & res$log2FoldChange < -1), ]
up.genes <- rownames(resup)
down.genes <- rownames(resdown)
bkgd.genes <- rownames(res)

#################edgeR#################################################

# Creating DGElist
genetable <- data.frame(gene.id = rownames(count_matrix),
                        stringsAsFactors = FALSE)
dge <- DGEList(counts = count_matrix, 
               samples = new_clincial_data, 
               genes = genetable)

# Normalization factors
dge <- calcNormFactors(dge)

# Differential Expression design
design_edge_CRABP2 <- model.matrix(~ grouped_staging + grouped_grades + CRABP2_val, data = dge$samples)

# Independent filtering of lowly expressed genes
keep_CRABP2 <- filterByExpr(dge, design_edge_CRABP2)
dge_CRABP2 <- dge[keep_CRABP2, ]
dge_CRABP2 <- estimateDisp(dge_CRABP2, design_edge_CRABP2)

# Carrying out DE for CRABP2 high vs low
fit_CRABP2 <- glmQLFit(dge_CRABP2, design_edge_CRABP2)
qlf_CRABP2 <- glmTreat(fit_CRABP2, coef = "CRABP2_valLow", lfc = 1)
tt.all_CRABP2 <- topTags(qlf_CRABP2, n = Inf, sort.by = "none") # all genes
tt_CRABP2 <- topTags(qlf_CRABP2, n = Inf, adjust.method = "BH", p.value = 0.01) # genes with adj.p<0.1
tt10_CRABP2 <- topTags(qlf_CRABP2) # just the tops 10 by default

#########limma#######limmaa##############################################
new_clincial_data$CRABP2_val <- factor(new_clincial_data$CRABP2_val)
new_clincial_data$grouped_staging <- factor(new_clincial_data$grouped_staging)
new_clincial_data$grouped_grades <- factor(new_clincial_data$grouped_grades)

#Diffierential Expression using Limma. data set creation
rownames(new_clincial_data) <- colnames(as.matrix(count_matrix))
rownames(genetable) <- rownames(as.matrix(count_matrix))
eset <- ExpressionSet(assayData = as.matrix(count_matrix),
                      phenoData = AnnotatedDataFrame(new_clincial_data),
                      featureData = AnnotatedDataFrame(genetable))

# creating design formula
design_limma_CRAP2 <- model.matrix( ~ grouped_staging + grouped_grades + CRABP2_val, data = pData(eset))
#Removing heteroscedascity from count data
voom_CRABP2 <- voom(eset, design_limma_CRAP2, plot=TRUE)

# fitting design with eset
fit_limma_crabp2 <- lmFit(voom_CRABP2, design_limma_CRAP2)
fit_limma_crabp2 <- contrasts.fit(fit_limma_crabp2, coefficients = 4)
fit_limma_crabp2 <- eBayes(fit_limma_crabp2)

# extracting results
tt_res_crabp2l <- topTable(fit_limma_crabp2, number = Inf, adjust.method = "BH", p.value = .01, lfc = 1)
str(tt_res_crabp2l)
# performing tests on results
results_limma_CRABP2 <- decideTests(fit_limma_crabp2, method = "separate", adjust.method = "BH", p.value = 0.01,
                                    lfc = 1)

################ Comparing Packages ########################

# Venn Diagram for CRABP2 packages
Deseq_combined = c(rownames(resdown),rownames(resup))
myCol <- brewer.pal(3, "Set1")
venn.diagram(
  x = list(Deseq_combined, c(rownames(tt_res_crabp2l)), c(rownames(tt_CRABP2))),
  category.names = c("DESeq" , "Limma " , "EdgeR"),
  filename = 'C:/Users/ThinkPad/Output_results/Comparing Packages_CRABP2.PNG', # exporting result
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = myCol
)

######### combining and exporting DEGs #####################

# combining results differential expression results into 1 Data Frame.
storing_p_columns <- as.vector(c("gene.id","adj.P.Val","logFC"))
df_padj_id <- tt_res_crabp2l[,storing_p_columns]

shared <- c(rownames(tt_res_crabp2l), rownames(tt_CRABP2), rownames(resdown),rownames(resup))
genes_in_2 <- shared[ave(seq_along(shared), shared, FUN = length) >= 2]
genes_in_2 <- unique(genes_in_2)
genes_in_2 <- df_padj_id[(df_padj_id$gene.id %in% genes_in_2),]
genes_in_2 <- genes_in_2 %>% arrange((logFC))
genes_in_2_copy <- shared_copy[ave(seq_along(shared_copy), shared_copy, FUN = length) >= 2]
genes_in_2_copy <- unique(genes_in_2_copy)
genes_in_2_copy <- df_padj_id[(df_padj_id$gene.id %in% genes_in_2_copy),]
genes_in_2_copy <- genes_in_2_copy %>% arrange((logFC))
all_inc_edge <-  shared[ave(seq_along(shared), shared, FUN = length) >= 2]
all_inc_edge <- unique(all_inc_edge)
excluded <- c(all_inc_edge, genes_in_2$gene.id)
excluded <- excluded[ave(seq_along(excluded), excluded, FUN = length) == 1]
edge_columns <- as.vector(c("gene.id","FDR", "logFC"))
df_fdr_edge <- tt_CRABP2[,edge_columns]
edge_deseq <- df_fdr_edge[(df_fdr_edge$table$gene.id %in% excluded),]
edge_deseq <- as.data.frame(edge_deseq)
edge_deseq <- edge_deseq %>% arrange((logFC))
colnames(edge_deseq)[2] <- "adj.P.Val"
all_3 <- rbind(edge_deseq, genes_in_2)
all_3 <- all_3 %>% arrange((logFC))

# exporting genes out to carry out ARACNE and GSEA and gene symbol
write.csv(all_3,"C:/Users/ThinkPad/Working_Directory/ranked_CRABP2_3_r.csv", row.names = FALSE)

############## Volcano Plot ###################################

# Note ensembl gene names converted to gene symbol. 5 genes removed (no gene symbol)
volcano_plot_3 <- read.csv("C:/Users/ThinkPad/Working_Directory/ranked_CRABP2_genenames.csv", header = TRUE)
names(volcano_plot_3)[3] <- "Log2FC"
volcano_plot_3$Log2FC <- (volcano_plot_3$Log2FC)*-1
EnhancedVolcano(volcano_plot_3,
                lab = volcano_plot_3$gene.id,
                x = 'Log2FC',
                y = 'adj.P.Val',
                FCcutoff = 1.5,
                pCutoff = 10e-16,
                transcriptPointSize = 3,
                xlim = c(-3, 3),
                ylim = c(0,30))
ggsave("volcano_plot_final.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600,width =9, height = 6)

############## GSEA ###################################

# importing results from GSEA
gsea_neg_res <- read.csv("C:/Users/ThinkPad/Working_Directory/gsea_neg.csv", header = TRUE)
gsea_pos_res <- read.csv("C:/Users/ThinkPad/Working_Directory/gsea_pos.csv", header = TRUE)

# plotting top 30 pathways that are upregulated in CRABP2 cohort. 
gsea_neg_res <- gsea_neg_res[1:30,]
gsea_neg_res$NES <- abs(gsea_neg_res$NES)

# plotting top 30 pathways that are down regualted pathways in CRABP2 cohort.
gsea_pos_res <- gsea_pos_res[1:30,]
gsea_pos_res$NES <- -abs(gsea_pos_res$NES)

# joning upregulated and downregulated
gsea_pos_neg <- rbind(gsea_pos_res[1:10,], gsea_neg_res[1:10,])
gsea_com_plot <- ggplot(gsea_pos_neg, aes(NES, NAME))
gsea_com_plot + geom_point(aes(colour=FDR.q.val, size=SIZE)) +
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.1)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  theme(panel.background=element_rect(fill="gray95", colour="gray95"),
        panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
        panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
        axis.title.y=element_blank()) +
  expand_limits(x=c(-5,5)) +
  scale_x_continuous(breaks=c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +
  scale_y_discrete(limits=rev(gsea_pos_neg$NAME)) 
ggsave("GSEA_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600, width =9, height = 6)

############## ARACNE ###################################

# importing significant results from ARACNE
high_vs_rest_sig <- read_tsv("C:/Users/ThinkPad/Working_Directory/CRAPB2.High_vs_rest.sig.tsv")
Low_vs_rest_sig <-  read_tsv("C:/Users/ThinkPad/Working_Directory/CRAPB2.Low_vs_rest.sig.tsv")

# Finding location in UCEC raw counts
RHCG_loc <- which(ucec_raw_counts$external_gene_name=="RHCG",arr.ind = FALSE)
TMEFF2_loc <- which(ucec_raw_counts$external_gene_name=="TMEFF2",arr.ind = FALSE)

#extracting row of TMEFF2 and RHCG
RHCG_values <- ucec_raw_counts[RHCG_loc,3:dim(ucec_raw_counts)[2]]
TMEFF2_values <-  ucec_raw_counts[TMEFF2_loc,3:dim(ucec_raw_counts)[2]]

# extracting values from row
RHCG_levels <- as.numeric(as.vector(RHCG_values[1,]))
TMEFF2_levels <- as.numeric(as.vector(TMEFF2_values[1,]))

# normalising values log2
log_RHCG <- log2(RHCG_levels)
log_TMEFF2 <- log2(TMEFF2_levels)

#visualising exp data for both genes
hist(log_RHCG)
hist(log_TMEFF2)

# Add genes to clinical data
new_clincial_data$LOG_RHCG <- log_RHCG
new_clincial_data$LOG_TMEFF2 <- log_TMEFF2

# tree analysis to determine RHCG threshold
tree_RHCG <- rpart(OS ~ LOG_RHCG, data= new_clinical_OS, method="anova")
fancyRpartPlot(tree_RHCG)
decision_values_RHCG <-tree_RHCG$splits
cut_off_RHCG <- decision_values_RHCG[1,4]

# tree analysis to determine TMEFF2 threshold
tree_TMEFF2 <- rpart(OS ~ LOG_TMEFF2, data= new_clinical_OS, method="anova")
fancyRpartPlot(tree_TMEFF2)
decision_values_TMEFF2 <-tree_TMEFF2$splits
cut_off_TMEFF2 <- decision_values_TMEFF2[1,4]

# Splitting patients into there associated cohorts
high_low_RHCG <- ifelse(log_RHCG >= cut_off_RHCG,"High","Low") 
new_clincial_data$RHCG_val <- high_low_RHCG
high_low_TMEFF2 <- ifelse(log_TMEFF2 >= cut_off_TMEFF2,"High","Low")
new_clincial_data$TMEFF2_val <- high_low_TMEFF2

# updating survival data frame
new_clinical_OS <- new_clincial_data[new_clincial_data$OS.time != 0 & !is.na(new_clincial_data$OS.time),]

########### kaplan Meier for RHCG and TMEFF2 ####################################

# Fit survival data using the Kaplan-Meier method for RHCG
surv_object_RHCG <- Surv(time = new_clinical_OS$OS.time, event = new_clinical_OS$OS)
fit1_RHCG <- survfit(surv_object_RHCG ~ RHCG_val, data = new_clinical_OS)
survdiff(surv_object_RHCG ~ RHCG_val, data = new_clinical_OS)  #log rank test
ggsurvplot(fit1_RHCG, data = new_clinical_OS, pval = TRUE, palette = "simpsons",
           legend = "bottom",
           xlab = "Time(Days)",
           legend.title = "RHCG expression",
           legend.labs = c("High Expression (N=323)", "Low Expression (N=203)"))
ggsave("RHCG_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600)

# Fit survival data using the Kaplan-Meier method for TMEFF2
surv_object_TMEFF2 <- Surv(time = new_clinical_OS$OS.time, event = new_clinical_OS$OS)
fit1_TMEFF2 <- survfit(surv_object_TMEFF2 ~ TMEFF2_val, data = new_clinical_OS)
survdiff(surv_object_TMEFF2 ~ TMEFF2_val, data = new_clinical_OS)  #log rank test
ggsurvplot(fit1_TMEFF2, data = new_clinical_OS, pval = TRUE, palette = "uchicago",
           legend = "bottom",
           xlab = "Time(Days)",
           legend.title = "TMEFF2 expression",
           legend.labs = c("High Expression (N=38)", "Low Expression (N=488)"))
ggsave("TMEFF2_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600)

################# Creating gene-panel #######################################

# new data frame with genes for gene panel
high_low_columns <- as.vector(c("TMEFF2_val","CRABP2_val", "RHCG_val"))                                   
high_low_df <- new_clinical_OS[,high_low_columns]

# assigning numeric value (1=high, 0=low)
to_number_TMEFF2 <- as.numeric(factor(high_low_df$TMEFF2_val, levels=c("Low","High")))
to_number_TMEFF2 <- to_number_TMEFF2 -1 
high_low_df$numeric_TMEFF2 <- to_number_TMEFF2

# assigning numeric value (1=high, 0=low)
to_number_RHCG <- as.numeric(factor(high_low_df$RHCG_val, levels=c("Low","High")))
to_number_RHCG <- to_number_RHCG -1
high_low_df$numeric_RHCG <- to_number_RHCG

# assigning numeric value (1=high, 0=low)
to_number_CRABP2 <- as.numeric(factor(high_low_df$CRABP2_val, levels=c("Low","High")))
to_number_CRABP2 <- to_number_CRABP2 -1
high_low_df$numeric_CRABP2 <- to_number_CRABP2

# adding numeric values and adding to gene panel data frame
new_clincial_data$grouped_grades <- ifelse(assign_numeric == -1 |assign_numeric == 0, "G1&G2", "G3")
new_clinical_OS <- new_clincial_data[new_clincial_data$OS.time != 0 & !is.na(new_clincial_data$OS.time),]
high_low_df$sum_of_genes <- rowSums(high_low_df[,c("numeric_TMEFF2", "numeric_RHCG","numeric_CRABP2")])
new_clinical_OS$sum_of_genes <- high_low_df$sum_of_genes

# kaplan meier based on the number of genes that are highly expressed
surv_object_sum <- Surv(time = new_clinical_OS$OS.time, event = new_clinical_OS$OS)
fit1_sum <- survfit(surv_object_sum ~ sum_of_genes, data = new_clinical_OS)
survdiff(surv_object_sum ~ sum_of_genes, data = new_clinical_OS)  #log rank test
ggsurvplot(fit1_sum, data = new_clinical_OS, pval = TRUE,
           legend = "bottom",
           xlab = "Time(Days)",
           legend.title = "Number of genes highly expressed",
           legend.labs = c("0","1","2","3"))
ggsave("sum_of_genes_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600)

#univariate analysis sum of genes 
cox_1_gene <- coxph(Surv(OS.time, OS) ~ as.factor(sum_of_genes) , data = new_clinical_OS)

# multivariate analyis based on sum of genes
cox_sum_genes <- coxph(Surv(OS.time, OS) ~ as.factor(sum_of_genes) + grouped_staging + grouped_grades + age_at_initial_pathologic_diagnosis + Endometrioid,
                       data=new_clinical_OS)

# relationship with clinicopathological features for gene signature
new_clinical_OS$signature_grouped <- as.numeric(gsub(3, 2, new_clinical_OS$sum_of_genes))
median_age_sig <- tapply(new_clinical_OS$age_at_initial_pathologic_diagnosis, new_clinical_OS$signature_grouped, median, na.rm=TRUE)
range_age_sig <- tapply(new_clinical_OS$age_at_initial_pathologic_diagnosis, new_clinical_OS$signature_grouped, range, na.rm=TRUE)
BMI_sig <- tapply(new_clinical_OS$BMI, new_clinical_OS$signature_grouped, mean, na.rm=TRUE)
BMI_SEM_sig <- tapply(new_clinical_OS$BMI, new_clinical_OS$signature_grouped, std.error, na.rm=TRUE)
table_histology_sig <- table(new_clinical_OS$signature_grouped,new_clinical_OS$histological_type)
hist_chi_test_sig <- chisq.test(table_histology_sig)
table_stage_sig <- table(new_clinical_OS$signature_grouped,new_clinical_OS$clinical_stage)
stage_chi_test_sig <- chisq.test(table_stage_sig)
new_clinical_OS$histological_grade <- gsub("High Grade", "G3", new_clinical_OS$histological_grade) 
table_grade_sig <- table(new_clinical_OS$signature_grouped,new_clinical_OS$histological_grade)
grade_chi_test_sig <- chisq.test(table_grade_sig)
table_tcga_subtype_sig <- table(new_clinical_OS$signature_grouped, new_clinical_OS$`TCGA Subtype`)
table_tcga_subtype_sig <- table_tcga_subtype_sig[ , (3:6)]
tcga_chi_test_sig <- chisq.test(table_tcga_subtype_sig)

################# CIBERSORT ##########################################

#Importing cibersort results. Generated using UCEC raw counts
cibersort <- read_tsv("C:/Users/ThinkPad/Working_Directory/CIBERSORT-Results.txt")
# perm = 100 

# appending CRABP2 values onto cibersort dataframe
cibersort$CRABP2_val <- new_clincial_data$CRABP2_val
cibersort <- cibersort %>% arrange(CRABP2_val)
cibersort <- cibersort[cibersort$`P-value` < 0.05, ]

# Changing column names in cibersort dataframe
names(cibersort)<-str_replace_all(names(cibersort), c(" " = "_" , "," = "" , "`" = ""))
colnames(cibersort)[10] <- "T_cells_regulatory"

# wilcoxon 
cibersort_names <- colnames(cibersort)

modelList<-list()
for(i in 2:23){
  fmla <- formula(paste(names(cibersort)[i]," ~ CRABP2_val"))
  modelList[[i]]<-wilcox.test(fmla, data = cibersort, paired = FALSE, alternative="two.sided")
}
capture.output(modelList, file = "C:/Users/ThinkPad/Output_results/Wilcox_test.txt")

# plotting 
cibersort_high <- cibersort[,c(27, 2:11)] # cibersort_high adds crabp2 val to first column and removes p-val etc
colnames(cibersort_high)[7] <- "CD4_memory_resting"
colnames(cibersort_high)[8] <- "CD4_memory_activated"
cibersort_melted <- melt(cibersort_high, id.var = "CRABP2_val")
boxplot_immune_innate <- ggplot(data = cibersort_melted, aes(x=variable, y=value)) + geom_boxplot(aes(fill=CRABP2_val)) +
  facet_wrap( ~ variable, scales="free") +
  xlab("Immune Cell Type") + ylab("Immune Cell Fraction") + ggtitle("Innate Immunity")
ggsave("Innate_Immunity_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600, width =9, height = 6)


cibersort_low <- cibersort[,c(27, 12:23)]
melted_adaptive <- melt(cibersort_low, id.var = "CRABP2_val")
boxplot_immune_adaptive <- ggplot(data = melted_adaptive, aes(x=variable, y=value)) + geom_boxplot(aes(fill=CRABP2_val)) +
  facet_wrap( ~ variable, scales="free") +
  xlab("Immune Cell Type") + ylab("Immune Cell Fraction") + ggtitle("Adaptive Immunity")
ggsave("Adaptive_immunity_plot.png", path = "C:/Users/ThinkPad/Output_results", dpi = 600, width =9, height = 6)

########## boxplots for expression levels of genes. #############
b_CRABP2 <- boxplot(LOG_CRABP2 ~ CRABP2_val,data = new_clincial_data, plot = 0)
boxplot(LOG_CRABP2 ~ CRABP2_val, ylab = "Relative CRABP2 expression (Log2)", xlab = "CRABP2 cohort", col = c('brown2', 'royalblue'), data = new_clincial_data,names=paste(b_CRABP2$names, "(n=", b_CRABP2$n, ")"))


b_TMEFF2 <- boxplot(LOG_TMEFF2 ~ TMEFF2_val,data = new_clincial_data, plot = 0)
boxplot(LOG_TMEFF2 ~ TMEFF2_val, ylab = "Relative TMEFF2 expression (Log2)", xlab = "TMEFF2 cohort", col = c('coral4', 'grey'), data = new_clincial_data,names=paste(b_TMEFF2$names, "(n=", b_TMEFF2$n, ")"))

b_RHCG <- boxplot(LOG_RHCG ~ RHCG_val,data = new_clincial_data, plot = 0)
boxplot(LOG_RHCG ~ RHCG_val, ylab = "Relative RHCG expression (Log2)", xlab = "RHCG cohort", col = c('gold1', 'dodgerblue'), data = new_clincial_data,names=paste(b_RHCG$names, "(n=", b_RHCG$n, ")") )



