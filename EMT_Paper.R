#################
### EMT PAPER ###
#################

  
# Loading libraries
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(shades)
  library(ggpubr)
  library(ggsignif)
  library(gridExtra)
  library(openxlsx)
  
# Loading functions
source("emtPaperFunctions.R")


# Phenotype designations table by study for TCGA and GTEx
abbreviationsTable <- read.table("study_abbreviations.csv", sep = ",", stringsAsFactors = FALSE, 
                                 header = TRUE)
# pdf("Final Figures/abbreviations_table.pdf", width = 10, height = 15.5)
abbTable <- grid.table(abbreviationsTable, rows = NULL)
# dev.off()

# Loading files with scores obtained from running on these methods on the same samples
gsva <- read.table("gsva_emt_scores.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
singscore <- read.table("singscore_emt_scores.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ssgsea <- read.table("ssgsea_emt_scores.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sspaths <- read.table("sspaths_emt_scores_runABS.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Loading sample information
samples <- read.table("all_samples.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Merging all the datasets with the samples df so each has all the information about each sample
gsva <- merge(gsva, samples, by = "sample_id")
singscore <- merge(singscore, samples, by = "sample_id")
ssgsea <- merge(ssgsea, samples, by = "sample_id")
sspaths <- merge(sspaths, samples, by = "sample_id")


##############################
# miRNA Correlation Analysis #
##############################

# Getting miRNA expression for TCGA samples (no expression available for GTEx samples)
mirna <- read.table("TCGA_miRNA_Expression.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)[, -c(2, 4)]

# Setting the TCGA IDs in mirna to match methodDFs 
tumorSamps <- gsva[gsva$Is_Normal == FALSE, "sample_id"]
tumorSamps <- strtrim(tumorSamps, 16)
# Selecting for common samples in the two dataframes
mirna <- mirna[mirna$sample_id %in% tumorSamps, ]
# Scaling by sample across all (1881) avalible mirna expression
mirna[, 3:ncol(mirna)] = t(scale(t(mirna[, 3:ncol(mirna)]))) 

# Selecting for important miRNA only
impMirna = c("hsa.mir.200a", "hsa.mir.200b", "hsa.mir.200c", "hsa.mir.429", "hsa.mir.141", "hsa.mir.214",
             "hsa.mir.199a.1", "hsa.mir.199a.2", "hsa.mir.199b")
mirna = mirna[, c(colnames(mirna)[1:2], impMirna)]

# Running a function that combines the score and miRNA expression dfs and creates a linear model 
# for tumor samples. Look at emtPaperFunctions.R for details
gsvaCorr <- lmScoreByMirna(methodDF = gsva[gsva$Is_Normal == FALSE, 1:2], mirnaDF = mirna, 
                           cancerType = TRUE, method = "GSVA")
singscoreCorr <- lmScoreByMirna(methodDF = singscore[singscore$Is_Normal ==  FALSE, 1:2], mirnaDF = mirna, 
                                cancerType = TRUE, method = "singscore")
ssgsesCorr <- lmScoreByMirna(methodDF = ssgsea[ssgsea$Is_Normal == FALSE, 1:2], mirnaDF = mirna, 
                             cancerType = TRUE, method = "ssGSEA")
sspathsCorr <- lmScoreByMirna(methodDF = sspaths[sspaths$Is_Normal == FALSE, 1:2], mirnaDF = mirna, 
                              cancerType = TRUE, method = "ssPATHS")

# Combine the results
allCorr <- rbind(gsvaCorr, singscoreCorr, ssgsesCorr, sspathsCorr)
# FDR correction on p-values
allCorr$pAdj <- p.adjust(allCorr$pValue, method = "fdr")
# Label using HGNC symbols
mirnaHGNC <- read.table("mirna_HGNC.tsv", header = TRUE, sep = "\t")
allCorr <- merge(allCorr, mirnaHGNC, by = "miRNA")
allCorr$miRNA_HGNC <- factor(allCorr$miRNA_HGNC, levels = c("MIR200A", "MIR200B", "MIR200C", "MIR429", "MIR141",
                                                            "MIR199A1", "MIR199A2", "MIR199B", "MIR214"))
# Creating a column that associates the miRNA with a phenotype (epithelial or mesenchymal)
allCorr$Phenotype <- ifelse(allCorr$miRNA_HGNC %in% c("MIR199A1", "MIR199A2", "MIR199B", "MIR214"), 
                            "Mesenchymal", "Epithelial") 
# Creating a column to show significance 
allCorr$Significance <- signifStars(allCorr$pAdj)
allCorr$Significance <- factor(allCorr$Significance, levels = c("ns", "*", "**", "***", "****"))

# Creating the forest plot with all important miRNAs
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")
# pdf("Final Figures/mirna_all_mirnas_tumor.pdf", width = 16, height = 6)
mirnaForest = ggplot(data = allCorr,
                     aes(x = Study_Abbr, y = Estimate, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col = Significance), size = 0.25)+ theme_bw() +
  geom_hline(aes(fill = Significance),yintercept =0, linetype=2)+
  xlab('TCGA Tumor Study')+ ylab("Coefficient Estimate")+ labs(color = "Significance", fill = "Significance") +
  geom_errorbar(aes(ymin = lower, ymax = upper,col=Significance),width=0.5,cex=0.5)+ 
  facet_grid(Method~miRNA_HGNC, scales = "free") + scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + coord_flip()

# Adding colors to the labels of the forest plot
blueSat = saturation(c("#0072B2"), 0.8)
colors = brightness(c("#D55E00", "#009E59", blueSat), 0.4)
scoreColors = brightness(colors, 0.85)
scoreColors = opacity(scoreColors, 0.4)

gp <- ggplotGrob(mirnaForest)

# The two for loops below were obtained from 
# https://stackoverflow.com/questions/52668135/how-to-color-facet-grid-by-group-in-ggplot2?rq=1
for(i in 1:5){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- scoreColors[1]
}
# Assign the next 4 right-side facet strips with red fill
for(i in 6:9){
  grob.i <- grep("strip-t", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- scoreColors[3]
}

grid::grid.draw(gp)

# dev.off()

# Saving the dataframe used to make the plot as a .tsv file
# write.table(allCorr[, -c(1, 9)], "Plot Creation Files/all_miRNA_forest_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


# Selecting only a two miRNAs for the main figure. Formating as above
# pdf("Final Figures/mirna_forest_plot_tumor.pdf", width = 7, height = 3)
selectedMirnas <- c("MIR429", "MIR214")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00")
mirnaForest = ggplot(data=allCorr[allCorr$miRNA_HGNC %in% selectedMirnas, ],
                     aes(x = Study_Abbr, y = Estimate, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col = Significance), size = 0.25)+ theme_bw() +
  geom_hline(aes(fill = Significance),yintercept =0, linetype=2)+
  xlab('TCGA Tumor Study')+ ylab("Coefficient Estimate")+ labs(color = "Significance", fill = "Significance") +
  geom_errorbar(aes(ymin = lower, ymax = upper,col=Significance),width=0.5,cex=0.5)+ 
  facet_grid(miRNA_HGNC~Method, scales = "free") + scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) + coord_flip()

blueSat = saturation(c("#0072B2"), 0.8)
colors = brightness(c("#D55E00", "#009E59", blueSat), 0.4)
scoreColors = brightness(colors, 0.85)
scoreColors = opacity(scoreColors, 0.4)

gp <- ggplotGrob(mirnaForest)

# The two for loops below were obtained from 
# https://stackoverflow.com/questions/52668135/how-to-color-facet-grid-by-group-in-ggplot2?rq=1
for(i in 1){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- scoreColors[1]
}
# assign the next 4 right-side facet strips with red fill
for(i in 2){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- scoreColors[3]
}

grid::grid.draw(gp)
# dev.off()

# Saving the dataframe used to make the plot as a .tsv file
# write.table(allCorr[allCorr$miRNA_HGNC %in% selectedMirnas, -c(1, 9)],
#             "Plot Creation Files/two_miRNA_forest_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


# Creating a table with the p-values of the miRNAs not considering cancer type (CT). Funtion details in 
# emtPaperFunctions.R
gsvaCorr <- lmScoreByMirna(methodDF = gsva[gsva$Is_Normal == FALSE, 1:2], mirnaDF = mirna , 
                           cancerType = FALSE, method = "GSVA")
singscoreCorr <- lmScoreByMirna(methodDF = singscore[singscore$Is_Normal == FALSE, 1:2], mirnaDF = mirna , 
                                cancerType = FALSE, method = "singscore")
ssgsesCorr <- lmScoreByMirna(methodDF = ssgsea[ssgsea$Is_Normal == FALSE, 1:2], mirnaDF = mirna , 
                             cancerType = FALSE, method = "ssGSEA")
sspathsCorr <- lmScoreByMirna(methodDF = sspaths[sspaths$Is_Normal == FALSE, 1:2], mirnaDF = mirna , 
                              cancerType = FALSE, method = "ssPATHS")

# Combine the results
allCorr <- rbind(gsvaCorr, singscoreCorr, ssgsesCorr, sspathsCorr)

# FDR correction on p-values
allCorr$pAdj <- p.adjust(allCorr$pValue, method = "fdr")

# Label using HGNC symbols
mirnaHGNC <- read.table("mirna_HGNC.tsv", header = TRUE, sep = "\t")
allCorr <- merge(allCorr, mirnaHGNC, by = "miRNA")
allCorr$miRNA_HGNC <- factor(allCorr$miRNA_HGNC, levels = c("MIR200A", "MIR200B", "MIR200C", "MIR429", "MIR141",
                                                            "MIR199A1", "MIR199A2", "MIR199B", "MIR214"))
# Creating a column that associates the miRNA with a phenotype (epithelial or mesenchymal)
allCorr$Phenotype <- ifelse(allCorr$miRNA_HGNC %in% c("MIR199A1", "MIR199A2", "MIR199B", "MIR214"), 
                            "Mesenchymal", "Epithelial") 
allCorr$Signif_Symbol <- signifStars(allCorr$pAdj)

# Order columns by miRNA and method
allCorr <- allCorr[order(allCorr[, "miRNA_HGNC"], allCorr[, "Method"]), ]

# pdf("Final Figures/miRNA_table.pdf", width = 12, height = 16)
tableCorr <- allCorr[, c("Phenotype", "miRNA_HGNC", "Estimate", "pValue", "pAdj", "Signif_Symbol", "Method")]
colnames(tableCorr)[3] <- c("Coeff_Estimate")
mirnaTable <- grid.table(tableCorr,  rows = NULL)
mirnaTable
# write.xlsx(tableCorr, "Supplementary Tables/mirna_method_score_correlation_table.xlsx")
# dev.off()


##########################
# Uterine Cancer Samples #
##########################
# Loading the samples 
uterineIDs <- read.table("Uterine_Cancers_TCGA_IDs.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
# Removing TCGA adj samples and selecting for TCGA tumor samples
uterineIDs <- uterineIDs[uterineIDs$Is_Normal == FALSE, ]
# Selecting for uterine samples in the method DFs
uGsva <- gsva[gsva$sample_id %in% uterineIDs$TCGA_ID, ]
uGsva$Method <- "GSVA"
uSingscore <- singscore[singscore$sample_id %in% uterineIDs$TCGA_ID, ]
uSingscore$Method <- "singscore"
uSsgsea <- ssgsea[ssgsea$sample_id %in% uterineIDs$TCGA_ID, ]
uSsgsea$Method <- "ssGSEA"
uSspaths <- sspaths[sspaths$sample_id %in% uterineIDs$TCGA_ID, ]
uSspaths$Method <- "ssPATHS"

# Combining all data frames
uterusDF <- rbind(uGsva, uSingscore, uSsgsea, uSspaths)
# Setting the order of the study levels
uterusDF$Study_Abbr <- factor(uterusDF$Study_Abbr, levels = c("UCEC", "UCS", "SARC"))

# Creating the desired dark colors
# boxplotColors = brightness(c("#F8766D", "#00BA38", "#619CFF"), 0.4)
# Colorblind palette
boxplotColors = brightness(c("#D55E00", "#009E73", "#0072B2"), 0.4)

# Plotting boxplots
# pdf("Final Figures/uterine_boxplots.pdf", width = 11, height = 10)
ggplot(uterusDF, aes(x = Study_Abbr, y = emt, fill = Tissue_Type, color = Tissue_Type)) +
  geom_boxplot(lwd = 0.8) + theme_bw()+ scale_fill_manual(values = c("#D55E00", "#009E73", "#0072B2")) +
  geom_jitter(aes(color = uterusDF$Tissue_Type, alpha = 0.2), width = 0.25) + xlab("TCGA Tumor Study") +
  ylab("Score") + labs(color = "Tissue\nType", fill = "Tissue\nType") + guides(alpha = FALSE) +
  theme(text = element_text(size = 20)) + facet_wrap(~Method, scales = "free") +
  geom_signif(comparisons = list(c("UCEC", "UCS"), c("UCS", "SARC"), c("UCEC", "SARC")),
              color = "black", map_signif_level = TRUE, step_increase = 0.2)
# dev.off()

# Saving the dataframe used to make the plot as a .tsv file
# write.table(uterusDF, "Plot Creation Files/uterine_boxplots.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


# Uterine cancer waterfall plots
# Getting the expression of genes of interest (GOI) for uterine TCGA samples
exprDF <- read.table("Uterine_Cancers_TCGA_GOI_Expr.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Narrowing down exprDF to samples included in analysis only (i.e. that we have scores for)
exprDF <- exprDF[exprDF$TCGA_ID %in% uGsva$sample_id, ]
colnames(exprDF)[3] <- "sample_id"

# Setting the order of the study levels in each method df so it plots properly below with UCEC as red,
# UCS green and SARC as blue
uGsva$Study_Abbr <- factor(uGsva$Study_Abbr, levels = c("UCEC", "UCS", "SARC"))
uSingscore$Study_Abbr <- factor(uSingscore$Study_Abbr, levels = c("UCEC", "UCS", "SARC"))
uSsgsea$Study_Abbr <- factor(uSsgsea$Study_Abbr, levels = c("UCEC", "UCS", "SARC"))
uSspaths$Study_Abbr <- factor(uSspaths$Study_Abbr, levels = c("UCEC", "UCS", "SARC"))

# Plotting the waterfall plot and maker expression heatmaps. Details of function in emtPaperFunctions.R
# pdf("Final Figures/uterine_waterfall_plots.pdf", width = 12, height = 10)
wGsva <- markersBySample(dfWscore = uGsva[, c(1:4)], dfWexpr = exprDF[, -c(1:2, 4)], 
                         xlabel = "Uterine Samples", ylabel = "GSVA", legendTitle = "Study",
                         scorePlotBreaks = 7, legendVisible = TRUE)
wGsva$Plot
wSingscore <- markersBySample(dfWscore = uSingscore[, c(1:4)], dfWexpr = exprDF[, -c(1:2, 4)],
                              xlabel = "Uterine Samples", ylabel = "singscore", legendTitle = "Study",
                              scorePlotBreaks = 9, markerYaxisSize = 15, legendVisible = TRUE)
wSingscore$Plot
wSsgsea <- markersBySample(dfWscore = uSsgsea[, c(1:4)], dfWexpr = exprDF[, -c(1:2, 4)], 
                            xlabel = "Uterine Samples", ylabel = "ssGSEA", legendTitle = "Study", 
                           markerYaxisSize = 15, legendVisible = TRUE)
wSsgsea$Plot
wSspaths <- markersBySample(dfWscore = uSspaths[, c(1:4)], dfWexpr = exprDF[, -c(1:2, 4)], 
                            xlabel = "Uterine Samples", ylabel = "ssPATHS", legendTitle = "Study",
                            scorePlotBreaks = 7, markerYaxisSize = 16.5, legendVisible = TRUE)
wSspaths$Plot
# dev.off()


# # Writing plot creation files for the uterine marker expression plots
# write.table(wGsva$exprDF, "Plot Creation Files/uterine_gsva_expr_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSingscore$exprDF, "Plot Creation Files/uterine_singscore_expr_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSsgsea$exprDF, "Plot Creation Files/uterine_ssgsea_expr_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSspaths$exprDF, "Plot Creation Files/uterine_sspaths_expr_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)

# # Plot creation files for uterine score waterfall plots
# write.table(wGsva$scoreDF, "Plot Creation Files/uterine_gsva_score_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSingscore$scoreDF, "Plot Creation Files/uterine_singscore_score_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSsgsea$scoreDF, "Plot Creation Files/uterine_ssgsea_score_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(wSspaths$scoreDF, "Plot Creation Files/uterine_sspaths_score_plot.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


# Getting correlation table from marker plot for supplementary data
wCorr <- rbind(wGsva$Corr_DF, wSingscore$Corr_DF, wSsgsea$Corr_DF, wSspaths$Corr_DF)
wCorr$Phenotype <- ifelse(wCorr$HGNC_Symbol %in% c("VIM", "FN1", "CDH2"), "Mesenchymal", "Epithelial")
wCorr <- wCorr[, c("HGNC_Symbol", "Corr_Coeff", "pValues", "pAdj", "Stars", "Phenotype", "Method")]
wCorr <- wCorr[order(wCorr[, "Method"], wCorr[, "Phenotype"]), ]
colnames(wCorr)[5] <- "Signif_Symbol"
# pdf("Final Figures/uterine_correlation_table.pdf", width = 12, height = 16)
wCorrTable <- grid.table(wCorr, rows = NULL)

# Creating an .xlsx file for the correlation of score and marker expression for all uterine cancers
# write.xlsx(wCorr, "Supplementary Tables/uterine_correlation_table.xlsx")
# dev.off()

# Testing correlation only on UCS
dfWexprOnlyUCS = exprDF[exprDF$Study_Abbr == "UCS", -c(1:2, 4)]
wUCSGsva <- markersBySample(dfWscore = uGsva[uGsva$Study_Abbr == "UCS", c(1:4)], dfWexpr = dfWexprOnlyUCS, 
                         ylabel = "GSVA", legendTitle = "Study",
                         scorePlotBreaks = 7, markerGradientTitle = 12, legendVisible = FALSE)
wUCSSingscore <- markersBySample(dfWscore = uSingscore[uSingscore$Study_Abbr == "UCS", c(1:4)], 
                                 dfWexpr = dfWexprOnlyUCS, ylabel = "singscore", legendTitle = "Study",
                              scorePlotBreaks = 9, markerYaxisSize = 15, markerGradientTitle = 12, legendVisible = FALSE)
wUCSSsgsea <- markersBySample(dfWscore = uSsgsea[uSsgsea$Study_Abbr == "UCS", c(1:4)], dfWexpr = dfWexprOnlyUCS, 
                           ylabel = "ssGSEA", legendTitle = "Study", markerYaxisSize = 15, 
                           markerGradientTitle = 12, legendVisible = FALSE)
wUCSSspaths <- markersBySample(dfWscore = uSspaths[uSspaths$Study_Abbr == "UCS", c(1:4)], dfWexpr = dfWexprOnlyUCS, 
                            ylabel = "ssPATHS", legendTitle = "Study", scorePlotBreaks = 7, 
                            markerYaxisSize = 16.5, markerGradientTitle = 12, legendVisible = TRUE)

# Creating a mulitplot of all UCS plots
# pdf("Final Figures/uterine_correlation_waterfall_plots_only_UCS.pdf", width = 30, height = 8)
ucsWaterfall <- ggarrange(wUCSGsva$Plot, wUCSSingscore$Plot, wUCSSsgsea$Plot, wUCSSspaths$Plot, 
                          ncol = 4, widths = c(0.7, 0.7, 0.7, 0.8))
ucsWaterfall <- annotate_figure(ucsWaterfall, bottom = text_grob("Uterine Samples", size = 16))
ucsWaterfall
# dev.off()

# Creating a table showing the correlation values of the marker expression vs UCS
wUCSCorr <- rbind(wUCSGsva$Corr_DF, wUCSSingscore$Corr_DF, wUCSSsgsea$Corr_DF, wUCSSspaths$Corr_DF)
wUCSCorr$Phenotype <- ifelse(wUCSCorr$HGNC_Symbol %in% c("VIM", "FN1", "CDH2"), "Mesenchymal", "Epithelial")
wUCSCorr <- wUCSCorr[, c("HGNC_Symbol", "Corr_Coeff", "pValues", "pAdj", "Stars", "Phenotype", "Method")]
wUCSCorr <- wUCSCorr[order(wUCSCorr[, "Method"], wUCSCorr[, "Phenotype"]), ]
colnames(wUCSCorr)[5] <- "Signif_Symbol"
# pdf("Final Figures/uterine_correlation_table_only_UCS.pdf", width = 12, height = 16)
wUCSCorrTable <- grid.table(wUCSCorr, rows = NULL)

# Creating an .xlsx file for the UCS marker expression and score correlation
# write.xlsx(wUCSCorr, "Supplementary Tables/uterine_correlation_table_only_UCS.xlsx")
# dev.off()



####################
# External Datasets#
####################

# HNSC Single Cell Data
gsvaHNSC <- read.table("gsva_hnscSC_emt_scores.tsv", sep = "\t", stringsAsFactors = FALSE,
                       header = TRUE)
singscoreHNSC <- read.table("singscore_hnscSC_emt_scores.tsv", sep = "\t", stringsAsFactors = FALSE,
                       header = TRUE)
ssgseaHNSC <- read.table("ssgsea_hnscSC_emt_scores.tsv", sep = "\t", stringsAsFactors = FALSE,
                       header = TRUE)
sspathsHNSC <- read.table("sspaths_hnscSC_emt_scores_runABS.tsv", sep = "\t", stringsAsFactors = FALSE,
                       header = TRUE)
# Applying a function that formats the dataframes for plotting. Details in emtPaperFunctions.R
gsvaHNSC <- hnscEData(gsvaHNSC, "GSVA")
singscoreHNSC <- hnscEData(singscoreHNSC, "singscore")
ssgseaHNSC <- hnscEData(ssgseaHNSC, "ssGSEA")
sspathsHNSC <- hnscEData(sspathsHNSC, "ssPATHS")

# Flipping the direction of singscore to match the others
singscoreHNSC$emt <- -singscoreHNSC$emt
singscoreHNSC$Scaled_emt <- -singscoreHNSC$Scaled_emt
# Combining all dataframes to plot
allHNSC <- rbind(gsvaHNSC, singscoreHNSC, ssgseaHNSC, sspathsHNSC)

# Creating the desired dark colors
# Colorblind palette
boxplotColors = brightness(c("#D55E00", "#0072B2"), 0.4)

# Signficance DF
# Applying FDR correction and setting up for plotting p-values
signifHNSCAll <- unique(allHNSC[, c("Method", "pValue")])
signifHNSCAll$Stars <- signifStars(signifHNSCAll$pValue)
signifHNSCAll$Start <- c(0.75, 1.75, 2.75, 3.75)
signifHNSCAll$End <- signifHNSCAll$Start + 0.5
signifHNSCAll$Group <- 1:nrow(signifHNSCAll)

# pdf("Final Figures/grouped_HNSC_scaled.pdf", width = 7, height = 5)
ggplot(allHNSC, aes(x = Method, y = Scaled_emt, color = Condition, fill = Condition)) + geom_boxplot() +
  theme_bw() +  scale_fill_manual(values = boxplotColors) + scale_color_manual(values = c("#D55E00", "#0072B2")) +
  ylab("Score") + labs(color = "Group", fill = "Group") +
    geom_signif(data = signifHNSCAll, aes(xmin = Start, xmax = End, annotations = Stars, group = Group,
                                       y_position = c(3.4, 2.3, 4.4, 3.4)),
                color = "black" , manual = TRUE, inherit.aes = FALSE)
# dev.off()


# Finding Medians for each method
methods = as.character(unique(allHNSC$Method))
conditions = as.character(unique(allHNSC$Condition))
dfHNSCmedians = matrix(NA, ncol = 4, nrow = 3)
for(j in 1:4){
  for(i in 1:2){
    dfHNSCmedians[i, j] = median(allHNSC[allHNSC$Method == methods[j] & allHNSC$Condition == conditions[i], "Scaled_emt"])
  }
}

dfHNSCmedians[3, ] = dfHNSCmedians[1, ] - dfHNSCmedians[2, ]
colnames(dfHNSCmedians) = methods
dfHNSCmedians = data.frame(Condition = c(conditions, "Difference"), dfHNSCmedians)

# pdf("Final Figures/hnsc_median_differences.pdf", width = 10, height = 4)
grid.table(dfHNSCmedians, rows = NULL)
# dev.off()

# Creating a .tsv file to recreate the above boxplots
# HNSCplotDF <- merge(allHNSC, signifHNSCAll[, -2], by = "Method")
# write.table(HNSCplotDF[, -c(5, 6, 12)], "Plot Creation Files/hnsc_boxplots.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)

# Comparing malignant and nonMalignant for all 
signifAllSamps <- compare_means(allHNSC, formula = emt ~ Condition, 
                                group.by =  c("SampleNo", "Method"), p.adjust.method = "fdr")
hnscSignifTable <- table(signifAllSamps[, c("p.signif", "Method")])
hnscSignifTable <- cbind(row.names(hnscSignifTable), hnscSignifTable)
hnscSignifTable <- as.data.frame(hnscSignifTable)
colnames(hnscSignifTable)[1] <- "Signif_Symbol"
hnscSignifTable$Signif_Symbol <- factor(hnscSignifTable$Signif_Symbol, levels = c("ns", "*", "**", "***", "****"))
hnscSignifTable <- hnscSignifTable[order(hnscSignifTable$Signif_Symbol), ]

# pdf("Final Figures/hnsc_sample_source_MvsNonM_signif_table.pdf", width = 10, height = 10)
grid.table(hnscSignifTable, rows = NULL)

# write.xlsx(hnscSignifTable, "Supplementary Tables/hnsc_sample_source_MvsNonM_signif_table.xlsx")
# dev.off()


# Comparing individual samples
# Randomly selecting 7 of the 30 unique SampleNos
# Removing samples that do not have both malignant and NonMalignant
# uniqueSamps <- unique(allHNSC[!allHNSC$SampleNo %in% c("HNSCCLN", "HNSCC8P", "HN23P"), "SampleNo"])
# randSamp <- sample(uniqueSamps, 7)
# Vector of seven samples I randomly went with is below
randSamp <- c("HNSCC12P", "HNSCC25LN", "HNSCC18P", "HNSCC5LN", "HNSCC24P", "HNSCC6P", "HN28LN")
selectedSampDF <- allHNSC[allHNSC$SampleNo %in% randSamp, ]

boxplotColors = brightness(c("#D55E00", "#0072B2"), 0.4)
# Creating a dataframe for manual plotting of significance with geom_signif
signifHNSC <- compare_means(data = selectedSampDF, formula = emt ~ Condition, 
                               group.by = c("SampleNo", "Method"), p.adjust.method = "fdr")
signifHNSC <- data.frame(signifHNSC)
signifHNSC <- signifHNSC[order(signifHNSC$Method, signifHNSC$SampleNo), ]
signifHNSC$Stars <- signifStars(signifHNSC$p.adj)
signifHNSC$Start <- c(rep(0.75, 7), rep(1.75, 7), rep(2.75, 7), rep(3.75, 7))
signifHNSC$End <- signifHNSC$Start + 0.5
signifHNSC$Group <- 1:nrow(signifHNSC)

# pdf("Final Figures/by_sample_HNSC.pdf", width = 10, height = 6)
ggplot(selectedSampDF, aes(x = Method, y = Scaled_emt, color = Condition, fill = Condition)) +
  geom_jitter(aes(color = Condition, alpha = 0.2), position = position_jitterdodge()) + geom_boxplot(alpha = 0.7) +
  theme_bw() + scale_color_manual(values = boxplotColors) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) + labs(color = "Group", fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Score") + facet_grid(~SampleNo, scales = "free", space = "free") + guides(alpha = FALSE) +
  geom_signif(data = signifHNSC, aes(xmin = Start, xmax = End, annotations = Stars, group = Group,
                               y_position = c(rep(3.2, 7), rep(2.2, 7), rep(4.2, 7), rep(3.2, 7))),
              color = "black" , manual = TRUE, inherit.aes = FALSE)
# dev.off()

#Boxplots df and Boxplots by sample signifcance df
# write.table(selectedSampDF[, -c(5, 8)], "Plot Creation Files/hnsc_by_sample_boxplots.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# write.table(signifHNSC[, -c(3, 8, 9, 14)], "Plot Creation Files/hnsc_by_sample_boxplots_signif.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


 
# Timeline Data
# Timeline expression data
timelineExpr <- read.table("GSE125365_Timeline_new.tsv", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
# Selecting for 267 genes (includes lowly expressed genes)
geneIDs <- read.table("All_EMT_HGNC_Ensembl ID.csv", header = TRUE, stringsAsFactors = FALSE,
                      sep = ",")
emtGenes <- colnames(timelineExpr)[colnames(timelineExpr) %in% geneIDs$Ensembl_ID]
timelineExpr <- timelineExpr[, c("sample_id", emtGenes)]
# Scaling by sample across all availible EMT genes (263)
timelineExpr[, 2:ncol(timelineExpr)] <- t(scale(t(timelineExpr[, 2:ncol(timelineExpr)])))

# Timeline scores
timelineGsva <- read.table("gsva_timeSeries_emt_scores.tsv", header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
timelineSingscore <- read.table("singscore_timeSeries_emt_scores.tsv", header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
timelineSsgsea <- read.table("ssgsea_timeSeries_emt_scores.tsv", header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
timelineSspaths <- read.table("sspaths_timeSeries_emt_scores_runABS.tsv", header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)

# Boxplots for visualization
timelineGsvaList<- timelineEData(timelineGsva, "GSVA")
timelineSingscoreList <- timelineEData(timelineSingscore, "singscore")
timelineSsgseaList <- timelineEData(timelineSsgsea, "ssGSEA")
timelineSspathsList <- timelineEData(timelineSspaths, "ssPATHS")

ggarrange(timelineGsvaList$Boxplot, timelineSingscoreList$Boxplot,
          timelineSsgseaList$Boxplot, timelineSspathsList$Boxplot, nrow = 4)

# Plotting the expression with LOESS of markers 
# Selecting for important genes
impGenes <- c("FN1", "CDH1", "CDH2", "VIM", "OCLN", "DSP", "SNAI1", "SNAI2", "ZEB1", "ZEB2")
impGenes <- geneIDs[geneIDs$HGNC_Symbol %in% impGenes, ]
tlExpr <- timelineExpr[, c("sample_id", impGenes$Ensembl_ID)]

# Creating a dataframe of mesenchymal and epithelial markers grouped
meltExpr <- tlExpr
meltExpr <- melt(meltExpr, id.vars = "sample_id", variable.name = "Ensembl_ID", value.name = "Expr")
meltExpr <- merge(meltExpr, geneIDs, by = "Ensembl_ID")
meltExpr$Time <- as.numeric(gsub(pattern = "^[^_]*_([^.]*).*", replacement = "\\1", x = meltExpr$sample_id)) 
mesMark <- c("FN1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2")
meltExpr$Phenotype <- ifelse(meltExpr$HGNC_Symbol %in% mesMark, "Grouped Mes.\nMarkers", "Grouped Epi.\nMarkers")

# Combining the data
allTLdfs <- rbind(timelineGsvaList$Data, timelineSingscoreList$Data, timelineSsgseaList$Data,
                  timelineSspathsList$Data)
# Scaling the scores by method so they are more comparable on the plot
allTLdfs <- dcast(allTLdfs[, c(1, 2, 5)], sample_id ~ Method, value.var = "emt")
allTLdfs[, 2:ncol(allTLdfs)] <- scale(allTLdfs[, -1])
notMeltedAllTLdfs <- allTLdfs
allTLdfs <- melt(allTLdfs, id.vars = "sample_id", variable.name = "Method", value.name = "emt")
allTLdfs$NumTime <- gsub(pattern = "^[^_]*_([^.]*).*", replacement = "\\1", x = allTLdfs$sample_id)
allTLdfs$NumTime <- as.numeric(allTLdfs$NumTime)

# Plotting timeline data
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Plotting the collected Mes and Epi markers
# As the samples were on the log10 scale they needed to be tweaked to remove zeros
allTLdfs$NumTimeNew <- allTLdfs$NumTime + 1
meltExpr$TimeNew <- meltExpr$Time + 1
# pdf("Final Figures/timeline_regression.pdf", width = 6, height = 4)
timelinePlot <- ggplot(allTLdfs, aes(x = NumTimeNew, y = emt, color = Method)) + 
  geom_smooth(se = FALSE, lwd = 1, method = "loess") +
  theme_bw() + geom_point() + scale_color_manual(values = cbPalette) + 
  labs(x = "Time (Hours)", y = "Score", color = "Method") + 
  geom_smooth(data = meltExpr, aes(x = TimeNew, y = Expr, color = Phenotype), 
              lwd = 1.5, se = FALSE, linetype = "twodash", method = "loess") +
  scale_y_continuous(sec.axis = dup_axis(name = "Expression")) + scale_x_log10()
timelinePlot 
# dev.off()

# Plot creation files
# write.table(meltExpr, "Plot Creation Files/timeline_marker_expression.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)
# 
# write.table(allTLdfs, "Plot Creation Files/timeline_method_scores.tsv", sep = "\t",
#             row.names = FALSE, col.names = TRUE)


# Combining expression and score DFs
tlExpr <- tlExpr[, c("sample_id", impGenes$Ensembl_ID)]
colnames(tlExpr) <- c("sample_id", impGenes$HGNC_Symbol)
finalCorrTLdf <- merge(tlExpr, notMeltedAllTLdfs, by = "sample_id")
epiMark <- c("CDH1", "OCLN", "DSP")
mesMark <- c("FN1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2")
# Calculating median expression values for epithelial and mesenchymal markers
finalCorrTLdf$GroupEpi <- apply(finalCorrTLdf[, epiMark], 1, median)
finalCorrTLdf$GroupMes <- apply(finalCorrTLdf[, mesMark], 1, median)

# Calculating the correlation of individual and grouped marker expression and score
tlCorrIndGenes <- tlCorrScaled(df = finalCorrTLdf, grouped = FALSE)
tlCorrGrouped <- tlCorrScaled(df = finalCorrTLdf, grouped = TRUE)

# pdf("Final Figures/timeline_correlation_table_individual_genes.pdf", width = 10, height = 14)
tlCorrIndGenes <- tlCorrIndGenes[order(tlCorrIndGenes[, "Method"], tlCorrIndGenes[, "Phenotype"]), ]
tlCorrIndGenesPdf <- grid.table(tlCorrIndGenes, rows = NULL)
# write.xlsx(tlCorrIndGenes, "Supplementary Tables/timeline_correlation_table_individual_genes.xlsx")
# dev.off()

# pdf("Final Figures/timeline_correlation_table_grouped_genes.pdf", width = 10, height = 10)
tlCorrGroupedPdf <- grid.table(tlCorrGrouped, rows = NULL)
# write.xlsx(tlCorrGrouped, "Supplementary Tables/timeline_correlation_table_grouped_genes.xlsx")
# dev.off()



