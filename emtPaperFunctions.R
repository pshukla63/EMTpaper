###########################
# Functions for EMT Paper #
###########################

# Function to find correlation of miRNA expression with method score using linear regression function lm(). 
# Inputs:
  # methodDF: data frame with columns "sample_id" and "emt" (score)
  # mirnaDF: data frame with columns "sample_id", "Study_Abbr", and expression by gene as the rest of the 
    # columns with Ensembl_ID column names
  # cancerType: Boolean. If TRUE calcualtes lm by cancer type listed 
    # i.e. lm(method score for cancer type ~ indivdual miRNA expession for that cancer type)
    # If FALSE ignores cancer type. Calculates lm(entire method score ~ corresponding individual miRNA expression)
  # method: string indicating the method to calculate correlation for e.g "GSVA"
# Returns: dataframe with columns "miRNA", "Method", "Study_Abbr", "pValue", "Estimate", "lower", "upper"
  # where "lower" and "upper" are the 0.95 confidence intervals 
lmScoreByMirna <- function(methodDF, mirnaDF, cancerType, method){
  
  # Trimming the sample_ids so they match in both dfs (mirnaDF has truncated IDs)
  methodDF$sample_id = strtrim(methodDF$sample_id, 16)
  mergedDFs = merge(methodDF, mirnaDF, by = "sample_id")
  # mergedDFs = mergedDFs[!duplicated(mergedDFs$sample_id), ]
  
  if(cancerType == TRUE){
    # Selecting for only eight cancers to make visualization easy
    impCanc = c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "STAD", "UCEC")
    mergedDFs = mergedDFs[mergedDFs$Study_Abbr %in% impCanc,]
    # Changing Study_Abbr variables to factors
    mergedDFs$Study_Abbr <- factor(
      mergedDFs$Study_Abbr,
      levels = c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", "STAD", "UCEC")
    )
    
    # Creating an empty vector to store all p-values
    pval = matrix(NA, nrow = 9, ncol = 8)
    # Creating an empty vector to store all coefficients estimates
    betas = matrix(NA, nrow = 9, ncol = 8)
    # Getting the standard deviation to calculate confidence intervals
    confInt = matrix(NA, nrow = 9, ncol = 16)
    # Iterating through the df to select for miRNA expression columns
    for (i in 4:ncol(mergedDFs)) {
      for (j in 1:8) {
        # Selecting for the cancer types one by one
        tempdf = mergedDFs[mergedDFs$Study_Abbr == impCanc[j],]
        # Calculating lm for expression for a single cancer type 
        tempLM = lm(emt ~ tempdf[, i], data = tempdf)
        # Extracting p-value, coefficent estimates, and confidence intervals
        pval[i - 3, j] = summary(tempLM)$coefficients[2, 4]
        betas[i - 3, j] = summary(tempLM)$coefficients[2, 1]
        confInt[i - 3, j] = confint(tempLM)[2, 1]
        confInt[i - 3, j + 8] = confint(tempLM)[2, 2]
      }
    }
    
    # Creating the dfs to retrun
    pvaldf = data.frame(miRNA = colnames(mergedDFs)[-c(1:3)], Method = method, pval)
    colnames(pvaldf) = c("miRNA", "Method", impCanc)
    betasdf = data.frame(miRNA = colnames(mergedDFs)[-c(1:3)], Method = method, betas)
    colnames(betasdf) = c("miRNA", "Method", impCanc)
    confLowerDf = data.frame(miRNA = colnames(mergedDFs)[-c(1:3)], Method = method, confInt[, 1:8])
    colnames(confLowerDf) = c("miRNA", "Method", impCanc)
    confUpperDf = data.frame(miRNA = colnames(mergedDFs)[-c(1:3)], Method = method, confInt[, 9:16])
    colnames(confUpperDf) = c("miRNA", "Method", impCanc)
    
    # Melting the above dfs to create columns and combining the data
    meltpvalDF = melt(
      pvaldf,
      id.vars = c("miRNA", "Method"),
      variable.name = "Study_Abbr",
      value.name = "pValue"
    )
    meltbetaDF = melt(
      betasdf,
      id.vars = c("miRNA", "Method"),
      variable.name = "Study_Abbr",
      value.name = "Estimate"
    )
    meltconfLowerDF = melt(
      confLowerDf,
      id.vars = c("miRNA", "Method"),
      variable.name = "Study_Abbr",
      value.name = "Lower"
    )
    
    meltconfUpperDF = melt(
      confUpperDf,
      id.vars = c("miRNA", "Method"),
      variable.name = "Study_Abbr",
      value.name = "Upper"
    )
    
    combDF = cbind(
      meltpvalDF,
      Estimate = meltbetaDF$Estimate,
      lower = meltconfLowerDF$Lower,
      upper = meltconfUpperDF$Upper
    )
  }
  
  # If cancerType == FALSE lm(entire method score ~ individual miRNA expression) is calculated
  else{
    pval = rep(NA, length(4:ncol(mergedDFs)))
    betas = rep(NA, length(4:ncol(mergedDFs)))
    
    for(i in 4:ncol(mergedDFs)){
      tempLM = lm(mergedDFs[, "emt"] ~ mergedDFs[, i], data = mergedDFs)
        pval[i-3] = summary(tempLM)$coefficients[2, 4]
        betas[i-3] = summary(tempLM)$coefficients[2, 1] 
    }
    
    # Creating output DF 
    combDF = data.frame(miRNA = colnames(mergedDFs)[4:ncol(mergedDFs)], 
                        Estimate = betas, pValue = pval, Method = method)
    
  }
  
  return(combDF)
}


# Function to create a waterfall plot (or whatever those are called) and heatmap of expression of 
# epithelial and mesenchymal markers 
# Input dataframes only with samples you want run!
# Inputs:
  # dfWscore: data frame with columns "sample_id", score and Tissue_Type
  # dfWexpr: data frame with "sample_id", expression and Ensembl_ID column names - expression is scaled by 
    # sample and log10(x + 1) transformed within the function
  # ylabel: string - y-axis label of waterfall plot
  # xlabel: string - x-axis label
  # legendTitle: sting - legend title (optional - "Tissue Type" if empty)
  # scorePlotBreaks: numeric - the number of breaks for the y-axis of the score plot (optional - default 5)
  # markerYaxisSize: numeric = the font size of the y-axis labels
  # markerGradientTitle: numeric - size of maker gradient legend titles 
  # legendVisible: logical - if FALSE, the legends of all three plots are removed. 
    # Keeps default settings otherwise
# Returns: a list with the final plot, correlation df, expression df, and score df
markersBySample <- function(dfWscore, dfWexpr, ylabel, xlabel, legendTitle, scorePlotBreaks, 
                            markerYaxisSize, markerGradientTitle, legendVisible){
  
  # Scaling and transforming expression by sample
  dfWexpr[, 2:ncol(dfWexpr)] = t(scale(log10(t(dfWexpr[, 2:ncol(dfWexpr)] + 1))))
  
  markers = data.frame(HGNC_Symbol = c("CDH1", "CDH2", "DSP", "FN1", "OCLN", "VIM"), 
                       Ensembl_ID = c("ENSG00000039068", "ENSG00000170558", "ENSG00000096696", 
                                      "ENSG00000115414", "ENSG00000197822", "ENSG00000026025"))
  epiMarkers = c("CDH1", "DSP", "OCLN")
  mesMarkers = c("CDH2", "VIM", "FN1")
  
  # Removing duplicate samples since they may have different scores
  dups = dfWscore[duplicated(dfWscore$sample_id), "sample_id"]
  dfScore = dfWscore[!(dfWscore$sample_id %in% dups), ]
  dfExpr = dfWexpr[!(dfWexpr$sample_id %in% dups), 
                   c(colnames(dfWexpr)[1], as.character(markers$Ensembl_ID))]
  
  # Ordering sample_id in both dfs to make them identical
  dfScore$sample_id = factor(dfScore$sample_id, levels = dfScore[order(dfScore$emt), "sample_id"])
  dfExpr$sample_id = factor(dfExpr$sample_id, levels = dfScore[order(dfScore$emt), "sample_id"])

  
  # Manipulating the expression df to be plottable
  meltExpr = melt(dfExpr, id.vars = colnames(dfExpr)[1], variable.name = "Ensembl_ID", 
                  value.name = "Expression")
  meltExpr = merge(meltExpr, markers, by = "Ensembl_ID")
  meltExpr$Marker_Type = ifelse(meltExpr$HGNC_Symbol %in% epiMarkers, 
                                "Epithelial Marker", "Mesenchymal Marker") 
  # Ordering sample_id in both dfs to make them identical
  meltExpr$sample_id = factor(meltExpr$sample_id, levels = dfScore[order(dfScore$emt), "sample_id"])
 
  
  # Plotting waterfall plot or whatever that column plot is called and the expression plots
  # Changing the colors based on the number of tissue types ("E", "E/M", or "M") are present
  if(length(unique(dfScore$Tissue_Type)) == 1){
    blueSat = saturation(c("#0072B2"), 0.8)
    colors = brightness(c("#D55E00", "#009E59", blueSat), 0.4)
    scoreColors = brightness(colors, 0.85)
    scoreColors = scoreColors[[2]]
    high = colors[[3]]
  }
  else if(length(unique(dfScore$Tissue_Type)) == 2){
    blueSat = saturation(c("#0072B2"), 0.8)
    colors = brightness(c("#D55E00", blueSat), 0.4)
    scoreColors = brightness(colors, 0.85)
    high = colors[[2]]

  }
  else{
    blueSat = saturation(c("#0072B2"), 0.8)
    colors = brightness(c("#D55E00", "#009E59", blueSat), 0.4)
    scoreColors = brightness(colors, 0.85)
    high = colors[[3]]
  }

  # Mid color for the red in the expression plot to make the distinctions clearer
  mid = brightness("#0072B2", 1)
  # Mid color for the blue
  colors1mid = brightness("#D55E00", 1)
  
  if(missing(ylabel)){
    ylabel = "Score"
  }
  
  if(missing(legendTitle)){
    legendTitle = "Tissue Type    "
  }

  if(missing(scorePlotBreaks)){
    scorePlotBreaks = 5
  }
  
  if(missing(markerYaxisSize)){
    markerYaxisSize = 16.5
  }
  
  if(missing(markerGradientTitle)){
    markerGradientTitle = 14.88
  }
  
  # Calculating and printing correlation significance on the plot
  df = dfWexpr[, c("sample_id", "ENSG00000039068", "ENSG00000170558", "ENSG00000096696", 
              "ENSG00000115414", "ENSG00000197822", "ENSG00000026025")]
  df = merge(dfWscore, df, by = "sample_id")
    corTemp = c()
    corrVec = rep(NA, 6)
    pvalues = rep(NA, 6)
    for(i in 1:6){
      corrTemp = cor.test(df$emt, df[, (i + 4)], method = "spearman", exact = FALSE)
      corrVec[i] = corrTemp$estimate
      pvalues[i] = corrTemp$p.value
    }
    
    padj = p.adjust(pvalues, method = "fdr")
    wstars = signifStars(padj)
    corrDF = data.frame(Genes = colnames(df)[-c(1:4)], Corr_Coeff = corrVec, pValues = pvalues, 
                        pAdj = padj, Stars = wstars, Method = ylabel)
    corrDF = merge(markers, corrDF, by.x = "Ensembl_ID", by.y = "Genes")
  
  # Creating a default y-axis based on the max and min scores for the score plot
  yaxisBreaks = round(seq(floor(min(dfScore$emt)), ceiling(max(dfScore$emt)), length.out = scorePlotBreaks), 2)
  
  scorePlot = ggplot(dfScore, aes(x = sample_id, y = emt, color = Study_Abbr, fill = Study_Abbr)) + 
    geom_col()+ scale_y_continuous(name = ylabel, breaks = yaxisBreaks) + theme_classic() + scale_fill_manual(values = scoreColors) + 
    scale_color_manual(values = scoreColors) +
    theme(axis.text.y = element_text(size = 12), axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 15),  
          text = element_text(size = 20)) +  guides(color = FALSE) + 
    labs(fill = legendTitle) 
  
  # Creating new y-axis labels for the marker plots 
  meltExpr = merge(meltExpr, corrDF[, c("Ensembl_ID", "Stars")], by = "Ensembl_ID")
  meltExpr$yLabs = paste(meltExpr$HGNC_Symbol, meltExpr$Stars, sep = "\n")
  # Epithelial Marker plot
  eMarkerPlot = ggplot(meltExpr[meltExpr$Marker_Type == "Epithelial Marker", ], 
                       aes(x = sample_id, y = yLabs, color = Expression, fill = Expression)) + geom_tile()  + theme_classic() +
    theme(axis.text.x = element_blank(), axis.title = element_blank(), line = element_blank(),
          text = element_text(size = 16.41), legend.title = element_text(size = markerGradientTitle),
          axis.text.y = element_text(size = markerYaxisSize)) + 
    labs(color = "Epi.\nMarker\nExpression", fill = "Epi.\nMarker\nExpression")  + 
    scale_color_gradientn(colors = c("white", colors1mid, colors[[1]])) + 
    scale_fill_gradientn(colors = c("white", colors1mid, colors[[1]]))
  
  # Mesenchymal Marker plot
  mMarkerPlot = ggplot(meltExpr[meltExpr$Marker_Type == "Mesenchymal Marker", ], 
                       aes(x = sample_id, y = yLabs, color = Expression, fill = Expression)) + geom_tile() + theme_classic() +
    theme(axis.text.x = element_blank(), axis.title.y = element_blank(), line = element_blank(),
          text = element_text(size = 16.41), legend.title = element_text(size = markerGradientTitle),
          axis.text.y = element_text(size = markerYaxisSize)) + 
    labs(color = "Mes.\nMarker\nExpression", fill = "Mes.\nMarker\nExpression") +
    scale_color_gradientn(colors = c("white", mid, high)) + scale_fill_gradientn(colors = c("white", mid, high))
  
  # Removes legends of all three plots 
  if(legendVisible == FALSE){
    mMarkerPlot = mMarkerPlot + theme(legend.position = "NONE")
    eMarkerPlot = eMarkerPlot + theme(legend.position = "NONE")
    scorePlot = scorePlot + theme(legend.position = "NONE")
  }

  if(missing(xlabel)){
    mMarkerPlot = mMarkerPlot + theme(axis.title.x = element_blank())
  }
  else{
    mMarkerPlot = mMarkerPlot + xlab(xlabel)
  }
  
  finalPlot = ggarrange(scorePlot, eMarkerPlot, mMarkerPlot, heights = c(1.2 ,0.6, 0.65), ncol = 1, nrow = 3)
  return(list(Plot = finalPlot, Corr_DF = corrDF, exprDF = meltExpr, scoreDF = dfScore))
}

# Function to format the HNSC Dataset to give a df with columns 
  # "sample_id", "emt", "Condition", "SampleNo", "LNorP", "Method", "Scaled_emt", "pValue"
# Input:
  # df: dataframe with columns "sample_id" and "emt" (score)
  # method: the method the scores are obtained from e.g. "GSVA"
hnscEData <- function(df, method){
  # Condition column selects for the "Malignant" or "NonMalignant" portion of the sample_id
  df$Condition = gsub(pattern = "^[^_]*_([^.]*).*", replacement = "\\1", x = df$sample_id) 
  # SampleNo column selects for the sample number e.g. "HN26LN" or "HN26P"
  df$SampleNo = gsub(pattern = "\\_.*", replacement = "", x = df$sample_id) 
  # LNorP column labels the sample origin as lymph node (LN) or primary tumor (P)
  df$LNorP = ifelse(df$SampleNo %like% "LN", "LN", "P")
  df$Method = method
  df$Scaled_emt = scale(df$emt)[, 1]
  # Function from Natalie Davidson to calculate p-value of scores of malignant vs non-malignant samples
  pval = mood_median_test(df[df$Condition == "Malignant", "emt"], df[df$Condition == "NonMalignant", "emt"])
  df$pValue = pval
  return(df)
}

# Function to plot timeline Dataset to give a list with a boxplot and a df with columns
  # "sample_id", "emt", "Time", "NumTime", "Method"
# Input:
  # df: dataframe with columns "sample_id" and "emt" (score)
  # method: the method the scores are obtained from e.g. "GSVA"
timelineEData <- function(df, method){
  # Time column selects the time points in hours as type factor
  df$Time <- gsub(pattern = "^[^_]*_([^.]*).*", replacement = "\\1", x = df$sample_id)
  df$Time <- factor(df$Time, levels = unique(df[order(as.numeric(df$Time)), "Time"]))
  # Same as Time column but with type as numeric for plotting
  df$NumTime <- as.numeric(as.character(df$Time))
  df$Method <- method
  plot1 = ggplot(df, aes(x = Time, y = emt, color = Time, fill = Time)) + geom_boxplot() + theme_bw() +
    theme(text = element_text(size = 20), legend.position = "NONE") + ylab(method) 
  g = ggplot_build(plot1)
  plotColors = brightness(unique(g$data[[1]]["fill"])[[1]], 0.4)
  finalPlot = plot1 + scale_fill_manual(values = plotColors) + 
    geom_jitter(aes(color = df$Time, alpha = 0.2), width = 0.25) + guides(alpha = FALSE) 
  return(list(Boxplot = finalPlot, Data = df))
}

# Creates a vector of stars to match the input adjsted p-values 
signifStars <- function(adjPvals){
  stars = rep(NA, length(adjPvals))
  for(i in 1:length(adjPvals)){
    if(adjPvals[i] > 0.05){
      stars[i] = "ns"
    }
    else if(adjPvals[i] <= 0.0001){
      stars[i] = "****"
    }
    else if(adjPvals[i] <= 0.001){
      stars[i] = "***"
    }
    else if(adjPvals[i] <= 0.01){
      stars[i] = "**"
    }
    else if(adjPvals[i] <= 0.05){
      stars[i] = "*"
    }
  }
  return(stars)
}

# Function from Natalie Davidson to calculate p-values using the Fisher test
mood_median_test <- function(x, y) {
  a <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(a)
  fisher.test(a < m, g)$p.value
}

# Function to calculate the correlation (using cor.test()) of timeline scaled values with score 
# Inputs: 
  # df: df with columns "sample_id", HGNC symbol genes with expression values, methods with their scores,
    # median of grouped epithelial marker expression, median of grouped mesenchymal marker expression
  # grouped: boolean. TRUE to use the median grouped expression columns in analysis. FALSE to use the 
    # individual gene expression in correlation analysis
# Returns: a df with columns "Genes", "Method", "Corr_Coeff", "pValue", "pAdj", "Signif_Symbols", "Phenotype"
tlCorrScaled <- function(df, grouped){
  markers = c("CDH1", "CDH2", "DSP", "FN1", "OCLN", 
              "SNAI1", "SNAI2", "VIM", "ZEB1", "ZEB2")
  groups = c("GroupEpi", "GroupMes")
  methods = c("GSVA", "singscore", "ssGSEA", "ssPATHS")
  
  
  if(grouped == TRUE){
    vecLength = length(groups)
    vec = groups
  }
  else{
    vecLength = length(markers)
    vec = markers
  }
  
  # Initializing p-value and coefficient estimate matrices
  pValues = matrix(NA, nrow = vecLength, ncol = 4)
  corr_Coeff = matrix(NA, nrow = vecLength, ncol = 4)
  for(i in 1:vecLength){
    for(j in 1:length(methods)){
      tempCorr = cor.test(df[, vec[i]], df[, methods[j]], method = "spearman", exact = FALSE)
      pValues[i, j] = tempCorr$p.value
      corr_Coeff[i, j] = tempCorr$estimate
    }
  }
  
  # Creating final dataframe by combining p-values and coefficients 
  pValues = cbind(vec, data.frame(pValues))
  colnames(pValues) = c("Genes", methods)
  corr_Coeff = cbind(vec, data.frame(corr_Coeff))
  colnames(corr_Coeff) = c("Genes", methods)
  
  epiMark <- c("CDH1", "OCLN", "DSP")
  
  pValues = melt(pValues, id.vars = "Genes", variable.name = "Method", value.name = "pValue")
  corr_Coeff = melt(corr_Coeff, id.vars = "Genes", variable.name = "Method", value.name = "Corr_Coeff")
  allValues = corr_Coeff
  allValues$pValue = pValues$pValue
  allValues$pAdj = p.adjust(allValues$pValue, method = "fdr")
  allValues$Signif_Symbols = signifStars(allValues$pAdj)
  allValues$Phenotype = ifelse(allValues$Genes %in% c(epiMark, "GroupEpi"), "Epithelial", "Mesenchymal")
  
  return(allValues)
}

