#' @importFrom stats coef cor fisher.test t.test
#' @importFrom utils install.packages
#' @title extract gene ratios with non-zero coefficient
#' @description
#' this function is used to extract important gene ratios from dataset for each subtype based on genes ratios selected from RRA() function
#'
#' @param data the original dataset we extract features
#' @param sampAnnote sampAnnot is the annotation file for each samples in data
#' @param all_rank_t_genes all_rank_t_genes is the differentially ranked genes selected in previous RCA() function
#' @param all_reversed_gp_genes all_reversed_gp_genes is the differentially ranked gene ratios selected in previous RRA() function
#'
#' @return a global object MB_RANK_GP
#' @export
#'
#' @examples MB_RANK_GP<-LaSelect(GSE85217, sampAnnote_GSE85217, all_rank_t_genes,all_reversed_gp_genes)
LaSelect<-function(data,sampAnnote, all_rank_t_genes, all_reversed_gp_genes){
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    install.packages("glmnet")
  }
  library("glmnet")
  set.seed(13)
  createRatio <- function(exprs, x) {
    # Extract the indices or names for g1 and g2
    g1 <- as.matrix(x[1])
    g2 <- as.matrix(x[2])
    # Calculate the difference
    difference <- exprs[g1,] - exprs[g2,]
    # Apply the condition: if difference > 0 then 1 else 0
    g1g2_ratio <- ifelse(difference > 0, 1, 0)
    # Return the result
    return(g1g2_ratio)
  }
  GSE85217<-data #the data format should be: each row is a sample, each column is a gene
  sampAnnot_GSE85217<-sampAnnote

  ranked_85217 <- apply(GSE85217, 2, function(x) rank(x, ties.method = "average"))
  SHH<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="SHH",1]
  WNT<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="WNT",1]
  Group3<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="Group3",1]
  Group4<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="Group4",1]

  SHH_df<-ranked_85217[,colnames(ranked_85217) %in% SHH]
  WNT_df<-ranked_85217[,colnames(ranked_85217) %in% WNT]
  Group3_df<-ranked_85217[,colnames(ranked_85217) %in% Group3]
  Group4_df<-ranked_85217[,colnames(ranked_85217) %in% Group4]


  SHH_rank_t_gene<-all_rank_t_genes$SHH_ranked_gene
  WNT_rank_t_gene<-all_rank_t_genes$WNT_ranked_gene
  Group3_rank_t_gene<-all_rank_t_genes$Group3_ranked_gene
  Group4_rank_t_gene<-all_rank_t_genes$Group4_ranked_gene

  SHH_vs_all_reversed_gp_gene<- all_reversed_gp_genes$SHH_ranked_gene
  WNT_vs_all_reversed_gp_gene<- all_reversed_gp_genes$WNT_ranked_gene
  Group3_vs_all_reversed_gp_gene<- all_reversed_gp_genes$Group3_ranked_gene
  Group4_vs_all_reversed_gp_gene<- all_reversed_gp_genes$Group4_ranked_gene

  corGenes <- cor(t(ranked_85217[SHH_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))

  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_85217), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_85217))

  SHH_generatioout<-geneRatioOut
  colnames(SHH_generatioout) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(SHH_generatioout) <- colnames(ranked_85217)


  #set.seed(43)
  SHH_gps<-SHH_vs_all_reversed_gp_gene
  SHH_generatioout<-SHH_generatioout[,colnames(SHH_generatioout) %in% SHH_gps]
  SHH_generatioout_SHH<-data.frame(SHH_generatioout[rownames(SHH_generatioout) %in% SHH, ])
  SHH_generatioout_SHH$group<-1
  SHH_generatioout_WNT<-data.frame(SHH_generatioout[rownames(SHH_generatioout) %in% WNT, ])
  SHH_generatioout_WNT$group<-0
  SHH_generatioout_Group3<-data.frame(SHH_generatioout[rownames(SHH_generatioout) %in% Group3, ])
  SHH_generatioout_Group3$group<-0
  SHH_generatioout_Group4<-data.frame(SHH_generatioout[rownames(SHH_generatioout) %in% Group4, ])
  SHH_generatioout_Group4$group<-0
  #####################################
  #####################################
  #####################################d1_SHH_gene_selection

  # SHH_VS_WNT, column is gene pair, row is sample
  SHH_VS_WNT<-rbind(SHH_generatioout_SHH,SHH_generatioout_WNT)

  X <- as.matrix(SHH_VS_WNT[, 1:(length(colnames(SHH_VS_WNT))-1)])  # Predictors (binary matrix)
  Y <- SHH_VS_WNT[, (length(colnames(SHH_VS_WNT)))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  SHH_WNT_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]


  # SHH_VS_Group3, column is gene pair, row is sample
  SHH_VS_Group3<-rbind(SHH_generatioout_SHH,SHH_generatioout_Group3)

  X <- as.matrix(SHH_VS_Group3[, 1:(length(colnames(SHH_VS_Group3))-1)])  # Predictors (binary matrix)
  Y <- SHH_VS_Group3[, length(colnames(SHH_VS_Group3))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  SHH_Group3_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]
  #275
  # SHH_VS_Group4, column is gene pair, row is sample
  SHH_VS_Group4<-rbind(SHH_generatioout_SHH,SHH_generatioout_Group4)
  X <- as.matrix(SHH_VS_Group4[, 1:(length(colnames(SHH_VS_Group4))-1)])  # Predictors (binary matrix)
  Y <- SHH_VS_Group4[, length(colnames(SHH_VS_Group4))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  SHH_Group4_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]
  #292
  SHH_mb_rank_gp<- Reduce(intersect, list(SHH_Group3_selected_genes_ridge$gene,SHH_WNT_selected_genes_ridge$gene,SHH_Group4_selected_genes_ridge$gene))
  SHH_mb_rank_gp <- SHH_mb_rank_gp[SHH_mb_rank_gp != "(Intercept)"]

  #211
  ##########################################################for wnt_mb_rank_gp
  corGenes <- cor(t(ranked_85217[WNT_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))


  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_85217), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_85217))

  WNT_generatioout<-geneRatioOut
  colnames(WNT_generatioout) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(WNT_generatioout) <- colnames(ranked_85217)

  WNT_gps<-WNT_vs_all_reversed_gp_gene
  WNT_generatioout<-WNT_generatioout[,colnames(WNT_generatioout) %in% WNT_gps]

  WNT_generatioout_WNT<-data.frame(WNT_generatioout[rownames(WNT_generatioout) %in% WNT, ])
  WNT_generatioout_WNT$group<-1
  WNT_generatioout_SHH<-data.frame(WNT_generatioout[rownames(WNT_generatioout) %in% SHH, ])
  WNT_generatioout_SHH$group<-0
  WNT_generatioout_Group3<-data.frame(WNT_generatioout[rownames(WNT_generatioout) %in% Group3, ])
  WNT_generatioout_Group3$group<-0
  WNT_generatioout_Group4<-data.frame(WNT_generatioout[rownames(WNT_generatioout) %in% Group4, ])
  WNT_generatioout_Group4$group<-0

  # Group4_vs_SHH, column is gene pair, row is sample
  WNT_vs_SHH<-rbind(WNT_generatioout_WNT,WNT_generatioout_SHH)

  X <- as.matrix(WNT_vs_SHH[, 1:(length(colnames(WNT_vs_SHH))-1)])  # Predictors (binary matrix)
  Y <- WNT_vs_SHH[, length(colnames(WNT_vs_SHH))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  WNT_SHH_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]


  # Group4_vs_WNT, column is gene pair, row is sample
  WNT_Group4<-rbind(WNT_generatioout_WNT,WNT_generatioout_Group4)

  X <- as.matrix(WNT_Group4[, 1:(length(colnames(WNT_Group4))-1)])  # Predictors (binary matrix)
  Y <- WNT_Group4[, length(colnames(WNT_Group4))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  WNT_Group4_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]


  # WNT_VS_Group3, column is gene pair, row is sample
  WNT_VS_Group3<-rbind(WNT_generatioout_WNT,WNT_generatioout_Group3)
  X <- as.matrix(WNT_VS_Group3[, 1:(length(colnames(WNT_VS_Group3))-1)])  # Predictors (binary matrix)
  Y <- WNT_VS_Group3[, length(colnames(WNT_VS_Group3))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  WNT_Group3_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]

  WNT_mb_rank_gp<- Reduce(intersect, list(WNT_Group3_selected_genes_ridge$gene,WNT_Group4_selected_genes_ridge$gene,WNT_SHH_selected_genes_ridge$gene))
  WNT_mb_rank_gp<-WNT_mb_rank_gp[WNT_mb_rank_gp != "(Intercept)"]


  ##########################################################FOR Group3_mb+rank_gp

  corGenes <- cor(t(ranked_85217[Group3_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))

  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_85217), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_85217))

  Group3_geneRatioOut<-geneRatioOut
  Group3_gps<-Group3_vs_all_reversed_gp_gene
  Group3_geneRatioOut<-Group3_geneRatioOut[,colnames(Group3_geneRatioOut) %in% Group3_gps]

  Group3_geneRatioOut_Group4<-data.frame(Group3_geneRatioOut[rownames(Group3_geneRatioOut) %in% Group4, ])
  Group3_geneRatioOut_Group4$group<-0
  Group3_geneRatioOut_SHH<-data.frame(Group3_geneRatioOut[rownames(Group3_geneRatioOut) %in% SHH, ])
  Group3_geneRatioOut_SHH$group<-0
  Group3_geneRatioOut_Group3<-data.frame(Group3_geneRatioOut[rownames(Group3_geneRatioOut) %in% Group3, ])
  Group3_geneRatioOut_Group3$group<-1
  Group3_geneRatioOut_WNT<-data.frame(Group3_geneRatioOut[rownames(Group3_geneRatioOut) %in% WNT, ])
  Group3_geneRatioOut_WNT$group<-0

  # Group4_vs_SHH, column is gene pair, row is sample
  Group3_vs_SHH<-rbind(Group3_geneRatioOut_Group3,Group3_geneRatioOut_SHH)

  X <- as.matrix(Group3_vs_SHH[, 1:(length(colnames(Group3_vs_SHH))-1)])  # Predictors (binary matrix)
  Y <- Group3_vs_SHH[, length(colnames(Group3_vs_SHH))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group3_SHH_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]

  # Group4_vs_Group3, column is gene pair, row is sample
  Group3_Group4<-rbind(Group3_geneRatioOut_Group3,Group3_geneRatioOut_Group4)

  X <- as.matrix(Group3_Group4[, 1:(length(colnames(Group3_Group4))-1)])  # Predictors (binary matrix)
  Y <- Group3_Group4[, length(colnames(Group3_Group4))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group3_Group4_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]



  # Group3_VS_WNT, column is gene pair, row is sample
  Group3_VS_WNT<-rbind(Group3_geneRatioOut_Group3,Group3_geneRatioOut_WNT)
  X <- as.matrix(Group3_VS_WNT[, 1:(length(colnames(Group3_VS_WNT))-1)])  # Predictors (binary matrix)
  Y <- Group3_VS_WNT[, length(colnames(Group3_VS_WNT))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group3_WNT_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]
  #32
  Group3_mb_rank_gp<- Reduce(intersect, list(Group3_Group4_selected_genes_ridge$gene,Group3_WNT_selected_genes_ridge$gene,Group3_SHH_selected_genes_ridge$gene))
  Group3_mb_rank_gp<-Group3_mb_rank_gp[Group3_mb_rank_gp !="(Intercept)"]


  ###################################################for Group4_mb_rank_gp


  corGenes <- cor(t(ranked_85217[Group4_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))


  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_85217), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_85217))

  Group4_generatioout<-geneRatioOut
  colnames(Group4_generatioout) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(Group4_generatioout) <- colnames(ranked_85217)

  Group4_gps<-Group4_vs_all_reversed_gp_gene
  Group4_generatioout<-Group4_generatioout[,colnames(Group4_generatioout) %in% Group4_gps]
  Group4_generatioout_SHH<-data.frame(Group4_generatioout[rownames(Group4_generatioout) %in% SHH, ])
  Group4_generatioout_SHH$group<-0
  Group4_generatioout_WNT<-data.frame(Group4_generatioout[rownames(Group4_generatioout) %in% WNT, ])
  Group4_generatioout_WNT$group<-0
  Group4_generatioout_Group3<-data.frame(Group4_generatioout[rownames(Group4_generatioout) %in% Group3, ])
  Group4_generatioout_Group3$group<-0
  Group4_generatioout_Group4<-data.frame(Group4_generatioout[rownames(Group4_generatioout) %in% Group4, ])
  Group4_generatioout_Group4$group<-1
  #####################################
  #####################################
  #####################################d1_SHH_gene_selection

  # SHH_VS_WNT, column is gene pair, row is sample
  Group4_VS_WNT<-rbind(Group4_generatioout_Group4,Group4_generatioout_WNT)

  X <- as.matrix(Group4_VS_WNT[, 1:(length(colnames(Group4_VS_WNT))-1)])  # Predictors (binary matrix)
  Y <- Group4_VS_WNT[, length(colnames(Group4_VS_WNT))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group4_WNT_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]


  # SHH_VS_Group3, column is gene pair, row is sample
  Group4_VS_Group3<-rbind(Group4_generatioout_Group4,Group4_generatioout_Group3)

  X <- as.matrix(Group4_VS_Group3[,  1:(length(colnames(Group4_VS_Group3))-1)])  # Predictors (binary matrix)
  Y <- Group4_VS_Group3[, length(colnames(Group4_VS_Group3))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group4_Group3_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]


  # Group4_vs_SHH, column is gene pair, row is sample
  Group4_vs_SHH<-rbind(Group4_generatioout_Group4,Group4_generatioout_SHH)
  X <- as.matrix(Group4_vs_SHH[, 1:(length(colnames(Group4_vs_SHH))-1)])  # Predictors (binary matrix)
  Y <- Group4_vs_SHH[, length(colnames(Group4_vs_SHH))]               # Response (binary outcome)
  # Fit model using cross-validation with an alpha value closer to ridge regression
  #set.seed(123)
  cvmod_ridge <- cv.glmnet(X, Y, family = "binomial", alpha = 0.1)  # Using alpha = 0.1 as an example
  plot(cvmod_ridge)
  # Extract lambda that gives one standard error away from the minimum (more regularized)
  lambda_1se_ridge <- cvmod_ridge$lambda.1se
  # Extract coefficients at lambda.1se
  coef_1se_ridge <- coef(cvmod_ridge, s = "lambda.1se")
  # Convert coefficients to a data frame for easier manipulation and viewing
  coef_1se_ridge_matrix <- as.matrix(coef_1se_ridge)
  coef_1se_ridge_df <- data.frame(coef_1se_ridge_matrix)
  # Add gene names or row names as a new column in the data frame
  coef_1se_ridge_df$gene <- rownames(coef_1se_ridge_df)
  # Rename the first column to 'Coefficient'
  names(coef_1se_ridge_df)[1] <- "Coefficient"
  # Filter to keep only non-zero coefficients (selected genes)
  Group4_SHH_selected_genes_ridge <- coef_1se_ridge_df[coef_1se_ridge_df$Coefficient != 0, ]

  Group4_mb_rank_gp<- Reduce(intersect, list(Group4_SHH_selected_genes_ridge$gene,Group4_Group3_selected_genes_ridge$gene,Group4_WNT_selected_genes_ridge$gene))
  Group4_mb_rank_gp<-Group4_mb_rank_gp[Group4_mb_rank_gp != "(Intercept)"]

  MB_RANK_GP <- list(
    SHH = SHH_mb_rank_gp,
    WNT = WNT_mb_rank_gp,
    Group3 = Group3_mb_rank_gp,
    Group4 = Group4_mb_rank_gp
  )

  return(MB_RANK_GP)

}
