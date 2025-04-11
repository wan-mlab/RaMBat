#' @importFrom stats coef cor fisher.test t.test
#' @importFrom utils install.packages
#' @title extract reversed gene ratio
#' @description
#' this function is used to extract reversed gene ratios (GERs) from dataset for each subtype based on genes selected from GRA() function
#'
#' @param data the training dataset we extract gene ratios (GERs)
#' @param sampAnnot sampAnnot is the annotation file for each samples in training dataset
#' @param all_rank_t_genes all_rank_t_genes is the differentially ranked genes selected in previous GRA() function
#'
#' @return  a global object 'all_reversed_gp_genes'
#' @export
#'
#' @examples all_reversed_gp_genes<-RRA(GSE85217, sampAnnote_GSE85217, all_rank_t_genes)
RRA<-function(data, sampAnnot, all_rank_t_genes){

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


  ranked_data <- apply(data, 2, function(x) rank(x, ties.method = "average"))
  SHH<-sampAnnot[sampAnnot$Subtype=="SHH",1]
  WNT<-sampAnnot[sampAnnot$Subtype=="WNT",1]
  Group3<-sampAnnot[sampAnnot$Subtype=="Group3",1]
  Group4<-sampAnnot[sampAnnot$Subtype=="Group4",1]


  SHH_rank_t_gene<-all_rank_t_genes$SHH_ranked_gene
  WNT_rank_t_gene<-all_rank_t_genes$WNT_ranked_gene
  Group3_rank_t_gene<-all_rank_t_genes$Group3_ranked_gene
  Group4_rank_t_gene<-all_rank_t_genes$Group4_ranked_gene

  ###for SHH
  corGenes <- cor(t(ranked_data[SHH_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))


  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_data), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(ranked_data)
  SHH_reversed_gp<- geneRatioOut

  #for SHH reversed analysis
  SHH_reversed_gp<-t(SHH_reversed_gp)
  ###

  SHH_df<-SHH_reversed_gp[,colnames(SHH_reversed_gp) %in% SHH]
  WNT_df<-SHH_reversed_gp[,colnames(SHH_reversed_gp) %in% WNT]
  Group3_df<-SHH_reversed_gp[,colnames(SHH_reversed_gp) %in% Group3]
  Group4_df<-SHH_reversed_gp[,colnames(SHH_reversed_gp) %in% Group4]
  #for SHH
  SHH_ab <- data.frame(rowSums(SHH_df))
  names(SHH_ab)[1]<-"a"
  SHH_ab$b<-length(colnames(SHH_df)) - SHH_ab$a

  #for WNT
  WNT_cd <- data.frame(rowSums(WNT_df))
  names(WNT_cd)[1]<-"c"
  WNT_cd$d<-length(colnames(WNT_df)) - WNT_cd$c

  #for Group3
  Group3_cd <- data.frame(rowSums(Group3_df))
  names(Group3_cd)[1]<-"c"
  Group3_cd$d<-length(colnames(Group3_df)) - Group3_cd$c

  #for Group4
  Group4_cd <- data.frame(rowSums(Group4_df))
  names(Group4_cd)[1]<-"c"
  Group4_cd$d<-length(colnames(Group4_df)) - Group4_cd$c

  ##########calculate reversal ratio and perform fisher exact test
  ###################for SHH_WNT
  SHH_WNT<-cbind(SHH_ab,WNT_cd)
  rownames(SHH_WNT)<-rownames(SHH_reversed_gp)
  SHH_WNT$result <- with(SHH_WNT, (a + d) / (a + b + c + d))
  SHH_WNT_reversed_gp<-SHH_WNT[SHH_WNT$result>0.95,]


  df<-data.frame(SHH_WNT[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  SHH_WNT_reversed_gp_pvalue<-df[df$p_value<0.001,]
  SHH_WNT_reversed_gp_gene<-intersect(rownames(SHH_WNT_reversed_gp),rownames(SHH_WNT_reversed_gp_pvalue))

  ###################for SHH_Group3
  SHH_Group3<-cbind(SHH_ab,Group3_cd)
  rownames(SHH_Group3)<-rownames(SHH_reversed_gp)
  SHH_Group3$result <- with(SHH_Group3, (a + d) / (a + b + c + d))
  SHH_Group3_reversed_gp<-SHH_Group3[SHH_Group3$result>0.95,]


  df<-data.frame(SHH_Group3[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  SHH_Group3_reversed_gp_pvalue<-df[df$p_value<0.001,]
  SHH_Group3_reversed_gp_gene<-intersect(rownames(SHH_Group3_reversed_gp),rownames(SHH_Group3_reversed_gp_pvalue))

  ###################for SHH_Group4
  SHH_Group4<-cbind(SHH_ab,Group4_cd)
  rownames(SHH_Group4)<-rownames(SHH_reversed_gp)
  SHH_Group4$result <- with(SHH_Group4, (a + d) / (a + b + c + d))
  SHH_Group4_reversed_gp<-SHH_Group4[SHH_Group4$result>0.95,]


  df<-data.frame(SHH_Group4[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  SHH_Group4_reversed_gp_pvalue<-df[df$p_value<0.001,]

  SHH_Group4_reversed_gp_gene<-intersect(rownames(SHH_Group4_reversed_gp),rownames(SHH_Group4_reversed_gp_pvalue))


  SHH_vs_all_reversed_gp_gene<-Reduce(intersect,list(SHH_Group4_reversed_gp_gene, SHH_Group3_reversed_gp_gene,SHH_WNT_reversed_gp_gene))
  ####above is the process for SHH rank reversed analysis



  ################################################################step2 for WNT
  corGenes <- cor(t(ranked_data[WNT_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))


  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_data), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_data))
  WNT_reversed_gp<- t(geneRatioOut)

  ###
  SHH_df<-WNT_reversed_gp[,colnames(WNT_reversed_gp) %in% SHH]
  WNT_df<-WNT_reversed_gp[,colnames(WNT_reversed_gp) %in% WNT]
  Group3_df<-WNT_reversed_gp[,colnames(WNT_reversed_gp) %in% Group3]
  Group4_df<-WNT_reversed_gp[,colnames(WNT_reversed_gp) %in% Group4]
  #for WNT
  WNT_ab <- data.frame(rowSums(WNT_df))
  names(WNT_ab)[1]<-"a"
  WNT_ab$b<-length(colnames(WNT_df)) - WNT_ab$a

  #for SHH
  SHH_cd <- data.frame(rowSums(SHH_df))
  names(SHH_cd)[1]<-"c"
  SHH_cd$d<-length(colnames(SHH_df)) - SHH_cd$c

  #for Group3
  Group3_cd <- data.frame(rowSums(Group3_df))
  names(Group3_cd)[1]<-"c"
  Group3_cd$d<-length(colnames(Group3_df)) - Group3_cd$c

  #for Group4
  Group4_cd <- data.frame(rowSums(Group4_df))
  names(Group4_cd)[1]<-"c"
  Group4_cd$d<-length(colnames(Group4_df)) - Group4_cd$c

  ##########calculate reversal ratio and fisher exact test
  ###################for WNT_SHH
  WNT_SHH<-cbind(WNT_ab,SHH_cd)
  rownames(WNT_SHH)<-rownames(WNT_reversed_gp)
  WNT_SHH$result <- with(WNT_SHH, (a + d) / (a + b + c + d))
  WNT_SHH_reversed_gp<-WNT_SHH[WNT_SHH$result>0.95,]


  df<-data.frame(WNT_SHH[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  WNT_SHH_reversed_gp_pvalue<-df[df$p_value<0.001,]
  WNT_SHH_reversed_gp_gene<-intersect(rownames(WNT_SHH_reversed_gp),rownames(WNT_SHH_reversed_gp_pvalue))

  ###################for WNT_Group3
  WNT_Group3<-cbind(WNT_ab,Group3_cd)
  rownames(WNT_Group3)<-rownames(WNT_reversed_gp)
  WNT_Group3$result <- with(WNT_Group3, (a + d) / (a + b + c + d))
  WNT_Group3_reversed_gp<-WNT_Group3[WNT_Group3$result>0.95,]


  df<-data.frame(WNT_Group3[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  WNT_Group3_reversed_gp_pvalue<-df[df$p_value<0.001,]
  WNT_Group3_reversed_gp_gene<-intersect(rownames(WNT_Group3_reversed_gp),rownames(WNT_Group3_reversed_gp_pvalue))

  ###################for WNT_Group4
  WNT_Group4<-cbind(WNT_ab,Group4_cd)
  rownames(WNT_Group4)<-rownames(WNT_reversed_gp)
  WNT_Group4$result <- with(WNT_Group4, (a + d) / (a + b + c + d))
  WNT_Group4_reversed_gp<-WNT_Group4[WNT_Group4$result>0.95,]


  df<-data.frame(WNT_Group4[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  WNT_Group4_reversed_gp_pvalue<-df[df$p_value<0.001,]
  WNT_Group4_reversed_gp_gene<-intersect(rownames(WNT_Group4_reversed_gp),rownames(WNT_Group4_reversed_gp_pvalue))

  WNT_vs_all_reversed_gp_gene<-Reduce(intersect,list(WNT_Group4_reversed_gp_gene, WNT_Group3_reversed_gp_gene,WNT_SHH_reversed_gp_gene))


  ###################################################step 2 for Group3
  corGenes <- cor(t(ranked_data[Group3_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))



  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_data), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_data))
  Group3_reversed_gp<- t(geneRatioOut)

  ###
  SHH_df<-Group3_reversed_gp[,colnames(Group3_reversed_gp) %in% SHH]
  WNT_df<-Group3_reversed_gp[,colnames(Group3_reversed_gp) %in% WNT]
  Group3_df<-Group3_reversed_gp[,colnames(Group3_reversed_gp) %in% Group3]
  Group4_df<-Group3_reversed_gp[,colnames(Group3_reversed_gp) %in% Group4]
  #for Group3
  Group3_ab <- data.frame(rowSums(Group3_df))
  names(Group3_ab)[1]<-"a"
  Group3_ab$b<-length(colnames(Group3_df)) - Group3_ab$a

  #for SHH
  SHH_cd <- data.frame(rowSums(SHH_df))
  names(SHH_cd)[1]<-"c"
  SHH_cd$d<-length(colnames(SHH_df)) - SHH_cd$c

  #for WNT
  WNT_cd <- data.frame(rowSums(WNT_df))
  names(WNT_cd)[1]<-"c"
  WNT_cd$d<-length(colnames(WNT_df)) - WNT_cd$c

  #for Group4
  Group4_cd <- data.frame(rowSums(Group4_df))
  names(Group4_cd)[1]<-"c"
  Group4_cd$d<-length(colnames(Group4_df)) - Group4_cd$c

  ##########calculate reversal ratio and fisher exact test
  ###################for Group3_SHH
  Group3_SHH<-cbind(Group3_ab,SHH_cd)
  rownames(Group3_SHH)<-rownames(Group3_reversed_gp)
  Group3_SHH$result <- with(Group3_SHH, (a + d) / (a + b + c + d))
  Group3_SHH_reversed_gp<-Group3_SHH[Group3_SHH$result>0.9,]


  df<-data.frame(Group3_SHH[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group3_SHH_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group3_SHH_reversed_gp_gene<-intersect(rownames(Group3_SHH_reversed_gp),rownames(Group3_SHH_reversed_gp_pvalue))
  #3097 gps

  ###################for Group3_WNT
  Group3_WNT<-cbind(Group3_ab,WNT_cd)
  rownames(Group3_WNT)<-rownames(Group3_reversed_gp)
  Group3_WNT$result <- with(Group3_WNT, (a + d) / (a + b + c + d))
  Group3_WNT_reversed_gp<-Group3_WNT[Group3_WNT$result>0.85,]


  df<-data.frame(Group3_WNT[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group3_WNT_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group3_WNT_reversed_gp_gene<-intersect(rownames(Group3_WNT_reversed_gp),rownames(Group3_WNT_reversed_gp_pvalue))
  #4540 gps
  ###################for Group3_Group4
  Group3_Group4<-cbind(Group3_ab,Group4_cd)
  rownames(Group3_Group4)<-rownames(Group3_reversed_gp)
  Group3_Group4$result <- with(Group3_Group4, (a + d) / (a + b + c + d))
  Group3_Group4_reversed_gp<-Group3_Group4[Group3_Group4$result>0.8,]
  #7705 gps

  df<-data.frame(Group3_Group4[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group3_Group4_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group3_Group4_reversed_gp_gene<-intersect(rownames(Group3_Group4_reversed_gp),rownames(Group3_Group4_reversed_gp_pvalue))

  Group3_vs_all_reversed_gp_gene<-Reduce(intersect,list(Group3_SHH_reversed_gp_gene, Group3_Group4_reversed_gp_gene,Group3_WNT_reversed_gp_gene))



  ###########################################step2 for Group4
  corGenes <- cor(t(ranked_data[Group4_rank_t_gene,]))
  corGenes[lower.tri(corGenes)] <- 1
  corGenes <- data.frame(reshape2::melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))

  geneRatioOut <- apply(corGenes,  function(x) createRatio(exprs = as.matrix(ranked_data), x = x), MARGIN=1)
  colnames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  rownames(geneRatioOut) <- colnames(as.matrix(ranked_data))
  Group4_reversed_gp<- t(geneRatioOut)

  ###
  SHH_df<-Group4_reversed_gp[,colnames(Group4_reversed_gp) %in% SHH]
  WNT_df<-Group4_reversed_gp[,colnames(Group4_reversed_gp) %in% WNT]
  Group3_df<-Group4_reversed_gp[,colnames(Group4_reversed_gp) %in% Group3]
  Group4_df<-Group4_reversed_gp[,colnames(Group4_reversed_gp) %in% Group4]
  #for Group4
  Group4_ab <- data.frame(rowSums(Group4_df))
  names(Group4_ab)[1]<-"a"
  Group4_ab$b<-length(colnames(Group4_df)) - Group4_ab$a

  #for SHH
  SHH_cd <- data.frame(rowSums(SHH_df))
  names(SHH_cd)[1]<-"c"
  SHH_cd$d<-length(colnames(SHH_df)) - SHH_cd$c

  #for WNT
  WNT_cd <- data.frame(rowSums(WNT_df))
  names(WNT_cd)[1]<-"c"
  WNT_cd$d<-length(colnames(WNT_df)) - WNT_cd$c

  #for Group3
  Group3_cd <- data.frame(rowSums(Group3_df))
  names(Group3_cd)[1]<-"c"
  Group3_cd$d<-length(colnames(Group3_df)) - Group3_cd$c

  ##########calculate reversal ratio and fisher exact test
  ###################for Group4_SHH
  Group4_SHH<-cbind(Group4_ab,SHH_cd)
  rownames(Group4_SHH)<-rownames(Group4_reversed_gp)
  Group4_SHH$result <- with(Group4_SHH, (a + d) / (a + b + c + d))
  Group4_SHH_reversed_gp<-Group4_SHH[Group4_SHH$result>0.9,]


  df<-data.frame(Group4_SHH[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group4_SHH_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group4_SHH_reversed_gp_gene<-intersect(rownames(Group4_SHH_reversed_gp),rownames(Group4_SHH_reversed_gp_pvalue))
  #3097 gps

  ###################for Group4_WNT
  Group4_WNT<-cbind(Group4_ab,WNT_cd)
  rownames(Group4_WNT)<-rownames(Group4_reversed_gp)
  Group4_WNT$result <- with(Group4_WNT, (a + d) / (a + b + c + d))
  Group4_WNT_reversed_gp<-Group4_WNT[Group4_WNT$result>0.85,]


  df<-data.frame(Group4_WNT[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group4_WNT_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group4_WNT_reversed_gp_gene<-intersect(rownames(Group4_WNT_reversed_gp),rownames(Group4_WNT_reversed_gp_pvalue))
  #4540 gps
  ###################for Group4_Group3
  Group4_Group3<-cbind(Group4_ab,Group3_cd)
  rownames(Group4_Group3)<-rownames(Group4_reversed_gp)
  Group4_Group3$result <- with(Group4_Group3, (a + d) / (a + b + c + d))
  Group4_Group3_reversed_gp<-Group4_Group3[Group4_Group3$result>0.85,]
  #7705 gps

  df<-data.frame(Group4_Group3[,1:4])
  # Calculate the Fisher's Exact Test p-value for each row
  df$p_value <- apply(df, 1, function(row) {
    # Convert the row into a matrix suitable for Fisher's test
    matrix_data <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
    # Conduct Fisher's Exact Test
    test_result <- fisher.test(matrix_data)
    # Return the p-value
    return(test_result$p.value)
  })
  Group4_Group3_reversed_gp_pvalue<-df[df$p_value<0.001,]
  Group4_Group3_reversed_gp_gene<-intersect(rownames(Group4_Group3_reversed_gp),rownames(Group4_Group3_reversed_gp_pvalue))

  Group4_vs_all_reversed_gp_gene<-Reduce(intersect,list(Group4_SHH_reversed_gp_gene, Group4_Group3_reversed_gp_gene,Group4_WNT_reversed_gp_gene))

  all_reversed_gp_genes <- list(
    SHH_ranked_gene = SHH_vs_all_reversed_gp_gene,
    WNT_ranked_gene = WNT_vs_all_reversed_gp_gene,
    Group3_ranked_gene = Group3_vs_all_reversed_gp_gene,
    Group4_ranked_gene = Group4_vs_all_reversed_gp_gene
  )

  return(all_reversed_gp_genes)


}
