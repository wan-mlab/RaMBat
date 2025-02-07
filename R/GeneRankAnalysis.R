#' @importFrom stats coef cor fisher.test t.test
#' @importFrom utils install.packages
#' @title extract differetially ranked genes for each subtype
#' @description
#' this function is used to extratc differetially ranked genes from dataset for each subtype
#'
#' @param data the original dataset we extract features
#' @param sampAnnot sampAnnot is the annotation file for each samples in data
#'
#' @return a global object all_rank_t_genes
#' @export
#'
#' @examples all_rank_t_genes<- RCA(GSE85217, sampAnnote_GSE85217)
GRA <- function(data,sampAnnot){
  if (!requireNamespace("purrr", quietly = TRUE)) {
    install.packages("purrr")
  }
  library('purrr')
  GSE85217<-data #the data format should be: each row is a sample, each column is a gene
  sampAnnot_GSE85217<-sampAnnot

  ranked_85217 <- apply(GSE85217, 2, function(x) rank(x, ties.method = "average"))

  SHH<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="SHH",1]
  WNT<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="WNT",1]
  Group3<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="Group3",1]
  Group4<-sampAnnot_GSE85217[sampAnnot_GSE85217$Subtype=="Group4",1]


  SHH_df<-ranked_85217[,colnames(ranked_85217) %in% SHH]
  WNT_df<-ranked_85217[,colnames(ranked_85217) %in% WNT]
  Group3_df<-ranked_85217[,colnames(ranked_85217) %in% Group3]
  Group4_df<-ranked_85217[,colnames(ranked_85217) %in% Group4]


  ########calculate rank change
  SHH_rank<-data.frame(rowMeans(SHH_df))
  names(SHH_rank)[1]<-"SHH_rank_average"
  WNT_rank<-data.frame(rowMeans(WNT_df))
  names(WNT_rank)[1]<-"WNT_rank_average"
  Group3_rank<-data.frame(rowMeans(Group3_df))
  names(Group3_rank)[1]<-"Group3_rank_average"
  Group4_rank<-data.frame(rowMeans(Group4_df))
  names(Group4_rank)[1]<-"Group4_rank_average"
  rank_change<-cbind(SHH_rank, WNT_rank, Group3_rank, Group4_rank)

  ########################################SHH rank change
  SHH_VS_WNT<-abs(data.frame(round(rank_change$SHH_rank_average - rank_change$WNT_rank_average)))
  rownames(SHH_VS_WNT)<-rownames(rank_change)
  names(SHH_VS_WNT)[1]<-"change"
  rows_with_high_change <- SHH_VS_WNT$change > 2314###############11% of total genes
  SHH_VS_WNT_gene <- rownames(SHH_VS_WNT)[rows_with_high_change]

  SHH_VS_Group3<-abs(data.frame(rank_change$SHH_rank_average - rank_change$Group3_rank_average))
  rownames(SHH_VS_Group3)<-rownames(rank_change)
  names(SHH_VS_Group3)[1]<-"change"
  rows_with_high_change <- SHH_VS_Group3$change > 2314
  SHH_VS_Group3_gene <- rownames(SHH_VS_Group3)[rows_with_high_change]

  SHH_VS_Group4<-abs(data.frame(rank_change$SHH_rank_average - rank_change$Group4_rank_average))
  rownames(SHH_VS_Group4)<-rownames(rank_change)
  names(SHH_VS_Group4)[1]<-"change"
  rows_with_high_change <- SHH_VS_Group4$change > 2314
  SHH_VS_Group4_gene <- rownames(SHH_VS_Group4)[rows_with_high_change]

  SHH_specific_rank_gene <- Reduce(intersect,list(SHH_VS_WNT_gene, SHH_VS_Group3_gene, SHH_VS_Group4_gene))#710genes

  ########################################WNT rank change
  WNT_VS_Group3<-abs(data.frame(rank_change$WNT_rank_average - rank_change$Group3_rank_average))
  rownames(WNT_VS_Group3)<-rownames(rank_change)
  names(WNT_VS_Group3)[1]<-"change"
  rows_with_high_change <- WNT_VS_Group3$change > 2314
  WNT_VS_Group3_gene <- rownames(WNT_VS_Group3)[rows_with_high_change]


  WNT_VS_Group4<-abs(data.frame(rank_change$WNT_rank_average - rank_change$Group4_rank_average))
  rownames(WNT_VS_Group4)<-rownames(rank_change)
  names(WNT_VS_Group4)[1]<-"change"
  rows_with_high_change <- WNT_VS_Group4$change > 2314
  WNT_VS_Group4_gene <- rownames(WNT_VS_Group4)[rows_with_high_change]

  WNT_specific_rank_gene<-Reduce(intersect,list(SHH_VS_WNT_gene, WNT_VS_Group3_gene, WNT_VS_Group4_gene)) #880genes

  ########################################Group3 rank change
  Group4_VS_Group3<-abs(data.frame(rank_change$Group4_rank_average - rank_change$Group3_rank_average))
  rownames(Group4_VS_Group3)<-rownames(rank_change)
  names(Group4_VS_Group3)[1]<-"change"
  rows_with_high_change <- Group4_VS_Group3$change > 2314
  Group4_VS_Group3_gene <- rownames(Group4_VS_Group3)[rows_with_high_change]

  Group3_specific_rank_gene<-Reduce(intersect,list(SHH_VS_Group3_gene, WNT_VS_Group3_gene, Group4_VS_Group3_gene)) #346genes

  ########################################Group4 rank change
  Group4_specific_rank_gene<-Reduce(intersect,list(SHH_VS_Group4_gene, WNT_VS_Group4_gene, Group4_VS_Group3_gene)) #345genes

  #######################t-test for each subtype
  SHH_WNT_t<-cbind(SHH_df, WNT_df)
  SHH_Group3_t<-cbind(SHH_df, Group3_df)
  SHH_Group4_t<-cbind(SHH_df, Group4_df)
  WNT_Group3_t<-cbind(WNT_df, Group3_df)
  WNT_Group4_t<-cbind(WNT_df, Group4_df)
  Group3_Group4_t<-cbind(Group3_df, Group4_df)

  ###################SHH_T_TEST
  # Perform t-tests for each gene and identify significant ones
  significant_genes <- apply(SHH_WNT_t, 1, function(x) {
    shh_values <- x[1:length(colnames(SHH_df))]
    wnt_values <- x[(1+length(colnames(SHH_df))):length(colnames(SHH_WNT_t))]
    t_test_result <- t.test(shh_values, wnt_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })

  SHH_WNT_t_genes <- rownames(SHH_WNT_t)[significant_genes]
  #########
  significant_genes <- apply(SHH_Group3_t, 1, function(x) {
    shh_values <- x[1:length(colnames(SHH_df))]
    group3_values <- x[(1+length(colnames(SHH_df))):length(colnames(SHH_Group3_t))]
    t_test_result <- t.test(shh_values, group3_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })

  SHH_Group3_t_genes <- rownames(SHH_Group3_t)[significant_genes]
  ########
  significant_genes <- apply(SHH_Group4_t, 1, function(x) {
    shh_values <- x[1:length(colnames(SHH_df))]
    group4_values <- x[(1+length(colnames(SHH_df))):length(colnames(SHH_Group4_t))]
    t_test_result <- t.test(shh_values, group4_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })
  SHH_Group4_t_genes <- rownames(SHH_Group4_t)[significant_genes]

  SHH_t_specific_genes <- Reduce(intersect, list(SHH_WNT_t_genes, SHH_Group3_t_genes, SHH_Group4_t_genes))

  ###################WNT_T_TEST
  significant_genes <- apply(WNT_Group3_t, 1, function(x) {
    wnt_values <- x[1:length(colnames(WNT_df))]
    group3_values <- x[(1+length(colnames(WNT_df))):length(colnames(WNT_Group3_t))]
    t_test_result <- t.test(wnt_values, group3_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })
  WNT_Group3_t_genes <- rownames(WNT_Group3_t)[significant_genes]

  significant_genes <- apply(WNT_Group4_t, 1, function(x) {
    wnt_values <- x[1:length(colnames(WNT_df))]
    group4_values <- x[(1+length(colnames(WNT_df))):length(colnames(WNT_Group4_t))]
    t_test_result <- t.test(wnt_values, group4_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })
  WNT_Group4_t_genes <- rownames(WNT_Group4_t)[significant_genes]

  WNT_t_specific_genes <- Reduce(intersect, list(SHH_WNT_t_genes, WNT_Group3_t_genes, WNT_Group4_t_genes))

  ###################Group3_T_TEST
  significant_genes <- apply(Group3_Group4_t, 1, function(x) {
    group3_values <- x[1:length(colnames(Group3_df))]
    group4_values <- x[(1+length(colnames(Group3_df))):length(colnames(Group3_Group4_t))]
    t_test_result <- t.test(group3_values, group4_values)
    if (t_test_result$p.value < 0.05) return(TRUE)
    else return(FALSE)
  })
  Group3_Group4_t_genes <- rownames(Group3_Group4_t)[significant_genes]
  Group3_t_specific_genes <- Reduce(intersect, list(SHH_Group3_t_genes, WNT_Group3_t_genes, Group3_Group4_t_genes))

  Group4_t_specific_genes <- Reduce(intersect, list(SHH_Group4_t_genes, WNT_Group4_t_genes, Group3_Group4_t_genes))

  ########################################################select differentially ranked genes for each subtype
  SHH_rank_t_gene<-intersect(SHH_specific_rank_gene,SHH_t_specific_genes)#710 genes
  WNT_rank_t_gene<-intersect(WNT_specific_rank_gene,WNT_t_specific_genes)#880 genes
  Group3_rank_t_gene<-intersect(Group3_specific_rank_gene,Group3_t_specific_genes)#346 genes
  Group4_rank_t_gene<-intersect(Group4_specific_rank_gene,Group4_t_specific_genes)#345 genes

  all_rank_t_genes <- list(
    SHH_ranked_gene = SHH_rank_t_gene,
    WNT_ranked_gene = WNT_rank_t_gene,
    Group3_ranked_gene = Group3_rank_t_gene,
    Group4_ranked_gene = Group4_rank_t_gene
  )

  return(all_rank_t_genes)
}


