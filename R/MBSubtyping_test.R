#' @importFrom stats coef cor fisher.test t.test
#' @importFrom utils install.packages
#' @title predict medulloblastoma (MB) subtype
#' @description
#' this function is used to subtype given samples based on gene ratios we selected from previous steps, for each sample, we predict its subtype
#'
#' @param data the example independent test dataset (all_13datasets) to evaluate the performance of RaMBat
#' @param MB_RANK_GP the unique gene ratios (GERs) selected from the training dataset
#' @param sampAnnot the sampAnnot is the annotation file for each sample in training dataset
#'
#' @return a global object 'Pred_performance'
#' @export
#'
#' @examples Pred_performance<-testMB(all_13datasets, MB_RANK_GP, samp_13)
testMB<-function(data,MB_RANK_GP,sampAnnot){
  medulloGeneSetsUp <- MB_RANK_GP
  signatureProbes <-c(medulloGeneSetsUp$WNT,medulloGeneSetsUp$SHH,
                      medulloGeneSetsUp$Group3,medulloGeneSetsUp$Group4)

  getGenes <- function(x) {
    out <- strsplit(x, split="_")
    output <- c(out[[1]])
  }

  signatureGenes <- unique(as.character(sapply(signatureProbes, FUN=getGenes)))

  #Filter Matrix
  exprs<-data
  exprs_SG <- exprs[intersect(rownames(exprs), signatureGenes),,drop = F]

  #Create Ratios
  createRatio <- function(x) {
    g1 <- x[1]
    g2 <- x[2]
    g1g2_ratio <- 2^(exprs[g1,,drop = F]-exprs[g2,,drop = F])
    return(g1g2_ratio)
  }
  if(ncol(exprs_SG) > 1){
    corGenes <- stats::cor(t(exprs_SG))
    corGenes <- data.frame(reshape2::melt(corGenes))
    corGenes <- corGenes[corGenes[,"value"]<.99,]
    print(dim(corGenes))
  } else {
    corGenes <- expand.grid(rownames(exprs_SG), rownames(exprs_SG))
    corGenes <- corGenes[which(corGenes$Var1 != corGenes$Var2),]
  }

  #exprs <- as.matrix(gse21140)
  exprs <- as.matrix(exprs)
  geneRatioOut <- apply(corGenes, MARGIN = 1, FUN = createRatio)
  if(ncol(exprs) == 1){
    geneRatioOut <- as.data.frame(geneRatioOut)
  } else {
    geneRatioOut <- data.frame(t(geneRatioOut))
  }

  # assign rownames and column names
  rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  colnames(geneRatioOut) <- colnames(exprs)
  geneRatioOut <- geneRatioOut[!is.infinite(rowSums(geneRatioOut)),,drop = F] # added

  ################################
  #Filter to signature ratios
  #Create Heatmap
  ################################

  geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),,drop = F]
  medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
  medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
  medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
  medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))

  calcScore <- function(myMat, mySetsUp) {
    getScoreSet <- function(x, myMat = myMat) {
      return(colMeans(myMat[rownames(myMat) %in% x,,drop = F]))
    }
    myMatUp <- data.frame(lapply(mySetsUp, FUN = getScoreSet, myMat))
    myMatUp$Pred_Subtype <- apply(myMatUp, 1, function(row) {
      colnames(myMatUp)[which.max(row)]
    })
    return( myMatUp)
  }

  myMat <- calcScore(myMat = geneRatioOut, mySetsUp = medulloGeneSetsUp)
  match<-match(samp_13$gsm,rownames(myMat))
  myMat$GroundTruth<-samp_13$Subtype[match]
  pred <- data.frame(cbind(myMat$Pred_Subtype,rownames(myMat)))
  names(pred)<-c("best.fit", "sample")
  rownames(pred)<-pred$sample
  sample.order <- rownames(pred)

  #pred <- pred[sample.order,]


  samp_13<-sampAnnot
  pred<-pred[samp_13$gsm,]# keeping the correct order
  best.fit<-pred[,1]
  #actual<-label_5_LBL202#label_2_validation

  pred<-cbind(samp_13, best.fit)

  data<-pred
  # Filter out rows where 'actual' is 'U'
  filtered_data <- data[data$Subtype != 'U', ]
  # Get unique studies
  unique_studies <- unique(filtered_data$dataset)
  # Initialize a vector to store accuracies
  accuracies <- numeric(length(unique_studies))
  # Loop through each study to calculate accuracy
  for (i in seq_along(unique_studies)) {
    study_name <- unique_studies[i]
    study_data <- filtered_data[filtered_data$dataset == study_name, ]
    accuracies[i] <- 100*mean(study_data$best.fit == study_data$Subtype)
  }
  # Combine the studies and their accuracies into a data frame
  accuracy_per_study <- data.frame(Study = unique_studies, Accuracy = accuracies)
  overal=100*mean(filtered_data$Subtype == filtered_data$best.fit)
  accuracy_per_study[(nrow(accuracy_per_study)+1),]<-c('Overal', overal)
  Pred_performance <- list(
    pred_score = myMat,
    pred_accuracy = accuracy_per_study
  )

  return(Pred_performance)
}
