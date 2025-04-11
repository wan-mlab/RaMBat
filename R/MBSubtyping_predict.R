#' @importFrom stats coef cor fisher.test t.test
#' @importFrom utils install.packages
#' @title predict medulloblastoma (MB) subtype
#' @description
#' this function is used to subtype given samples. Based on gene ratios we selected from previous steps, for each sample from independent test dataset, we predict its subtype
#'
#' @param data the independent test dataset we need to predit subtype information based on extracted gene ratios (GERs)
#'
#' @param MB_RANK_GP the unique gene ratios (GERs) selected from the training dataset
#'
#' @return a global object 'Predicted_subtype'
#' @export
#'
#' @examples Predicted_subtype<-predMB(data, MB_RANK_GP)
predMB<-function(data,MB_RANK_GP){
  medulloGeneSetsUp <- MB_RANK_GP #MB_RANK_GP#GB_RANK_GP_reversalratio#medulloSetsUp#MB_RANK_GP#MB_RANK_GP_reversalratio_alpha0.06
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

  Predicted_subtype <- calcScore(myMat = geneRatioOut, mySetsUp = medulloGeneSetsUp)


  return(Predicted_subtype)
}
