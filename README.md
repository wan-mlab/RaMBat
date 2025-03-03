# RaMBat : A gene-_r̲a̲n̲k̲i̲n̲g̲_ based framework for identifying Medulloblasotma (<ins>MB</ins>) subtypes with severe _batch_ effects

## Workflow
![Workflow of RaMBat](workflow.png)

## Installation
To install the package in RStudion:
```bash
devtools::install_github('wan-mlab/RaMBat')
library(RaMBat)
```

## Train your own data to extract features
```bash
1. use GRA() function to extract differentially ranked genes
all_rank_t_genes<- RCA(data, sampAnnot)

2. use RRA() function to extract gene ratios (GERs):
all_reversed_gp_genes<-RRA(data, sampAnnot, all_rank_t_genes)

3. use LaSelect() function to select non-zero coefficient GERs:
MB_RANK_GP<-LaSelect(data, sampAnnot, all_rank_t_genes,all_reversed_gp_genes)

4. use predMB() function to predict the MB subtype information
myMat<- predMB(data)
```

## Predict the MB subtype information directly
```bash
data(MB_RANK_GP)
myMat<-preMB(data)
```

## Example
1. Extrat features from GSE85217 then predict the subtype information for give data
```bash
data(GSE85217)
data(sampAnnot_GSE85217)
1. all_rank_t_genes<- RCA(GSE85217, sampAnnot_GSE85217)
2. all_reversed_gp_genes<-RRA(GSE85217, sampAnnot_GSE85217, all_rank_t_genes)
3.MB_RANK_GP<-LaSelect(GSE85217, sampAnnot_GSE85217, all_rank_t_genes,all_reversed_gp_genes)
4. myMat<- predMB(data)
```
2. Predict the subtype information for give data directly
```bash
data(MB_RANK_GP)
myMat<-preMB(GSE85217)
```
## Evaluate the performance of RaMBat
```bash
data(MB_RANK_GP)
data(samp_13)
data(all_13datasets)
pred_result<-testMB(all_13datasets, MB_RANK_GP, samp_13)
```
