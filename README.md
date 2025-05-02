# RaMBat : A gene-<ins>Ra</ins>nking based framework for identifying Medulloblasotma (<ins>MB</ins>) subtypes with severe b<ins>at</ins>ch effects

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
1. use GRA() function to extract differentially ranked genes (gene rank analysis)
all_rank_t_genes<- GRA(data, sampAnnot)

2. use RRA() function to extract gene ratios (GERs) (Reversal ratio analysis)
all_reversed_gp_genes<-RRA(data, sampAnnot, all_rank_t_genes)

3. use LaSelect() function to select non-zero coefficient GERs (LASSO feature selection)
MB_RANK_GP<-LaSelect(data, sampAnnot, all_rank_t_genes,all_reversed_gp_genes)

4. use predMB() function to predict the MB subtype information (MB subtype identification)
Predicted_subtype<- predMB(test_data,MB_RANK_GP)
```

## Predict the MB subtype information directly
```bash
data(MB_RANK_GP)
Predicted_subtype<-preMB(test_data,MB_RANK_GP)
```

## Example
1. Extrat features from training dataset (GSE85217) then predict the subtype information for test dataset (GSE21140)
```bash
data(GSE85217)
data(sampAnnot_GSE85217)
1. all_rank_t_genes<- GRA(GSE85217, sampAnnot_GSE85217)
2. all_reversed_gp_genes<-RRA(GSE85217, sampAnnot_GSE85217, all_rank_t_genes)
3. MB_RANK_GP<-LaSelect(GSE85217, sampAnnot_GSE85217, all_rank_t_genes,all_reversed_gp_genes)
4. myMat<- predMB(GSE21140,MB_RANK_GP)
```
2. Predict the subtype information for given dataset based on features selected from training dataset (GSE85217)
```bash
data(MB_RANK_GP)
myMat<-preMB(GSE21140, MB_RANK_GP)
```
## Evaluate the performance of RaMBat
Load the necessary data. **all_13datasets** is the normalized combination of all 13 independent test datasets. **samp_13** is the annotation file for all_13datasets
```bash
data(MB_RANK_GP)
data(samp_13)
data(all_13datasets)
Pred_performance<-testMB(all_13datasets, MB_RANK_GP, samp_13)
```
