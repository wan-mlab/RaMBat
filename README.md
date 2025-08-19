# RaMBat : A gene-<ins>Ra</ins>nking based framework for identifying Medulloblasotma (<ins>MB</ins>) subtypes with severe b<ins>at</ins>ch effects

## Publication
Accurate identification of medulloblastoma subtypes from diverse data sources with severe batch effects by RaMBat. Mengtao Sun, Jieqiong Wang and Shibiao Wan*.  BioRxiv 2025.02.24.640010; DOI: https://doi.org/10.1101/2025.02.24.640010

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
1. Extrat features from your own training dataset (training_dataset), then predict the subtype information for the new test MB samples (new_MB_samples). in this part, a samp_Annot file need to be prepared by the user (check the format of 'sampAnnot_GSE85217' in data folder).
```bash
1. all_rank_t_genes<- GRA(training_dataset, samp_Annot)
2. all_reversed_gp_genes<-RRA(training_dataset, samp_Annot, all_rank_t_genes)
3. MB_RANK_GP<-LaSelect(training_dataset, samp_Annot, all_rank_t_genes,all_reversed_gp_genes)
4. myMat<- predMB(new_MB_samples,MB_RANK_GP)
```
2. Predict the subtype information for given new test MB sample (new_MB_samples) based on features selected by RaMBat from  dataset GSE85217.
```bash
data(MB_RANK_GP)
myMat<-preMB(new_MB_samples, MB_RANK_GP)
```
## Evaluate the performance of RaMBat
Load the necessary data. **all_13datasets** is the normalized combination of all 13 independent test datasets. **samp_13** is the annotation file for all_13datasets
```bash
data(MB_RANK_GP)
data(samp_13)
data(all_13datasets)
Pred_performance<-testMB(all_13datasets, MB_RANK_GP, samp_13)
```
