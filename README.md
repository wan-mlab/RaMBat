# RaMBat : A gene-ranking based framework for identifying Medulloblasotma (MB) subtypes with severe batch effects 

## Workflow
![Workflow of RaMBat](workflow.png)

## Installation
To install the package in RStudion:
```bash
devtools::install_github('wan-mlab/RaMBat')
library(RaMBat)
```

##Train your own data to extract features
```bash
1. use GRA() function to extract differentially ranked genes
all_rank_t_genes<- RCA(data, sampAnnot)

2. use RRA() function to extract gene ratios (GERs):
all_reversed_gp_genes<-RRA(data, sampAnnot, all_rank_t_genes)

3. use LaSelect() function to select non-zero coefficient GERs:
MB_RANK_GP<-LaSelect(data, sampAnnot, all_rank_t_genes,all_reversed_gp_genes)

4. use predMB() function to predict the MB subtype information
myMat<- (data)


##Predict the MB subtype information directly
data(MB_RANK_GP)
myMat<-preMB(data)

