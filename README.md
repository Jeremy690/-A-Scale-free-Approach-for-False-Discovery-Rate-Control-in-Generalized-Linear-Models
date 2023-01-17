# A-Scale-free-Approach-for-False-Discovery-Rate-Control-in-Generalized-Linear-Models
Reproduce all the simulations for the paper "A Scale-free Approach for False Discovery Rate Control in Generalized Linear Models"

# Data

**Abstract**

The real data set in our paper is an scRNAseq data set, which contains expressions of a total of 400 T47D A1-2 human breast cancer cells. This data set contains 400 samples in the control group, in which the cells are vehicle-treated, and 400 samples in the treatment group, in which the cells are treated with 100 nM synthetic glucocorticoid dexamethasone (Dex). After proper normalization, the final scRNAseq data consists of 800 samples, each with  expressions of 32,049 genes.

**Availability**

The scRNAseq data is available at
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834

**Description**

The dataset is explained in Section 5.3 of the paper.

# Code

**Abstract**

The code to implement our method is summarized [here](https://github.com/Jeremy690/-A-Scale-free-Approach-for-False-Discovery-Rate-Control-in-Generalized-Linear-Models/tree/main/code) in this repository.

**Description

The code folder contains two sub folders: functions and simulations. The functions folder contains utils.R, in which there are some common utility functions used to make selections and calculate the fdp and the power. These utility functions will be used in the reproduction of every figure in our paper. The simulations folder contains all the codes needed to reproduce the simulations in the main text of the paper. Every method is separated into its own python or R file



# Instructions for Use

Please see the repo [here](https://github.com/Jeremy690/False-Discovery-Rate-via-Data-Splitting).
