# A-Scale-free-Approach-for-False-Discovery-Rate-Control-in-Generalized-Linear-Models
Reproduce all the results for the paper `` A Scale-free Approach for False Discovery Rate Control in Generalized Linear Models''

# Data

**Abstract**

We use two real data sets in our paper. The first one is an scRNAseq data set, which contains expressioons of a total of 400 T47D A1-2 human breast cancer cells. This data set contains 400 samples in the control group, in which the cells are vehicle-treated, and 400 samples in the treatment group, in which the cells are treated with 100 nM synthetic glucocorticoid dexamethasone (Dex). After proper normalization, the final scRNAseq data consist of 800 samples, each with  expressions of 32,049 genes.

**Availability**

The scRNAseq data is available at
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834

**Description**

The datasets are explained in Section 5.3 of the paper.

# Code

**Abstract**

The code to implement our method is summarized [here](https://github.com/Jeremy690/-A-Scale-free-Approach-for-False-Discovery-Rate-Control-in-Generalized-Linear-Models/tree/main/code/functions) in this repository.

**Description

The code folder contains all the scripts we used to reproduce the results in the main text of the paper. It consists of two parts: functions and simulations. The functions folder contains the implementation of our methods and the methods considered for comparisons. It also contains some utility functions that are used to calculate the fdp and power. The simulation folder contains all the codes needed to reproduce the simulations in the main text of the paper. Each sub-folder of the simulation fold contains the R code, python code and the sh file. The sh file is used to run the simulations in a parallel way via the Odyssey system at Harvard. 

# Instructions for Use

Please see the repo [here](https://github.com/Jeremy690/False-Discovery-Rate-via-Data-Splitting).
