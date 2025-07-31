# EnsembleAge Clock

## Web service

Our web service is available at: https://ensemble.epiclock.app/

## How to prepare your sample files

Single Sample Upload: 
Upload a csv file with your DNA methylation data for a specific sample. The first column should be the cpg labels that correspond to the DNA methylation values in the file. The second column should be your DNA methylation values. Make sure the top of the column with methylation data values is labeled with the sample name (e.g. Sample 1). Please refer to [individual01.csv](https://github.com/comp-epi-lee-lab/EnsembleAgeClock/blob/main/examples/individual01.csv)as an example. 


![Alt text](examples/singlesampleEX.png)


Multi-Sample Upload: 
Upload a csv file with the DNA methylation data of multiple samples as multiple sequential columns in the file. The first column should be the cpg labels that correspond to the DNA methylation values in the file. The rest of the columns should be the DNA methylation data per sample. Make sure the top of each column is labeled with the sample name (e.g. Sample 1, Sample 2, etc). 


![Alt text](examples/multisampleEX.png)


Note about Grim Age: 
For the single sample upload for grim age, it asks for an input of age and gender. If the user doesn't have this, simply press submit and the predictions for the other clocks will still be displayed. 


