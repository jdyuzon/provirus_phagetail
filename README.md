# provirus_phagetail
Machine Learning modules to distinguish between proviruses and phage tail-like particles from genes or operons identified in bacterial genomes

## 1. Software requirements:
channels:
  - bioconda
  - conda-forge

dependencies:
  - snakemake-minimal >=5.24.1
  - geNomad >=1.5.0
  - mmseqs2 >=13.45111
  - checkv >=1.0.1
  - blast >=2.6.0
  - prokka =1.11
  - pandas
  - bedtools =2.27.1
  - trnascan-se >= 2.0.2
 
Set up the environment:
```
conda env create --name snakemake-bacteriophage --file environment.yaml
```

## 2. config.yaml list of bacterial genomes to query:
```
less config.yaml
```
The current list shows the training and test set of bacterial genomes with phages, T6SS, and eCIS. For example the config file can be configured:

>samples:
>>Pseudomonas_simiae_WCS417_IMG_2585427642
>>Pseudomonas_putida_KT2440_ASM756v2


## 3.SL_bacteriophage.sh shows the pipeline:
Samples are processed in the following order:
Bacterial genomes are annotated and proviruses are identified by genomad<br>
CheckV returns complete and clean phage genomes<br>
MMseqs queries the viral proteins against the Phrog database<br>
Identified proviruses and phage genes are passed to bacteriophage_edit.py<br>

An example command:
```
sbatch SL_bacteriophage.sh Pseudomonas_simiae_WCS417_IMG_2585427642
```

## 4. bacteriophage_edit.py gathers information for each sample from checkv, genomad, mmseqs/phrog into tables

## 5. bacteriophage.sh submits jobs to slurm
The pipeline can also be run on slurm with all the samples/bacterial genomes listed in config.yaml
```
bash bacteriophage.sh
```

## 6. bacteriophage_filter.py compares identified proviruses and phage gene results to known known proviruses and phage tail-like structures in bacterial genomes
gathers checkv, phrog, and plasmid results from all samples into tables<br>
annotated genes from genomad output are put into a bed file format
```
prediction_out/all_genes.phage_tail.bed
```
bed file of annotated genes are intersected with all_genes.phage_tail.bed that have operon information of phages, T6SS, and eCIS. A table listing operons and all gene annotations associated with operon are returned
```
../prediction_out/genes_phage_tail_variables.csv
```

## 7. NaturalLanguageProcessing.ipynb identifies genes as either from provirus or phage tail-like structure using gene annotation descriptions
prediction_out/all_genes.phage_tail.bed annotation descriptions are the X variable and type are the y/target variables. X variables are strings that are standardized, tokenized, and vectorized for the recurrent neural network (RNN). The y variables are coded as binary with phages as 1 and phage tail-like structures as 0. Data is split into train, validation, and test with no overlaps between the datasets. A model with Long short-term memory network (LSTM) layers and sigmoid activations is used.

## 8. RandomForest_phage_trail.ipynb and XGBoost_phage_trail.ipynb identifies operons as either provirus genomes or phage tail-like structures
In RandomForest_phage_trail.ipynb, prediction_out/genes_phage_tail_variables.csv annotations are X variables and y variables are phages or phage tail like structures. Presence of an X variable/annotation is encoded as 1 and absence is 0. y variables are phages encoded as 1 and phage tail-like structures as 0. Data is split into train and test datasets. Model performance is shown as accuracy, precision, recall, and F1 score. Feature importance described by GINI/Mean Decrease in Impurity but better metrics will be implemented.<br>

The approach for XGBoost_phage_trail.ipynb dataset uses two approaches: 1) use the same training and test dataset as Random Forest model, or 2) reduce dataset using the first PCAs that describe most of the variation. The second method also uses optune for hyperperameter optimization.










