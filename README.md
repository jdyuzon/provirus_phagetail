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

## 6. bacteriophage_filter.py 

```

```



