#!/usr/bin/env python
# coding: utf-8

from glob import glob
import pandas as pd
import numpy as np
import sys
import os
import pybedtools
import importlib

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
norm_e = StandardScaler()

np.set_printoptions(threshold=sys.maxsize)
from functools import reduce
##from pybedtools import BedTool
##pybedtools.set_bedtools_path('~/anaconda3/pkgs/bedtools-2.26.0-0/bin/')
##pybedtools = importlib.reload(pybedtools)


### set work directory
os.chdir("/pscratch/sd/j/jdyuzon/snakemake-bacteriophage2")

### Gather checkv, phrog, and plasmid results
## predicted phage genomes
os.system("cat checkv_output/*_completeness.name.tsv > checkv_output/all_genome_completeness.tsv")
os.system("cat checkv_output/*/complete_genomes.tsv > checkv_output/all_complete_genomes.tsv")
os.system("cat phrog_output/*_virus_proteins_mmseqs_Annotated.csv > phrog_output/all_virus_proteins.csv")
## predicted mobile element genes & proteins (virus/plasmid)
os.system("cat phrog_output/*_proteins_mmseqs_Annotated.csv > phrog_output/all_proteins.csv")
os.system("cat genomad_output/*/*_annotate/*_genes_Annotated.tsv> phrog_output/all_genes.csv")
## predicted plasmids
os.system("cat plasmid_output/*_plasmid_summary.name.tsv > plasmid_output/all_plasmids.tsv")


### Format file of annotated genes from genomad
genes=pd.read_csv("phrog_output/all_genes.csv", sep='\t')
genes['#contig_gene']=genes['gene'].str.replace('\..*','')
genes=genes.loc[(genes['gene'] != 'gene') & (genes['gc_content'] != 'gc_content')]
proteins=pd.read_csv("phrog_output/all_proteins.csv")
proteins=proteins.loc[(proteins['#phrog'] != '#phrog') & (proteins['host_seq'] != 'host_seq') & ~proteins['host_seq'].str.contains("provirus_")]
cols = list(genes)
cols.insert(0, cols.pop(cols.index('#contig_gene')))
cols.insert(1, cols.pop(cols.index('start')))
cols.insert(2, cols.pop(cols.index('end')))
genes = genes.loc[:, cols]
genes=genes.rename(columns = {'start':'start_gene'})
genes=genes.rename(columns = {'end':'stop_gene'})
genes=genes.drop(['Unnamed: 0'], axis=1)
genes.to_csv("prediction_out/all_genes.bed",sep='\t',index=False,header=True)


### Gather Tables of TIGER/Phage1&2, eCIS, T6SS coordinates
os.system("mkdir prediction_out")
os.system("mkdir prediction_out/random_forest")

### Format bedfile of phages and tails to remove the extension on the contig name
phage_tail=pd.read_csv("prediction_out/phage_tail.bed", sep='\t',header=None)
phage_tail[0]=phage_tail[0].str.replace('\..*','',regex=True)
phage_tail[4]=phage_tail[3].str.replace('_.*','',regex=True)
phage_tail[4]=phage_tail[4].str.replace('T6SS.*','T6SS',regex=True)
phage_tail.columns=['#contig_db','start_db','stop_db','ID','type']
phage_tail.to_csv("prediction_out/phage_tail.bed",sep='\t',index=False,header=True)
phage_tail['type'].value_counts()


### intersect Phage/T6SS/eCIS/Tailocins table with genomad_annotation_genes
os.system("bedtools-2.27.1-hd03093a_6/bin/sortBed -i prediction_out/all_genes.bed -header> prediction_out/all_genes.sort.bed ")
os.system("bedtools-2.27.1-hd03093a_6/bin/sortBed -i prediction_out/phage_tail.bed -header> prediction_out/phage_tail.sort.bed ")
os.system("bedtools-2.27.1-hd03093a_6/bin/intersectBed -wa -wb -a prediction_out/phage_tail.sort.bed -b prediction_out/all_genes.bed -header -sortout|sed '1 s_$_\tcontiggene\tstartgene\tstopgene\tgene\tlength\tstrand\tgccontent\tgeneticcode\trbsmotif\tmarker\tevalue\tbitscore\tuscg\tplasmidhallmark\tvirushallmark\ttaxid\ttaxname\tannotationconjscan\tannotationamr\tannotationaccessions\tannotationdescription\thost_'> prediction_out/all_genes.phage_tail.bed")


### output should have coordinates ID, type, gene, annotation (Pfam,COG,etc), host, virus_taxa
genes_phage_tail=pd.read_csv("prediction_out/all_genes.phage_tail.bed", sep='\t')
genes_phage_tail_col=genes_phage_tail.columns
new_col=genes_phage_tail['annotationaccessions'].str.split(';',expand=True).columns
genes_phage_tail[new_col]=genes_phage_tail['annotationaccessions'].str.split(';',expand=True)
genes_phage_tail
### create a new line for each annotation (Pfam,COG,etc) into long format
genes_phage_tail=genes_phage_tail.melt(id_vars=genes_phage_tail_col,ignore_index=False).reset_index().sort_values(['#contig_db', 'annotationaccessions'], ascending=[True, True])
genes_phage_tail=genes_phage_tail[genes_phage_tail["value"].str.contains("None")==False]
genes_phage_tail=genes_phage_tail[genes_phage_tail["value"].str.contains("NaN")==False]
genes_phage_tail=genes_phage_tail.rename(columns={'value': 'annotation'})
genes_phage_tail['variable']=1
genes_phage_tail['viral_ID']=genes_phage_tail['ID']+':'+genes_phage_tail['#contig_db']+':'+genes_phage_tail['start_db'].astype(str)+'_'+genes_phage_tail['start_db'].astype(str)
### reshape dataframe into wide format based on viral_ID (Phage/T6SS/tailocin/eCIS) for each row
genes_phage_tail_variables=genes_phage_tail.reset_index().pivot_table(index=['viral_ID','type'], columns='annotation', values='variable').fillna(0)
genes_phage_tail_variables=genes_phage_tail_variables.reset_index()
genes_phage_tail_variables.to_csv('prediction_out/genes_phage_tail_variables.csv')







