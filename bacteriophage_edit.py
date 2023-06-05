from glob import glob
import pandas as pd
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

sample=str(sys.argv[1])
print('sample:', sample)

### checkv_edit
def checkv_edit(sample_id):
    try:
        checkv_table=pd.read_csv("checkv_output/"+sample_id+"_completeness.tsv",sep = '\t', encoding_errors='ignore')
        checkv_table['host']=sample_id
        genomad_table=pd.read_csv("genomad_output/"+sample_id+"/"+sample+"_summary/"+sample_id+"_virus_summary.tsv", sep = '\t', encoding_errors='ignore')
        genomad_table.columns=["contig_id","length","topology","coordinates","n_genes","genetic_code","virus_score","fdr","n_hallmarks"," marker_enrichment","taxonomy"]
        checkv_table=pd.merge(checkv_table,genomad_table, how = 'left', on = 'contig_id').fillna('NA')
        checkv_table.to_csv("checkv_output/"+sample_id+"_completeness.name.tsv", sep='\t')
    except:
        pass
    
checkv_edit(sample)


### phrog_edit
def phrog_edit(file_in,file_out):
    try:
        phrog_table=pd.read_csv(file_in,sep = '\t', encoding_errors='ignore',header=None)
        phrog_table.columns=['#phrog','host_seq','alnScore','seqIdentity','eVal','qStart','qEnd','qLen','tStart','tEnd','tLen']
        phrog_table[['#phrog','phrog_seq']] =phrog_table['#phrog'].str.split(' ## ',expand=True)
        df_index = pd.read_csv('phrogs_mmseqs_db/PHROG_index.csv')
        df_Bins_Index = pd.merge(phrog_table,df_index, how = 'left', on = '#phrog').fillna('NA')
        df_Bins_Index['host']=sample
        np.set_printoptions(threshold=sys.maxsize)
        df_Bins_Index.to_csv(file_out)
    except:
        pass
    
def gene_edit(file_in,file_out):
    try:
        gene_table=pd.read_csv(file_in,sep = '\t', encoding_errors='ignore')
        gene_table['host']=sample
        gene_table.to_csv(file_out, sep='\t')
    except:
        pass

        
#make virus file
phrog_edit("mmseqs_target_seq/"+sample+"_virus_proteins_mmseqs.tsv","phrog_output/"+sample+"_virus_proteins_mmseqs_Annotated.csv")

#make general protein and gene files
phrog_edit("mmseqs_target_seq/"+sample+"_proteins_mmseqs.tsv","phrog_output/"+sample+"_proteins_mmseqs_Annotated.csv")

gene_edit("genomad_output/"+sample+"/"+sample+"_annotate/"+sample+"_genes.tsv","genomad_output/"+sample+"/"+sample+"_annotate/"+sample+"_genes_Annotated.tsv")

### plasmid_edit
plasmid_table=pd.read_csv("genomad_output/"+sample+"/"+sample+"_summary/"+sample+"_plasmid_summary.tsv",sep = '\t', encoding_errors='ignore')
print(sample)
plasmid_table['host']=sample
plasmid_table.to_csv("plasmid_output/"+sample+"_plasmid_summary.name.tsv", sep='\t')











