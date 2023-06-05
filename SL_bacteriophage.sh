#!/bin/bash
#SBATCH --account m342
#SBATCH --qos regular
#SBATCH --nodes 1
#SBATCH --constraint cpu
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH -t 0:10:00
#SBATCH -J bacteriophage

#####################
#####rule genomad:###
#####################
source /global/homes/j/jdyuzon/.bashrc
conda activate snakemake-bacteriophage


sample=$1

### Genomad
mkdir genomad_output/
mkdir genomad_output/${sample}
genomad end-to-end --min-score 0.7 --cleanup --splits 8 bacterial_genomes/${sample}.fna genomad_output/${sample} genomad_db

### CheckV
mkdir checkv_output/
mkdir checkv_output/${sample}/
checkv end_to_end genomad_output/${sample}/${sample}_summary/${sample}_virus.fna checkv_output/${sample} -d checkv-db-v1.5 -t 16
cp checkv_output/${sample}/completeness.tsv checkv_output/${sample}_completeness.tsv

### MMseqs2
conda activate mmseqs2
mkdir mmseqs_target_seq/
mkdir mmseqs_target_seq/${sample}
mkdir phrog_output/
cp genomad_output/${sample}/${sample}_summary/${sample}_virus_proteins.faa mmseqs_target_seq/${sample}/${sample}_virus_proteins.faa 
mmseqs createdb mmseqs_target_seq/${sample}/${sample}_virus_proteins.faa mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq 

### MMseqs2/Phrogs
mmseqs search phrogs_mmseqs_db/phrogs_profile_db \
mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs \
mmseqs_target_seq/${sample}/tmp -s 7

mmseqs createtsv phrogs_mmseqs_db/phrogs_profile_db \
mmseqs_target_seq/${sample}/${sample}_virus_proteins.target_seq \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs \
mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs.tsv --full-header 

cp mmseqs_target_seq/${sample}/${sample}_virus_proteins_mmseqs.tsv mmseqs_target_seq
echo "file: mmseqs_target_seq/${sample}_virus_proteins_mmseqs.tsv"

mkdir plasmid_output

conda activate snakemake-bacteriophage
python3 bacteriophage_edit.py ${sample}

