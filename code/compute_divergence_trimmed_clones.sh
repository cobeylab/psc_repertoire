#!/bin/bash

# Read dataset name from 1st input argument
dataset=$1

echo $dataset    
# ---------------------------- Create sbatch file for job array --------------------------
sbatch_file=compute_divergence_trimmed_$dataset.sbatch

echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=divergence_trimmed_"$dataset >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH -o out_err_files/divergence_trimmed_"$dataset"_%A_%a.out" >> $sbatch_file       
echo "#SBATCH -e out_err_files/divergence_trimmed_"$dataset"_%A_%a.err" >> $sbatch_file         
echo "#SBATCH --time=100:00:00" >> $sbatch_file
echo "#SBATCH --mem-per-cpu=2000" >> $sbatch_file

echo input_file_wholeseq=../results/partis/${dataset}_clones/${dataset}_clone_'${SLURM_ARRAY_TASK_ID}'_alignment_TRIMMED.fasta >> $sbatch_file
echo input_file_cdr3=../results/partis/${dataset}_clones/${dataset}_clone_'${SLURM_ARRAY_TASK_ID}'_CDR3_alignment_TRIMMED.fasta >> $sbatch_file

echo python calculate_aa_divergence.py '$input_file_wholeseq' NAIVE >> $sbatch_file
echo python calculate_aa_divergence.py '$input_file_cdr3' NAIVE >> $sbatch_file

# ----------------------------------------------------------------------------------------


# Get all clone identifiers
CLONE_IDS=$(ls ../results/partis/${dataset}_clones/*alignment_TRIMMED* | grep -v '.csv'| grep -o 'clone_[0-9]*'|tr -d [a-z_] | uniq) 

#Concatenate into single variable with ',' for sbatch array:
CLONE_IDS=$(echo $CLONE_IDS | tr ' ' ',')
    



# Run sbatch command passing retained clone ids to --array option
sbatch --array=$CLONE_IDS $sbatch_file     
# Remove sbatch file
rm $sbatch_file
