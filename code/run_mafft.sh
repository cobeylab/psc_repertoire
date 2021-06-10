#!/bin/bash

# Read dataset name from 1st input argument
dataset=$1

# ---------------------------- Create sbatch file for job array --------------------------
sbatch_file=run_macse_$dataset.sbatch

echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=align_"$dataset >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH -o out_err_files/align_"$dataset"_%A_%a.out" >> $sbatch_file       
echo "#SBATCH -e out_err_files/align_"$dataset"_%A_%a.err" >> $sbatch_file         
echo "#SBATCH --time=100:00:00" >> $sbatch_file
echo "#SBATCH --mem-per-cpu=2000" >> $sbatch_file

echo module load mafft/7.310  >> $sbatch_file


echo INPUT_FILE=../results/partis/${dataset}_clones/${dataset}_clone_'${SLURM_ARRAY_TASK_ID}'.fasta >> $sbatch_file
echo OUTPUT_FILE=../results/partis/${dataset}_clones/${dataset}_clone_'${SLURM_ARRAY_TASK_ID}'_alignment.fasta >> $sbatch_file

echo mafft --op 5 --thread -1 '$INPUT_FILE' '>' '$OUTPUT_FILE' >> $sbatch_file


# ----------------------------------------------------------------------------------------


# Align only clones with at least MINSEQ sequences:
MINSEQ=4

# Get all clone identifiers
ls ../results/partis/${dataset}_clones/${dataset}_clone_* | grep -v 'alignment' | grep -v '.csv'| grep -o 'clone_[0-9]*'|tr -d [a-z_] > clone_ids_${dataset}.tmp

# Split clone ids in files with 500 ids (so that sbatch commands do not become too long)
split -l 500 clone_ids_${dataset}.tmp "clone_ids_${dataset}_set_" -a 1
rm clone_ids_${dataset}.tmp

# Loop over aggregated clone id files, submitting array jobs
for f in clone_ids_${dataset}*
do
    #Get clone identifiers, concatenate into single variable with ',' for sbatch array:
    CLONE_IDS=$(cat $f | tr '\n' ',')
    
    # Remove trailing comma
    CLONE_IDS=${CLONE_IDS::-1}
    
    # Retain only clones with at least MINSEQ sequences
    RETAINED_IDS=''
    
    # For each clone id
    for ID in ${CLONE_IDS//,/ }
    do
      # Number of sequences in clone
      NUMSEQ=$(grep '>' ../results/partis/${dataset}_clones/${dataset}_clone_${ID}.fasta | wc -l)

      # Naive sequence doesn't count towards minimum number of sequences
      NUMSEQ=$(expr $NUMSEQ - $'1')
      
      # Retain if clone has minimum number of sequences
      if [ "$NUMSEQ" -ge "$MINSEQ" ]
        then
        RETAINED_IDS=${RETAINED_IDS},$ID
      
      fi
 
    done

    # Remove leading comma
    RETAINED_IDS=${RETAINED_IDS:1} 
    
    # If at least one clone was retained, align
    if [ ! -z "$RETAINED_IDS" ]
    then
      # Run sbatch command passing retained clone ids to --array option
      sbatch --array=$RETAINED_IDS $sbatch_file     
    fi
    
done

# Remove temporary files
rm clone_ids_${dataset}*
# Remove sbatch file
rm $sbatch_file
