#!/bin/bash
# Runs igphyml analysis for a single clone (EXCLUDES CDR3)

dataset=$1
clone_id=$2

sbatch_file=run_igphyml_${dataset}_clone_${clone_id}.sbatch

seqfile=../results/partis/${dataset}_igphyml/${dataset}_clone_${clone_id}_noCDR3_igphyml.fasta
partfile=../results/partis/${dataset}_igphyml/${dataset}_clone_${clone_id}_noCDR3_part.txt

echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=run_igphyml_${dataset}_clone_${clone_id}" >> $sbatch_file
echo "#SBATCH --output=out_err_files/run_igphyml"_${dataset}_"clone_"$clone_id.out >> $sbatch_file
echo "#SBATCH --error=out_err_files/run_igphyml"_${dataset}_"clone_"$clone_id.err >> $sbatch_file
echo "#SBATCH --time=5:00:00" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=1" >> $sbatch_file

echo igphyml -i $seqfile -m HLP --root NAIVE --partfile $partfile --oformat tab --omega ce,ce  >> $sbatch_file

echo "mkdir -p ../results/igphyml/${dataset}" >> $sbatch_file

# Pull tree out of igphyml stat.tab file
echo "tail -1 ${seqfile}_igphyml_stats.tab | grep -o '(.*;' > ../results/igphyml/${dataset}/${dataset}_clone_${clone_id}_noCDR3_igphyml.tree" >> $sbatch_file
 
echo mv ${seqfile}_igphyml_stats.tab ../results/igphyml/${dataset}/${dataset}_clone_${clone_id}_noCDR3_igphyml_stats.tab >> $sbatch_file
echo rm ${seqfile}_igphyml_CIlog.txt >> $sbatch_file


sbatch $sbatch_file
rm $sbatch_file

