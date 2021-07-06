#!/bin/bash
fasta_file=$1
isotype_info_file=$2


file_id=$(echo $fasta_file | grep -o '[^\/]*.fasta' | cut -d '.' -f 1)

# Automatically generated igblast output file path
igblast_output_file=${file_id}_igblast.fmt7 

# Automatically generated file path with makedb output (igblast run processed by changeO)
makedb_output_file=${file_id}_igblast_db-pass.tab


sbatch_file=run_partis_${file_id}.sbatch


echo '#!/bin/bash' > $sbatch_file
echo "#SBATCH --job-name=run_partis"_$file_id >> $sbatch_file
echo "#SBATCH --output=out_err_files/run_partis"_$file_id.out >> $sbatch_file
echo "#SBATCH --error=out_err_files/run_partis"_$file_id.err >> $sbatch_file
echo "#SBATCH --time=200:00:00" >> $sbatch_file
echo "#SBATCH --partition=cobey" >> $sbatch_file
echo "#SBATCH --nodes=1" >> $sbatch_file
echo "#SBATCH --ntasks-per-node=16" >> $sbatch_file


# Run partis partition:
echo "/project2/cobey/partis/bin/partis partition --n-procs 15 --infname" $fasta_file "--outfname ../results/partis/"$file_id".yaml" "--extra-annotation-columns regional_bounds:cdr3_seqs:cdr3_seqs:seqs_aa:naive_seq_aa" >> $sbatch_file

# Igblast run to identify FRs and CDRs.
echo "module load python/anaconda-2020.02" >> $sbatch_file
echo "AssignGenes.py igblast -s" $fasta_file "-b ~/share/igblast --organism human --loci ig --format blast --outdir ." >> $sbatch_file

# Processing Igblast output into presto-format table
echo "MakeDb.py igblast -i" $igblast_output_file "-s" $fasta_file "-r ~/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta ~/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta --extended --outdir ../results/partis/" >> $sbatch_file

echo "rm" $igblast_output_file  >> $sbatch_file # Removes raw igphyml output

# Processes partis output (partitioning clones into separate fasta files and preparing igphyml input)
echo "module load R/3.6.1" >> $sbatch_file
echo "Rscript process_partis_yaml.R ../results/partis/"$file_id".yaml" "../results/partis/"$makedb_output_file $isotype_info_file >> $sbatch_file


#echo "rm -r _output" >> $sbatch_file

# Run and then remove sbatch file
sbatch $sbatch_file
rm $sbatch_file
