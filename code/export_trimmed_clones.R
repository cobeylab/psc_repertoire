# Using pre-generated alignments, exports alignments for selecting clones that have been trimmed
# (sequences that don't share the light-chain genes of the majority)

library(dplyr)
library(readr)
library(stringr)
library(seqinr)

selected_clones <- read_csv('selected_trimmed_clones.csv')

get_trimmed_clones <- function(dataset, selected_clones){
  dataset_selected_seqs <- selected_clones %>% filter(sample_code == !!dataset) %>%
    select(sample_code, clone, seq_id)
  
  for(clone in unique(dataset_selected_seqs$clone)){
    print(c(dataset, clone))
    
    clone_selected_seqs <- selected_clones %>% filter(sample_code == !!dataset, clone == !!clone) %>% pull(seq_id)
    
    clone_alignment_whole_seq = read.fasta(paste0('../results/partis/', dataset, '_HC_clones/',
                                              dataset, '_HC_clone_', clone, '_alignment.fasta'))
    clone_alignment_CDR3 = read.fasta(paste0('../results/partis/', dataset, '_HC_clones/',
                                             dataset, '_HC_clone_', clone, '_CDR3_alignment.fasta'))
    
    trimmed_clone_whole_seq <- clone_alignment_whole_seq[names(clone_alignment_whole_seq) %in% c('NAIVE', clone_selected_seqs)]
    trimmed_clone_CDR3 <- clone_alignment_CDR3[names(clone_alignment_CDR3) %in% c('NAIVE', clone_selected_seqs)]
    
    write.fasta(trimmed_clone_whole_seq,
                names = names(trimmed_clone_whole_seq),
                file.out = paste0('../results/partis/', dataset, '_HC_clones/',
                                     dataset, '_HC_clone_', clone, '_alignment_TRIMMED.fasta'))
    
    write.fasta(trimmed_clone_CDR3,
                names = names(trimmed_clone_CDR3),
                file.out = paste0('../results/partis/', dataset, '_HC_clones/',
                                  dataset, '_HC_clone_', clone, '_CDR3_alignment_TRIMMED.fasta'))
    
  }
  
}


lapply(as.list(unique(selected_clones$sample_code)),
       FUN = get_trimmed_clones,
       selected_clones = selected_clones)
