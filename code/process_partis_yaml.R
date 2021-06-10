library(dplyr)
library(tidyr)
library(yaml)
library(stringr)
library(msa)
library(seqinr)
library(Biostrings)
source('imgt_region_bounds.R')

args = commandArgs(trailingOnly = T)
yaml_file_path = args[1] # yaml_file_path = '../results/partis/PSC50B_HC.yaml'
igblast_annotation_path = args[2] # Path to igblast annotation, e.g. igblast_annotation_path = "../results/partis/PSC50B_HC_igblast_db-pass.tab"
isotype_info_path = args[3] # isotype_info_path = '../data/PSC50B_HC_isotypes.csv'

if(is.na(isotype_info_path)){
  isotype_info_path <- NULL
}

# Converts a range of positions in an ungapped sequence to a range in the same sequence with gaps
# For CDR3 bounds in alignment given bounds in germline sequence

convert_position_range <- function(ungapped_seq, gapped_seq, ungapped_range){
  ungapped_seq = strsplit(ungapped_seq, '')[[1]]
  gapped_seq = strsplit(gapped_seq, '')[[1]]
  adjusted_range = ungapped_range
  
  n_gaps = 0
  
  for(i in 1:length(gapped_seq)){
    if((gapped_seq[i]) %in% c('-','.')){
      n_gaps <- n_gaps + 1
      position_in_ungapped_seq = i - n_gaps
      # If gap occurs before the start of the range in the ungapped seq, entire range is pushed by 1
      if(position_in_ungapped_seq <  ungapped_range[1]){
        adjusted_range <- adjusted_range + 1
      }
      # If gap is within the range in the ungapped seq, end of range is pushed by 1
      if(position_in_ungapped_seq >=  ungapped_range[1] & position_in_ungapped_seq < ungapped_range[2]){
        adjusted_range[2] <- adjusted_range[2] + 1
      }
    }
  }
  return(adjusted_range)
} 

# Some tests for the conversion function, including cases with gaps at either end of range
stopifnot(all(convert_position_range('AACTAAAG','A-ACTAA-AG',c(5,7)) == c(6,9)))
stopifnot(all(convert_position_range('AACTAAAG','AACTAAAG---',c(5,7)) == c(5,7)))
stopifnot(all(convert_position_range('AACTAAAG','AACTA-AAG',c(5,7)) == c(5,8)))
stopifnot(all(convert_position_range('AACTAAAG','AACT-AAAG',c(5,7)) == c(6,8)))
stopifnot(all(convert_position_range('AACTAAAG','AACTAAA-G',c(5,7)) == c(5,7)))
stopifnot(all(convert_position_range('AACTAAAG','A-ACTAAAG',c(5,7)) == c(6,8)))

# Pulls sequence of a FR or CDR from a matrix-object alignment
pull_region_seq <- function(aln, imgt_region_bounds, region){
  start <- imgt_region_bounds %>% filter(region == !!region) %>% pull(start)
  end <- imgt_region_bounds %>% filter(region == !!region) %>% pull(end)
  
  if(nrow(aln) > 1){
    seqs <- apply(aln[, start:end], 1, FUN = paste, collapse = '')
  }else{
    seqs <- paste(aln[, start:end], collapse = '')
  }
  return(seqs)
}

# Given a character vector of aligned sequences, removes sites that are all gaps
# Can be used to remove sites that are all 'N' if gap_character set to 'N'
remove_allgap_sites <- function(seqs, gap_character = '-'){
  if(length(seqs) >1){
    seqmatrix <- c()
    for(s in seqs){
      seqmatrix <- rbind(seqmatrix, str_split(s,'')[[1]])
    }
    rownames(seqmatrix) <- names(seqs)
    allgap_sites <- apply(seqmatrix, 2, FUN = function(x){return(all(x == gap_character))})
    seqmatrix <- seqmatrix[,allgap_sites == F]
    adjusted_seqs <- apply(seqmatrix,1,paste,collapse ='')
  }else{
    adjusted_seqs <- seqs
  }
  return(adjusted_seqs)
}
stopifnot(remove_allgap_sites(c('--AAA','---BB')) == c("AAA", "-BB"))
stopifnot(remove_allgap_sites(c('NNAAA','NNNBB'), gap_character = 'N') == c("AAA", "NBB"))

# Gets number of uniq seqs, v/d/j, n uniqu igg seqs. Also returns vector of seqs (incl. germline) with Ns removed
get_basic_clone_info <- function(clone, isotype_info = NULL){
  clone_info <- tibble(n_seqs = length(clone$unique_ids) + 
                         sum(sapply(clone$duplicates, FUN = length)),
                       n_unique_seqs =  length(clone$unique_ids),
                       v_gene = clone$v_gene, d_gene = clone$d_gene, j_gene = clone$j_gene)
  
   clone_seqs <- tibble(seq_id = clone$unique_ids, in_frame = clone$in_frames, has_stop_codon = clone$stops,
                       key_residues_conserved = !clone$mutated_invariants, seq = clone$input_seqs)
  if(!is.null(isotype_info)){
    clone_seqs <- left_join(clone_seqs, isotype_info %>% mutate(seq_id = Name, isotype = Isotype) %>% 
                              select(seq_id, isotype))
    n_unique_igg <- length(clone_seqs %>% filter(grepl('IgG', isotype)) %>% pull(seq_id))
    n_unique_iga <- length(clone_seqs %>% filter(grepl('IgA', isotype)) %>% pull(seq_id))
    
    clone_info <- clone_info %>% mutate(n_unique_igg, n_unique_iga) %>%
      select(n_seqs, n_unique_seqs, n_unique_igg, n_unique_iga, everything())
  }
  
  # Add germline sequence and remove 'N' characters introduced by partis for sites outside of VDJ region
  clone_seqs <- clone_seqs %>%
    mutate(seq = str_replace_all(seq,'N','')) %>% pull(seq)
  
  clone_seqs <- c(str_replace_all(clone$naive_seq, 'N',''), clone_seqs)
  names(clone_seqs) <- c('NAIVE', clone$unique_ids)
    
  return(list('info' = clone_info, 'seqs' = clone_seqs))
}

# Replaces codons containing a gap with 'NNN'
mask_gapped_codons <- function(seqs){
  mask_function <- function(s){
    masked_s <- ''
    n_codon_sites <- ceiling(nchar(s)/3)
    codon_start_positions <- c(1, (2:n_codon_sites)*3 - 2)
    for(codon_start in codon_start_positions){
      # If sequence ends before codon is complete
      if(codon_start + 2 > nchar(s)){
        codon <- 'NNN'
      }else{
        codon <- substr(s, codon_start , codon_start +2)
        if(grepl('-',codon)){
          codon <- 'NNN'
        }
      }
      masked_s <- paste0(masked_s, codon)
    }
    return(masked_s)
  }
  # Apply internal function to all sequences
  masked_seqs <- sapply(as.list(seqs), FUN = mask_function)
  
  # If a sequence is shorter than the others, fill difference with masked codons
  max_seq_length = max(nchar(masked_seqs))
  for(i in 1:length(masked_seqs)){
    if(nchar(masked_seqs[i]) < max_seq_length){
      stopifnot((max_seq_length - nchar(masked_seqs[i])) %% 3 == 0)
      masked_seqs[i] <- paste(c(masked_seqs[i], rep('N', max_seq_length - nchar(masked_seqs[i]))), collapse = '')
    }
  }
  names(masked_seqs) <- names(seqs)
  return(masked_seqs)
}

# Given a character vector of aligned sequences, returns whether each has a stop codon
check_stop_codons <- function(seqs){
  # Replace gaps with 'Ns' for translation
  seqs <- str_replace_all(seqs, '-','N')
  
  has_stop_codons <- as.character(Biostrings::translate(DNAStringSet(seqs),
                                                        if.fuzzy.codon = 'solve'))
  has_stop_codons <- grepl('\\*', has_stop_codons)
  return(has_stop_codons)
}
stopifnot(check_stop_codons(c('AAAGGG','TAGAAG','TT----')) == c(F,T,F))

# Takes partis yaml and exports clone-specific fasta files and igphyml input
process_partis_yaml <- function(yaml_file_path, igblast_annotation_path, output_dir, isotype_info_path = NULL){
  yaml_object <- read_yaml(yaml_file_path)
  igblast_annotation <- as_tibble(read.table(igblast_annotation_path, sep = '\t', header = T))
  dataset_name <- str_replace(basename(yaml_file_path),'.yaml','')
  
  master_input_file <-  paste0(output_dir, 'igphyml_input_', dataset_name,'.tsv')
  master_input_file_noCDR3 <-  paste0(output_dir,'igphyml_input_', dataset_name,'_noCDR3.tsv')
  # Excluded clones (ones in which inferred FR/CDR region has a length that's not a multiple of 3)
  excluded_clones_path <-  paste0(output_dir, dataset_name,'_clones_excluded_from_igphyml.csv')
  clone_info_path <- paste0(output_dir, dataset_name,'_clone_info.csv')
  
  master_input <- c()
  clone_number <- 0
  excluded_clones <- c()
  clone_info <- c() # Tibble with summary information for each clone
  
  # Read isotype information, if available
  if(!is.null(isotype_info_path)){
    isotype_info <- as_tibble(read.csv(isotype_info_path))
  }else{
    isotype_info <- NULL
  }
  
  parent_clone_dir <- paste0(output_dir, dataset_name, '_clones/')
  dir.create(parent_clone_dir, showWarnings = F)
  
  parent_igphyml_dir <- paste0(output_dir, dataset_name, '_igphyml/')
  dir.create(parent_igphyml_dir, showWarnings = F)
  
  for(clone in yaml_object$events){
    print(paste('Processing clone',clone_number +1))
    # Path to igphyml input files for this clone
    igphyml_seq_path <- paste0(parent_igphyml_dir, dataset_name,'_clone_', clone_number,'_igphyml.fasta')
    igphyml_seq_path_noCDR3 <- paste0(parent_igphyml_dir, dataset_name,'_clone_', clone_number,'_noCDR3_igphyml.fasta')
    igphyml_partition_path <-  paste0(parent_igphyml_dir, dataset_name,'_clone_', clone_number,'_part.txt')
    igphyml_partition_path_noCDR3 <-  paste0(parent_igphyml_dir, dataset_name,'_clone_',
                                             clone_number,'_noCDR3_part.txt')
    
    # Clone fasta file (with 'N' characters removed but not processed for igphyml):
    raw_clone_fasta_file <- paste0(parent_clone_dir, dataset_name,'_clone_',
                                   clone_number,'.fasta') 
    
    # Update tibble with general clone information
    clone_info <- bind_rows(clone_info,
                            get_basic_clone_info(clone, isotype_info)$info %>%
                              mutate(clone_id = clone_number)) %>%
      select(clone_id, everything())
    
    # Write raw clone seqs + germline (with 'N's removed but without processing for igphyml)
    clone_seqs <- get_basic_clone_info(clone, isotype_info)$seqs
    write.fasta(as.list(clone_seqs), names = names(clone_seqs), file.out =  raw_clone_fasta_file)
    
    # Constrain igphyml input to productive sequences
    is_productive <- (clone$stops == F)&(clone$in_frames)&(clone$mutated_invariants == F)
    
    productive_seqs <- clone$unique_ids[is_productive]
    
    if(length(productive_seqs)>0){
      
      # Check all clone seqs. in yaml file are also in presto annotation file. 
      stopifnot(all(clone$unique_ids %in% igblast_annotation$SEQUENCE_ID))
      
      # Initial processing of inferred naive sequence (remove 'N' characters at each end of germline sequence)
      naive_seq <- clone$naive_seq
      
      germline_Ns <- str_locate_all(naive_seq,'N+')[[1]]
      n_start_Ns <- 0
      terminal_Ns_start = nchar(naive_seq) +1
      if(nrow(germline_Ns) > 0){
        germline_start_Ns <- germline_Ns[1,]
        if(germline_start_Ns[1] == 1){
          n_start_Ns <- germline_start_Ns[2] 
        }
        if(nrow(germline_Ns) > 1){
          terminal_N_sites = germline_Ns[nrow(germline_Ns),]
          if(terminal_N_sites[2] == nchar(naive_seq))
            terminal_Ns_start = terminal_N_sites[2]
        }
      }
      
      # Find naive sequence and IMGT-aligned observed sequences from this clone
      naive_seq <- substr(naive_seq, n_start_Ns +1, terminal_Ns_start -1)
      imgt_aligned_seqs <- left_join(tibble(SEQUENCE_ID = productive_seqs),
                                     igblast_annotation %>% filter(SEQUENCE_ID %in% productive_seqs),
                                     by = 'SEQUENCE_ID') 
      # Check all imgt-aligned regions have the same length (except possibly CDR3 and FR4)
      lengths_match <- imgt_aligned_seqs %>% summarise_at(.vars =  vars(c("FWR1_IMGT","FWR2_IMGT", "FWR3_IMGT",
                                                                          "CDR1_IMGT","CDR2_IMGT")),
                                                          .funs = function(x){return(length(unique(nchar(as.character(x)))))}) %>%
        unlist() == 1
      stopifnot(all(lengths_match))
      
      # Final sequence tibble
      final_seqs <- imgt_aligned_seqs %>% select(SEQUENCE_ID,FWR1_IMGT, CDR1_IMGT, FWR2_IMGT,
                                                 CDR2_IMGT, FWR3_IMGT, CDR3_IMGT, FWR4_IMGT) %>%
        mutate_at(.vars = vars(matches('IMGT')), .funs = str_replace_all, pattern = '\\.', replacement = '-') %>%
        # Mask gapped codons in each region
        mutate_at(.vars = vars(matches('IMGT')), .funs = mask_gapped_codons) %>%
        # Remove sites that are 'N' in all sequences
        mutate_at(.vars = vars(matches('IMGT')), .funs = remove_allgap_sites, gap_character = 'N')
      
      # Get final lengths of each region (after removing all-gap sites)
      region_lengths <- final_seqs %>% summarise_at(vars(matches('IMGT')),
                                                    .funs = function(x){return(unique(nchar(x)))}) %>%
        unlist()
      names(region_lengths) <- str_replace(names(region_lengths),'_IMGT','')
      
      # Full sequences (including CDR3)
      full_seqs <- final_seqs %>%
        rowwise() %>%
        mutate(full_seq = paste(c(FWR1_IMGT, CDR1_IMGT, FWR2_IMGT,
                                  CDR2_IMGT, FWR3_IMGT, CDR3_IMGT, FWR4_IMGT), collapse = '')) %>%
        ungroup() %>%
        select(SEQUENCE_ID, full_seq)
      
      # Align germline against IMGT-aligned observed sequences to identify its CDR3 boundaries in that alignment
      alignment_input <- c(naive_seq, as.character(full_seqs$full_seq))
      names(alignment_input) <- c('NAIVE', full_seqs$SEQUENCE_ID)
      
      aln <- msa(alignment_input, type = 'dna')
      aln <- as.matrix(aln)
      aln <- aln[c('NAIVE', full_seqs$SEQUENCE_ID), ]
      
      # Get CDR3 positions for un-alingned germline sequence (provided by partis)
      # Positions of conserved cysteine in V and conserved Try/Phe in J, bounds of CDR3
      # Positions are in NT, but correspond to first NT in the respective codon (convert to R numbering)
      germline_CDR3_positions = c(clone$codon_positions$v+ 1, clone$codon_positions$j + 3)
      
      # From the germline CDR3 positions, get corresponding positions in the alignment (accounting for N removal)
      CDR3_aln_bounds <- convert_position_range(naive_seq, paste0(aln['NAIVE',], collapse = ''),
                                                ungapped_range = germline_CDR3_positions - n_start_Ns)
      
      # Conserved cysteine is considered part of FR3 in IMGT, conserved tryp is not. Correct for that
      CDR3_aln_bounds <- CDR3_aln_bounds + c(3,-3)
      
      # Get aligned CDR3 sequences (including aligned NAIVE)
      aln_CDR3_seqs <- pull_region_seq(aln, tibble(region = 'CDR3', start = CDR3_aln_bounds[1],
                                                   end = CDR3_aln_bounds[2]), 'CDR3')
      
      # Check assignment of germline CDR3 by checking aln positions identify correct CDR3 of observed sequences.
      aln_consistent = all(aln_CDR3_seqs[-1] == str_replace_all(as.character(imgt_aligned_seqs %>% pull(CDR3_IMGT)), '\\.','-'))
      
      # Get aligned naive sequence, and its versions with NP1-D-NP2 masked and without CDR3
      aligned_naive_seq = paste(aln['NAIVE',], collapse = '')
      aligned_naive_seq_noCDR3 <- aln['NAIVE',]
      aligned_naive_seq_noCDR3 <- paste(aligned_naive_seq_noCDR3[-seq(CDR3_aln_bounds[1],CDR3_aln_bounds[2])],
                                        collapse = '')

      # Masking NP1-D-NP2 region
      germline_v_bounds <- clone$regional_bounds$v + c(1,0) # Converts from python slicing convention to R numbering
      germline_j_bounds <- clone$regional_bounds$j + c(1,0)
      germline_v_bounds <- germline_v_bounds - n_start_Ns
      germline_j_bounds <- germline_j_bounds - n_start_Ns
      
      # Adjusts NP1-D-NP2 bounds for alignment
      aln_v_bounds <- convert_position_range(naive_seq, aligned_naive_seq, germline_v_bounds)
      aln_j_bounds <- convert_position_range(naive_seq, aligned_naive_seq, germline_j_bounds)
      
      dmasked_germline_seq = str_split(aligned_naive_seq,'')[[1]]
      dmasked_germline_seq[(aln_v_bounds[2] + 1):(aln_j_bounds[1] - 1)] <- 'N'
      dmasked_germline_seq <- paste(dmasked_germline_seq, collapse = '')
      
      # For full sequence tibble, use naive sequence with masked NP1-D-NP2
      full_seqs <- bind_rows(tibble(SEQUENCE_ID = 'NAIVE', full_seq = dmasked_germline_seq),
                             full_seqs) %>%
        # Trim any excess terminal nucleotides introduced by germline sequence
        mutate(full_seq = mask_gapped_codons(full_seq)) %>%
        mutate(full_seq = remove_allgap_sites(full_seq, gap_character = 'N')) %>%
        # Check for stop codons
        mutate(has_stop_codons = check_stop_codons(full_seq)) %>%
        filter(has_stop_codons == F)
      
      # If at least one sequence (plus the naive sequence) has been retained, export clone igphyml input
      # And if alignment of germline sequence is consistent
      if(nrow(full_seqs) >= 2 & aln_consistent){
        # If a single observed sequence is present, igphyml requires it to be duplicated
        if(sum(full_seqs$SEQUENCE_ID != 'NAIVE') == 1){
          full_seqs <- bind_rows(full_seqs,
                                 full_seqs %>% filter(SEQUENCE_ID != 'NAIVE') %>% 
                                   mutate(SEQUENCE_ID = paste0(SEQUENCE_ID,'_duplicate')))
        }
        
        write.fasta(as.list(full_seqs$full_seq), full_seqs$SEQUENCE_ID, igphyml_seq_path)
        writeLines(c(paste('2', sum(region_lengths)/3),
                     'FWR:IMGT',
                     'CDR:IMGT',
                     clone$v_gene,
                     clone$j_gene,
                     paste(c(rep('13', region_lengths['FWR1']/3), rep('30', region_lengths['CDR1']/3),
                             rep('45',  region_lengths['FWR2']/3), rep('60', region_lengths['CDR2']/3),
                             rep('80', region_lengths['FWR3']/3), rep('108', region_lengths['CDR3']/3),
                             rep('120', region_lengths['FWR4']/3)), collapse = ',')
        ),
        con = file(igphyml_partition_path))
        
        # Input fasta and partition files excluding CDR3
        noCDR3_seqs <- final_seqs %>%
          rowwise() %>%
          mutate(noCDR3_seq = paste(c(FWR1_IMGT, CDR1_IMGT, FWR2_IMGT, CDR2_IMGT, FWR3_IMGT,
                                      FWR4_IMGT), collapse = '')) %>%
          ungroup() %>%
          select(SEQUENCE_ID, noCDR3_seq)
        
        # Add germline sequence with CDR3 removed
        noCDR3_seqs <- bind_rows(tibble(SEQUENCE_ID = 'NAIVE', noCDR3_seq = aligned_naive_seq_noCDR3),
                                 noCDR3_seqs) %>%
          # Trim any excess terminal nucleotides introduced by germline sequence
          mutate(noCDR3_seq = mask_gapped_codons(noCDR3_seq)) %>%
          mutate(noCDR3_seq = remove_allgap_sites(noCDR3_seq, gap_character = 'N')) %>%
          # Check for stop codons
          mutate(has_stop_codons = check_stop_codons(noCDR3_seq)) %>%
          filter(has_stop_codons == F)
        
        # Again if a single observed sequence is present, igphyml requires it to be duplicated
        if(sum(noCDR3_seqs$SEQUENCE_ID != 'NAIVE') == 1){
          noCDR3_seqs <- bind_rows(noCDR3_seqs,
                                   noCDR3_seqs %>% filter(SEQUENCE_ID != 'NAIVE') %>% 
                                     mutate(SEQUENCE_ID = paste0(SEQUENCE_ID,'_duplicate')))
        }
        
        write.fasta(as.list(noCDR3_seqs$noCDR3_seq), noCDR3_seqs$SEQUENCE_ID, igphyml_seq_path_noCDR3)
        writeLines(c(paste('2',sum(region_lengths[names(region_lengths)!='CDR3'])/3),
                     'FWR:IMGT',
                     'CDR:IMGT',
                     clone$v_gene,
                     clone$j_gene,
                     paste(c(rep('13', region_lengths['FWR1']/3), rep('30', region_lengths['CDR1']/3),
                             rep('45', region_lengths['FWR2']/3), rep('60', region_lengths['CDR2']/3),
                             rep('80', region_lengths['FWR3']/3), rep('120', region_lengths['FWR4']/3)),
                           collapse = ',')
        ),
        con = file(igphyml_partition_path_noCDR3))
        
        master_input <- bind_rows(master_input,
                                  tibble(cdr3_excluded = c(F,T),
                                         fasta_file = c(igphyml_seq_path, igphyml_seq_path_noCDR3), start_tree = 'N',
                                         partition_file = c(igphyml_partition_path, igphyml_partition_path_noCDR3)))
        
      }else{
        # All productive sequences have stop codons
        excluded_clones <- bind_rows(excluded_clones,
                                     tibble(clone = paste0(dataset_name, '_clone_', clone_number),
                                            n_seqs = length(!!clone$unique_ids)))
      }
    }else{
      # clone nas no productive sequences
      excluded_clones <- bind_rows(excluded_clones,
                                   tibble(clone = paste0(dataset_name, '_clone_', clone_number),
                                          n_seqs = length(!!clone$unique_ids)))
    }
   clone_number <- clone_number + 1
  
  }
  # Write master igphyml input files
  writeLines(c(sum(master_input$cdr3_excluded == F),
               master_input %>% filter(cdr3_excluded == F) %>%
                 mutate(line = paste(fasta_file, start_tree, 'NAIVE', partition_file, sep = '\t')) %>% 
                 pull(line)),
             con = file(master_input_file))
  
  writeLines(c(sum(master_input$cdr3_excluded == 1),
               master_input %>% filter(cdr3_excluded == T) %>%
                 mutate(line = paste(fasta_file, start_tree, 'NAIVE', partition_file, sep = '\t')) %>% 
                 pull(line)),
             con = file(master_input_file_noCDR3))
  
  # Export clone info
  write.csv(clone_info, clone_info_path, row.names = F)
  write.csv(excluded_clones, excluded_clones_path, row.names = F)
  
}

process_partis_yaml(yaml_file_path, igblast_annotation_path, output_dir = '../results/partis/', isotype_info_path)
