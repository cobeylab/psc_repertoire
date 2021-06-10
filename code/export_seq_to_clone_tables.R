library(tidyverse)
library(readr)
library(yaml)

yaml_files <- list.files('../results/partis', pattern = 'yaml', full.names = T)

link_seqs_to_clones <- function(yaml_object){
  bind_rows(lapply(as.list(yaml_object$events),
         FUN = function(clone){
           seq_ids <- c(clone$unique_ids, unlist(clone$duplicates))
           return(tibble(seq_id = seq_ids))
          
         }), .id = 'clone') %>%
    mutate(clone = as.numeric(clone) - 1)
}



seq_to_clones <- bind_rows(
  lapply(as.list(yaml_files), 
       FUN = function(file_path){
         print(file_path)
         dataset <- str_replace(basename(file_path),'.yaml','')
         yaml_object <- read_yaml(file_path)
         link_seqs_to_clones(yaml_object) %>%
           mutate(dataset = dataset) %>%
           select(dataset, everything())})
  )

write_csv(seq_to_clones, '../results/seqs_in_each_clone.csv')
