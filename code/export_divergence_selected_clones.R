library(dplyr)
library(readr)
library(stringr)

divergence_trimmed_clones_files <- list.files('../results/aa_divergence/', pattern = 'TRIMMED', full.names = T)

divergence_results <- lapply(as.list(divergence_trimmed_clones_files),
                             FUN = function(path){
                               dataset <- str_split(path,'\\/')[[1]]
                               dataset <- dataset[length(dataset)]
                               clone_number <- str_extract(dataset, 'clone_[0-9]+')
                               dataset <- str_remove(dataset, '_clone_[0-9]+_.*')
                               
                               if(grepl('CDR3', path)){
                                 region = 'CDR3'
                               }else{
                                 region = 'whole sequence'
                               }
                               
                               read_csv(path) %>%
                                 mutate(dataset = dataset, clone_number = clone_number, region = region) %>%
                                 select(dataset, clone_number, region, everything())
                               
                             })

divergence_results <- bind_rows(divergence_results)

write_csv(divergence_results,
          '../results/divergence_selected_trimmed_clones.csv')
