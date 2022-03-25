library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())


igphyml_files <- list.files('../results/igphyml/', recursive = T, pattern = 'stats', full.names = T)

divergence_files <- list.files('../results/aa_divergence/', full.names = T)

# Import only divergence results for clones with IgPhyml results 
# (so we only look at the top clones, excluding any that have fewer than 2 seqs.)
igphyml_clone_ids <- unlist(lapply(as.list(igphyml_files),
                            FUN = function(x){
                              clone_id = rev(strsplit(x,'/')[[1]])[1]
                              clone_id = str_remove(clone_id, '_noCDR3_igphyml_stats.tab')
                            }))


divergence_files_clone_ids <- unlist(lapply(as.list(divergence_files),
                                     FUN = function(x){
                                       clone_id = rev(strsplit(x,'/')[[1]])[1]
                                       clone_id = str_remove(clone_id, '_pairwise_divergence.csv')
                                       clone_id = str_remove(clone_id,'_CDR3')
                                     }))

divergence_files <- divergence_files[divergence_files_clone_ids %in% igphyml_clone_ids]
                                     
                                     

igphyml_results <- lapply(as.list(igphyml_files), 
       FUN = function(path){
         dataset <- str_split(path,'\\/')[[1]]
         dataset <- dataset[length(dataset)]
         clone_number <- str_extract(dataset, 'clone_[0-9]+')
         dataset <- str_remove(dataset, '_clone_[0-9]+_.*')

         as_tibble(read.table(path, sep = '\t', header = T)) %>%
           slice(1) %>%
           select(-CLONE) %>%
           mutate(dataset = dataset, clone_number = clone_number) %>%
           select(dataset, clone_number, everything())
       })

igphyml_results <- bind_rows(igphyml_results) %>%
  arrange(dataset)

write_csv(igphyml_results,'../results/combined_igphyml_results.csv')
  
long_format_tibble <- bind_rows(igphyml_results %>% select(dataset, clone_number, OMEGA_FWR_MLE, OMEGA_FWR_LCI,
                                                           OMEGA_FWR_UCI) %>%
                                  dplyr::rename(omega_estimate = OMEGA_FWR_MLE,
                                                omega_lower = OMEGA_FWR_LCI,
                                                omega_upper = OMEGA_FWR_UCI) %>% 
                                  mutate(region = 'FWR') %>% select(dataset, clone_number, region, everything()),
                                igphyml_results %>% select(dataset, clone_number, OMEGA_CDR_MLE, OMEGA_CDR_LCI,
                                                           OMEGA_CDR_UCI) %>%
                                  dplyr::rename(omega_estimate = OMEGA_CDR_MLE,
                                                omega_lower = OMEGA_CDR_LCI,
                                                omega_upper = OMEGA_CDR_UCI) %>% 
                                  mutate(region = 'CDR') %>% select(dataset, clone_number, region, everything()))

dNdS_plot <- long_format_tibble %>%
  ggplot(aes(x = dataset, y = omega_estimate, color = region)) +
  geom_pointrange(size = 1, aes(ymin = omega_lower, ymax = omega_upper), 
                  position = position_jitter(width = 0.2)) +
  ylab('dN/dS estimate from IgPhyML for the biggest clone(s)') +
  xlab('Patient') +
  scale_x_discrete(labels = function(x){str_remove(x,'_HC')}) +
  geom_hline(yintercept = 1, linetype = 2) +
  background_grid() +
  theme(legend.position = c(0.9,0.9)) +
  scale_y_log10()

save_plot('../figures/IgPhyML_dNdS.pdf',dNdS_plot,
          base_width = 12, base_height = 6)
  

divergence_results <- lapply(as.list(divergence_files),
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

# Export divergence results for the same clones we did IgPhyml 
divergence_results <- bind_rows(divergence_results) %>%
  mutate(region = factor(region, levels = c('whole sequence','CDR3')))

#divergence_results <- left_join(igphyml_results %>% select(dataset, clone_number), divergence_results) %>%
#  filter(!is.na(mean_divergence))

write_csv(divergence_results, '../results/divergence_results.csv')

divergence_results_pl <- divergence_results %>%
  filter(reference == 'NAIVE') %>%
  ggplot(aes(x = dataset, y = mean_divergence)) +
  geom_pointrange(size = 1, aes(ymin = min_divergence, ymax = max_divergence), 
                  position = position_jitter(width = 0.2)) +
  ylab('Mean divergence from naive ancestor\n(bar shows min and max)') +
  xlab('Patient') +
  scale_x_discrete(labels = function(x){str_remove(x,'_HC')}) +
  background_grid() +
  theme(legend.position = c(0.9,0.9),
        axis.text.x = element_text(size = 9)) +
  facet_wrap('region')
save_plot('../figures/divergence_from naive.pdf', divergence_results_pl,
          base_width = 20, base_height = 8)
  

  

