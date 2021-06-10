# Plot clone size distributions
library(dplyr)
library(seqinr)
library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(viridis)
library(forcats)

clone_info_files <- list.files('../results/partis/', pattern = 'clone_info.csv')
clone_info_files <- clone_info_files[!grepl('PSC67D_11521_top_clone',clone_info_files)]
clone_info_files <- paste0('../results/partis/', clone_info_files)

read_clone_info <- function(clone_info_file){
  dataset <- str_split(clone_info_file,'\\/')[[1]]
  dataset <- str_replace(dataset[length(dataset)],'_clone_info.csv','')

  #dataset_type <- ifelse(grepl('plate',dataset),'plate','10x')
  
  clone_info_tibble <- as_tibble(read.csv(clone_info_file)) %>%
    mutate(dataset = dataset) %>%
    #mutate(dataset = dataset, dataset_type = dataset_type) %>% 
    mutate(patient = str_replace(str_replace(dataset,'_plate',''),'_HC','')) %>%
    select(dataset, patient, everything())
  return(clone_info_tibble)
}

plot_rank_abundance <- function(clone_info_files, output_dir = '../figures/'){

  clone_info <- lapply(clone_info_files, FUN = read_clone_info)
  clone_info <- bind_rows(clone_info)
  
  size_dist <- clone_info %>% 
    group_by(dataset, patient) %>%
    mutate(rel_frequency = n_seqs / sum(n_seqs)) %>%
    mutate(abundance_rank = rank(-rel_frequency,ties.method = 'random')) %>%
    ungroup()
  
  # Order patients by absolute size of largest clone 
  patient_order <- as.character(size_dist %>% filter(abundance_rank == 1) %>%
    arrange(desc(n_unique_seqs)) %>% pull(patient))
  
  # Patients missing from 10x data will come at the end in whatever order
  #patient_order <- c(patient_order, unique(size_dist$patient)[(unique(size_dist$patient) %in% patient_order) == F])
  
  size_dist$patient = factor(size_dist$patient, levels = patient_order)
  
  # Fraction of sequences represented by largest clone in each data set
  size_dist %>% group_by(dataset) %>%
    summarise(largest_clone_fraction = max(rel_frequency))
  
  # Identity and VDJ use of largest clone in each dataset
  size_dist %>% group_by(dataset) %>%
    filter(n_seqs == max(n_seqs))
  
  # *Number* of sequences represented by largest clone in each data set
  size_dist %>% group_by(dataset) %>%
    summarise(largest_clone_size = max(n_seqs))
  
  cumsum_pl <- size_dist %>%
    group_by(dataset) %>%
    arrange(dataset, abundance_rank) %>%
    mutate(cumulative_sum = cumsum(rel_frequency)) %>%
    ungroup() %>%
    ggplot(aes(x = abundance_rank, y = cumulative_sum, color = patient)) +
    geom_point(alpha = 0.3) +
    xlab('Rank abundance') +
    ylab('Cumulative fraction of sequences')  +
    theme(legend.position = c(0.7,0.2))
  save_plot(paste0(output_dir, 'cumulative_sum.pdf'),
            cumsum_pl,
            base_width = 6)
  
  # Plot with rank abundance in terms of frequency
  rank_abundance_freq_pl <- size_dist %>%
    ggplot(aes(x = abundance_rank, y = rel_frequency, color = patient)) + 
    geom_line() +
    geom_point(alpha = 0.5) +
    xlab('Rank abundance') +
    ylab('Fraction of dataset sequences in clone') +
    scale_x_continuous(breaks = c(1,seq(0,20,5)[-1]),limits = c(1,20)) +
    #scale_y_log10(breaks = c(0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.2),
    #              limits = c(0.005,0.25)) +
    #ylim(-2.5,-1.3) +
    #+ facet_grid(.~dataset_type)
    theme(legend.position = 'none',
          axis.text.y = element_text(size = 8)) 
  
  # Plot with rank abundance in terms of abs numbers
  rank_abundance_abs_pl <- size_dist %>%
    ggplot(aes(x = abundance_rank, y = n_seqs, color = patient)) + 
    geom_line() +
    geom_point(alpha = 0.5) +
    xlab('Rank abundance') +
    ylab('Number of sequences in the clone') +
    scale_x_continuous(breaks = c(1,seq(0,20,5)[-1]),limits = c(1,20)) +
    scale_y_log10(breaks = c(1,5,10,15,30,60, 300)) +
    #ylim(-2.5,-1.3) +
    theme(#legend.position = c(0.83,0.6),
          legend.position = 'left', 
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)) +
    #facet_grid(.~dataset_type) +
    scale_color_discrete(name = 'Patient')
  
  save_plot(paste0(output_dir, 'rank_abundance.pdf'),
            plot_grid(rank_abundance_abs_pl,
                      rank_abundance_freq_pl),
            base_width = 12, base_height = 5)
  
  V_gene_use <- size_dist %>%
    group_by(dataset, v_gene) %>%
    summarise(n_clones = n(), n_seqs = sum(n_seqs)) %>%
    mutate(fraction_clones = n_clones / sum(n_clones),
           fraction_seqs = n_seqs/sum(n_seqs)) %>%
    ungroup()
  stopifnot(all(abs(V_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_clones)) %>% pull(S)) - 1 < 1e-07))
  stopifnot(all(abs(V_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_seqs)) %>% pull(S)) - 1 < 1e-07))
  
  J_gene_use <- size_dist %>%
    group_by(dataset, j_gene) %>%
    summarise(n_clones = n(), n_seqs = sum(n_seqs)) %>%
    mutate(fraction_clones = n_clones / sum(n_clones),
           fraction_seqs = n_seqs/sum(n_seqs)) %>%
    ungroup()
  stopifnot(all(abs(J_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_clones)) %>% pull(S)) - 1 < 1e-07))
  stopifnot(all(abs(J_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_seqs)) %>% pull(S)) - 1 < 1e-07))
  
  D_gene_use <- size_dist %>%
    group_by(dataset, d_gene) %>%
    summarise(n_clones = n(), n_seqs = sum(n_seqs)) %>%
    mutate(fraction_clones = n_clones / sum(n_clones),
           fraction_seqs = n_seqs/sum(n_seqs)) %>%
    ungroup()
  stopifnot(all(abs(D_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_clones)) %>% pull(S)) - 1 < 1e-07))
  stopifnot(all(abs(D_gene_use %>% group_by(dataset) %>% summarise(S = sum(fraction_seqs)) %>% pull(S)) - 1 < 1e-07))
  
  
  # V gene usage in terms of number of clones
  V_gene_clones_pl <- V_gene_use %>% 
    ggplot(aes(x = v_gene, y = fraction_clones, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('V gene') +
    ylab('Fraction of clones in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'V_gene_usage_clones.pdf'),
            V_gene_clones_pl, base_height = 15, base_width = 12)
  
  # V gene usage in terms of number of sequences
  V_gene_seqs_pl <- V_gene_use %>%
    ggplot(aes(x = v_gene, y = fraction_seqs, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('V gene') +
    ylab('Fraction of sequences in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'V_gene_usage_seqs.pdf'),
            V_gene_seqs_pl, base_height = 15, base_width = 12)
  
  # D gene usage in terms of number of clones
  D_gene_clones_pl <- D_gene_use %>% group_by(dataset) %>%
    ungroup() %>%
    ggplot(aes(x = d_gene, y = fraction_clones, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('D gene') +
    ylab('Fraction of clones in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'D_gene_usage_clones.pdf'),
            D_gene_clones_pl, base_height = 15, base_width = 12)
  
  # D gene usage in terms of number of clones
  D_gene_seqs_pl <- D_gene_use %>% group_by(dataset) %>%
    ungroup() %>%
    ggplot(aes(x = d_gene, y = fraction_seqs, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('D gene') +
    ylab('Fraction of sequences in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'D_gene_usage_seqs.pdf'),
            D_gene_seqs_pl, base_height = 15, base_width = 12)
  
  
  # J gene usage in terms of number of clones
  J_gene_clones_pl <- J_gene_use %>% group_by(dataset) %>%
    ungroup() %>%
    ggplot(aes(x = j_gene, y = fraction_clones, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('J gene') +
    ylab('Fraction of clones in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'J_gene_usage_clones.pdf'),
            J_gene_clones_pl, base_height = 10, base_width = 7)
  
  # J gene usage in terms of number of clones
  J_gene_seqs_pl <- J_gene_use %>% group_by(dataset) %>%
    ungroup() %>%
    ggplot(aes(x = j_gene, y = fraction_seqs, fill = dataset)) +
    geom_bar(stat = 'identity') +
    facet_grid(dataset~.) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) + 
    xlab('J gene') +
    ylab('Fraction of sequences in dataset') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'J_gene_usage_seqs.pdf'),
            J_gene_seqs_pl, base_height = 10, base_width = 7)
  
  
  VJ_pair_use <- size_dist %>%
    mutate(v_gene = str_replace(v_gene,'IGHV',''),
           j_gene = str_replace(j_gene,'IGHJ','')) %>%
    group_by(dataset, v_gene, j_gene) %>%
    summarise(n_clones = n(), n_seqs = sum(n_seqs)) %>%
    ungroup() %>% group_by(dataset) %>%
    mutate(fraction_clones = n_clones / sum(n_clones),
           fraction_seqs = n_seqs/sum(n_seqs)) %>%
    ungroup()
  stopifnot(all(abs(VJ_pair_use %>% group_by(dataset) %>% summarise(S = sum(fraction_clones)) %>% pull(S)) - 1 < 1e-07))
  stopifnot(all(abs(VJ_pair_use %>% group_by(dataset) %>% summarise(S = sum(fraction_seqs)) %>% pull(S)) - 1 < 1e-07))
  
  
  VJ_pair_use_clones_pl <- VJ_pair_use %>%
    ggplot(aes(x = v_gene, y = j_gene, fill = fraction_clones)) +
    geom_tile(color = 'black') +
    facet_grid(dataset~.) +
    scale_fill_viridis(name = 'Fraction of clones in dataset') +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
    xlab("V gene") +
    ylab("J gene") +
    theme(legend.position = 'top',
          legend.text = element_text(size = 10)) +
    guides(fill = guide_colourbar(barwidth = 15, barheight = 0.5,
                                  title.position = 'top',
                                  title.hjust = 0.5))
  save_plot(paste0(output_dir, 'VJ_pair_usage_clones.pdf'),
            VJ_pair_use_clones_pl, base_height = 8, base_width = 10)
  
  VJ_pair_use_seqs_pl <- VJ_pair_use %>%
    ggplot(aes(x = v_gene, y = j_gene, fill = fraction_seqs)) +
    geom_tile(color = 'black') +
    facet_grid(dataset~.) +
    scale_fill_viridis(name = 'Fraction of sequences in dataset') +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
    xlab("V gene") +
    ylab("J gene") +
    theme(legend.position = 'top',
          legend.text = element_text(size = 10)) +
    guides(fill = guide_colourbar(barwidth = 15, barheight = 0.5,
                                  title.position = 'top',
                                  title.hjust = 0.5))
  save_plot(paste0(output_dir, 'VJ_pair_usage_seqs.pdf'),
            VJ_pair_use_seqs_pl, base_height = 8, base_width = 10)
  
  # Write summary statistics on clone size to a file
  clone_size_stats <- clone_info %>% group_by(dataset) %>%
    summarise(mean_clone_size = mean(n_seqs),
              clone_size_sd = sd(n_seqs),
              median_clone_size = median(n_seqs),
              max_clone_size = max(n_seqs)) %>%
  write.csv('../results/clone_size_stats.csv', row.names = F)
  
  # Box plot showing clone size distributions
  clone_size_boxplots <- clone_info %>%
    ggplot(aes(x = dataset, y = n_seqs, color = patient)) +
    geom_boxplot() +
    xlab('Dataset') +
    ylab('Clone size') +
    theme(legend.position = 'none')
  save_plot(paste0(output_dir, 'clone_size_boxplot.pdf'),
            clone_size_boxplots, base_height = 5, base_width = 6)
  
  write.csv(size_dist, '../results/clone_freqs.csv', row.names = F)
  
}




plot_rank_abundance(clone_info_files)


