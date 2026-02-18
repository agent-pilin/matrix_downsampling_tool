suppressMessages({
  library(yaml)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
})

# parameters <- read_yaml('/home/stanislas.lipin/ownCloud/R&D Scipio/BioInfo/gitlab_scripts/stanislas/Benchmark_3_scipio/mix_samples/matrix_downsampling_tool/config_matrix_downsampling_mix_samples.yml')

### Read config parameters
config_file <- commandArgs(trailingOnly = T)
parameters <- read_yaml(config_file)

count_table_folders <- parameters$count_table_folders
count_prefix <- parameters$count_prefix
downsampling_folders <- parameters$downsampling_folders
downsampling_prefix <- parameters$downsampling_prefix
whitelist_folders <- parameters$whitelist_folders
whitelist_prefix <- parameters$whitelist_prefix
sample_list <- parameters$samples
output_folder <- parameters$output_folder
output_prefix <- parameters$output_prefix
desired_RPIB <- as.integer(parameters$desired_RPIB)
add_RPIB_prefix <- as.logical(as.integer(parameters$add_RPIB_prefix))

### Create palette
getPalette <- colorRampPalette(brewer.pal(9, 'Set1')[-6])
if (length(sample_list) <= 10) {
  my_palette <- c(brewer.pal(8, 'Set1')[c(3, 1, 2, 4, 5, 7, 8)], brewer.pal(8, 'Set2')[c(1, 3, 8)])
  my_palette <- my_palette[1:length(sample_list)]
} else {
  my_palette <- getPalette(length(sample_list))
}
names(my_palette) <- sample_list

### Read data
downsampling_filelist <- outer(downsampling_folders, sample_list, function(x, y) paste0(x, y, downsampling_prefix)) %>% as.vector
downsampling_filelist <- downsampling_filelist[file.exists(downsampling_filelist)]

whitelist_filelist <- outer(whitelist_folders, sample_list, function(x, y) paste0(x, y, whitelist_prefix)) %>% as.vector
whitelist_filelist <- whitelist_filelist[file.exists(whitelist_filelist)]

count_filelist <- outer(count_table_folders, sample_list, function(x, y) paste0(x, y, count_prefix)) %>% as.vector
count_filelist <- count_filelist[file.exists(count_filelist)]

### Load downsampling data 
downsampling_list <- lapply(downsampling_filelist, fread, data.table = F)
whitelist <- lapply(whitelist_filelist, fread, data.table = F)

#### Load count tables, convert to matrices
count_list <- mapply(function(current_filename, current_whitelist) {
  current_table <- fread(current_filename, data.table = F)
  rownames(current_table) <- current_table$gene
  current_table$gene <- NULL
  current_table <- as.matrix(current_table)
  current_table <- current_table[, colnames(current_table) %in% current_whitelist$V1]
  return(current_table)
}, count_filelist, whitelist, SIMPLIFY = F)

### Filter barcodes in downsampling_list by whitelist and calculate median/mean per barcode
downsampling_df <- mapply(function(current_df, current_whitelist) {
  current_df <- current_df %>% filter(barcode %in% current_whitelist$V1) %>%
    dplyr::group_by(ds_coef) %>%
    dplyr::summarize(tr_median_per_barcode = median(transcripts_per_barcode_total),
                     tr_mean_per_barcode = mean(transcripts_per_barcode_total),
                     gene_median_per_barcode = mean(genes_per_barcode_total),
                     reads_per_barcode = mean(reads_per_barcode))
  return(current_df)
}, downsampling_list, whitelist, SIMPLIFY = F)


### Build fitting curve
fit_curve_mean_list <- lapply(downsampling_df, function(x) {
  x <- tryCatch({
    x %>% nls(tr_mean_per_barcode ~ SSlogis(log(reads_per_barcode), Asym, xmid, scal), .)
  },
  error = function(e) {
    x %>% nls(tr_mean_per_barcode ~ SSlogis(reads_per_barcode, Asym, xmid, scal), .)
  })
  return(x)
})

fit_curve_median_list <- lapply(downsampling_df, function(x) {
  x <- tryCatch({
    x %>% nls(tr_median_per_barcode ~ SSlogis(log(reads_per_barcode), Asym, xmid, scal), .)
  },
  error = function(e) {
    x %>% nls(tr_median_per_barcode ~ SSlogis(reads_per_barcode, Asym, xmid, scal), .)
  })
  return(x)
})

fit_curve_gene_median_list <- lapply(downsampling_df, function(x) {
  x <- tryCatch({
    x %>% nls(gene_median_per_barcode ~ SSlogis(log(reads_per_barcode), Asym, xmid, scal), .)
  },
  error = function(e) {
    x %>% nls(gene_median_per_barcode ~ SSlogis(reads_per_barcode, Asym, xmid, scal), .)
  })
  return(x)
})

#### Calculate min reads per barcode among samples (separately by species)
mean_reads_per_barcode <- sapply(downsampling_df, function(x) x$reads_per_barcode[x$ds_coef == 1])

if (desired_RPIB == 0) {
  predict_RPIB <- min(mean_reads_per_barcode)
} else if (desired_RPIB <= min(mean_reads_per_barcode)) {
  predict_RPIB <- desired_RPIB
} else {
  cat(paste0('\nParameter desired_RPIB (', desired_RPIB, ') is higher than minimal RPIB among samples (', round(min(mean_reads_per_barcode), 1), ').\nRPIB per sample:\n\n'))
  print(data.frame(sample = sample_list, RPIB = mean_reads_per_barcode))
  quit()
}

if (add_RPIB_prefix == T) {
  RPIB_prefix <- paste0('_rpib', round(predict_RPIB, 0))
} else {
  RPIB_prefix <- ''
}

transcript_mean_ds <- sapply(fit_curve_mean_list, function(x) unname(predict(x, data.frame(reads_per_barcode = predict_RPIB))))
transcript_median_ds <- sapply(fit_curve_median_list, function(x) unname(predict(x, data.frame(reads_per_barcode = predict_RPIB))))
gene_median_ds <- sapply(fit_curve_gene_median_list, function(x) unname(predict(x, data.frame(reads_per_barcode = predict_RPIB))))

#### Calculate transcript ds coefficient for count tables
matrix_ds_coef <- mapply(function(current_mean_ds, current_table) round(current_mean_ds/mean(colSums(current_table)), 3), transcript_mean_ds, count_list)
matrix_ds_coef <- ifelse(matrix_ds_coef <= 1, matrix_ds_coef, 1)

### Write count tables in files
if (dir.exists(output_folder) == F)
  dir.create(output_folder)

### Build plot
df_plot <- mapply(function(x, y, z, current_sample) {
  df <- data.frame(reads_per_barcode =  seq(1, max(z$reads_per_barcode), 100),
                   transcript_estimation = unname(predict(x, data.frame(reads_per_barcode = seq(1, max(z$reads_per_barcode), 100)))),
                   gene_estimation = unname(predict(y, data.frame(reads_per_barcode = seq(1, max(z$reads_per_barcode), 100)))),
                   sample = current_sample)
  return(df)
}, fit_curve_median_list, fit_curve_gene_median_list, downsampling_df, sample_list, SIMPLIFY = F) %>% rbindlist
df_plot$sample <- factor(df_plot$sample, levels = sample_list)

downsampling_result <- mapply(function(x, current_sample) {
  x$sample <- current_sample
  return(x)
}, downsampling_df, sample_list, SIMPLIFY = F) %>% rbindlist
downsampling_result$sample <- factor(downsampling_result$sample, levels = sample_list)

transcript_df <- data.frame(transcript_median_ds = transcript_median_ds,
                            gene_median_ds = gene_median_ds,
                            sample = sample_list)
ds_plot <- df_plot %>%
  ggplot(aes(reads_per_barcode, transcript_estimation))+
  geom_point(data = downsampling_result, aes(x = reads_per_barcode, y = tr_median_per_barcode, colour = sample, group = ds_coef), alpha = 0.7)+
  geom_line(aes(colour = sample))+
  geom_vline(xintercept = predict_RPIB)+
  geom_hline(data = transcript_df, aes(yintercept = transcript_median_ds, colour = sample), linetype = 'dashed')+
  scale_colour_manual(values = my_palette)+
  labs(x = 'mean RRPC', y = 'ESTIMATED median transcripts per barcode', title = 'Downsampling curves fitted with NLS method')

suppressMessages(ggsave(paste0(output_folder, 'ds_plot', output_prefix, RPIB_prefix, '.png'), ds_plot, device = 'png'))

ds_gene_plot <- df_plot %>%
  ggplot(aes(reads_per_barcode, gene_estimation))+
  geom_point(data = downsampling_result, aes(x = reads_per_barcode, y = gene_median_per_barcode, colour = sample, group = ds_coef), alpha = 0.7)+
  geom_line(aes(colour = sample))+
  geom_vline(xintercept = predict_RPIB)+
  geom_hline(data = transcript_df, aes(yintercept = gene_median_ds, colour = sample), linetype = 'dashed')+
  scale_colour_manual(values = my_palette)+
  labs(x = 'mean RRPC', y = 'ESTIMATED median genes per barcode', title = 'Downsampling curves fitted with NLS method')

suppressMessages(ggsave(paste0(output_folder, 'ds_gene_plot', output_prefix, RPIB_prefix, '.png'), ds_gene_plot, device = 'png'))

#### Downsample count tables (rearrangement)
invisible({
  count_list_ds <- mapply(function(current_table, current_ds_coef, current_sample) {
    cat(current_sample, '\n')
    current_table_ds <- current_table * 0
    current_values <- sample(rep(1:length(current_table), current_table), sum(current_table) * current_ds_coef)
    current_table_ds[as.integer(names(table(current_values)))] <- table(current_values)
    
    index_list <- split(1:nrow(current_table_ds), sort(rep_len(1:100, nrow(current_table_ds))))
    for (i in 1:length(index_list)) {
      cat(current_sample, paste0(i, '%\n'))
      current_df <- current_table_ds[index_list[[i]], ] %>% as.matrix %>% as.data.frame
      current_df <- rownames_to_column(current_df, var = 'gene')
      if (i == 1) {
        fwrite(current_df, paste0(output_folder, current_sample, output_prefix, RPIB_prefix, '_counts.tsv'),
               quote = F, sep = '\t', row.names = F, append = F)
      } else {
        fwrite(current_df, paste0(output_folder, current_sample, output_prefix, RPIB_prefix, '_counts.tsv'),
               quote = F, sep = '\t', row.names = F, append = T)
      }
    }
    return(current_table_ds)
  }, count_list, matrix_ds_coef, sample_list, SIMPLIFY = F)
})

### Summary table
real_transcript_median <- sapply(count_list_ds, function(x) median(colSums(x)))
summary_check <- data.frame(sample = sample_list,
                            estimated_transcript_median = round(transcript_median_ds, 1),
                            real_transcript_median = real_transcript_median)

fwrite(summary_check, paste0(output_folder, 'summary_check', output_prefix, RPIB_prefix, '.tsv'), quote = F, sep = '\t')

### Statistics summary
df_initial_stats <- mapply(function(current_sample, current_table) {
  df <- data.frame(transcript_count = apply(current_table, 2, sum),
                   sample = current_sample)
  return(df)
}, sample_list, count_list, SIMPLIFY = F) %>% rbindlist %>%
  mutate(sample = factor(sample, levels = sample_list))

df_ds_stats <- mapply(function(current_sample, current_table) {
  df <- data.frame(transcript_count = apply(current_table, 2, sum),
                   sample = current_sample)
  return(df)
}, sample_list, count_list_ds, SIMPLIFY = F) %>% rbindlist %>%
  mutate(sample = factor(sample, levels = sample_list))

rbind(df_initial_stats %>% group_by(sample) %>%
  dplyr::summarize(transcript_per_barcode_min = min(transcript_count),
                   transcript_per_barcode_Q1 = quantile(transcript_count)[2],
                   transcript_per_barcode_median = median(transcript_count),
                   transcript_per_barcode_mean = mean(transcript_count),
                   transcript_per_barcode_Q3 = quantile(transcript_count)[4],
                   transcript_per_barcode_max = max(transcript_count),
                   transcript_per_barcode_IQR = IQR(transcript_count),
                   transcript_per_barcode_stdev = sd(transcript_count)) %>%
    mutate(transcript_per_barcode_CV = transcript_per_barcode_stdev/transcript_per_barcode_mean,
           transcript_per_barcode_rIQR = transcript_per_barcode_IQR/transcript_per_barcode_median,
           RPIB = mean_reads_per_barcode,
           transcript_per_barcode_state = 'before downsampling'),
df_ds_stats %>% group_by(sample) %>%
  dplyr::summarize(transcript_per_barcode_min = min(transcript_count),
                   transcript_per_barcode_Q1 = quantile(transcript_count)[2],
                   transcript_per_barcode_median = median(transcript_count),
                   transcript_per_barcode_mean = mean(transcript_count),
                   transcript_per_barcode_Q3 = quantile(transcript_count)[4],
                   transcript_per_barcode_max = max(transcript_count),
                   transcript_per_barcode_IQR = IQR(transcript_count),
                   transcript_per_barcode_stdev = sd(transcript_count)) %>%
  mutate(transcript_per_barcode_CV = transcript_per_barcode_stdev/transcript_per_barcode_mean,
         transcript_per_barcode_rIQR = transcript_per_barcode_IQR/transcript_per_barcode_median,
         RPIB = predict_RPIB,
         transcript_per_barcode_state = 'after downsampling')) %>%
  mutate_if(is.numeric, round, 3) %>%
  fwrite(paste0(output_folder, 'summary_stats', output_prefix, RPIB_prefix, '.tsv'), quote = F, sep = '\t')

print('Mission completed.')
