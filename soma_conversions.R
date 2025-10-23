## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------
library(rjson)
library(biomaRt)
library(clusterProfiler)
library(DESeq2)
library(dplyr)
library(org.Mm.eg.db)
library(tidyr)
library(pheatmap)
library(ggplot2)



## ---------------------------------------------------------------------------------------------------------------
input <- "C:/Users/prv892/Documents/GitHub/slam/metadata.json"

file_path_ref <- fromJSON(file = input)
#source(file.path(file_path_ref$project_directory, file_path_ref$sequencing$setup_script))
#source(file.path(file_path_ref$project_directory, file_path_ref$plots$setup_script))

metadata <- read.csv(file.path(file_path_ref$project_directory, file_path_ref$design))

background_colors <- colorRampPalette(c("white", "black"))(10)
bkg_sample_colors <- background_colors[c(3,5,7)]

accent_colors <- c("#a50026",
"#d73027",
"#f46d43",
"#fdae61",
"#fee090",
"#e0f3f8",
"#abd9e9",
"#74add1",
"#4575b4",
"#313695")

qc_peaks <- list(actb="peak_12290")

# create ensembl biomart instance 
grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")


## ---------------------------------------------------------------------------------------------------------------
read_counts_table <- function(counts_file) {
  
  res <- read.table(counts_file, sep="\t", 
                    skip=1,
                    col.names=c("Geneid", "TC_count_total", "T_count_total", 
                                "Conversions_total", "Coverage_total",
                                sprintf("reads_%s_TC_conversions", 0:10),
                                sprintf("reads_%s_any_conversions", 0:10),
                                "invalid_reads"))
  res$tc_geq1_reads_total <- rowSums(res[,sprintf("reads_%s_TC_conversions", 1:10)])
  res$tc_geq2_reads_total <- rowSums(res[,sprintf("reads_%s_TC_conversions", 2:10)])

  res$reads_total <- rowSums(res[,sprintf("reads_%s_TC_conversions", 0:10)])

  res$sample_id <- gsub(".conversions.txt", "", basename(counts_file))
  
  return(res)
}

get_soma_counts_data <- function() {
  
  meta_soma <- subset(metadata, location == "Soma")
  
  counts_files <- file.path(file_path_ref$project_directory, 
                            file_path_ref$results$conversion_counts_dir,
                            sprintf("%s.conversions.txt", meta_soma$sample_id))
  
  counts_soma <- data.frame(do.call(rbind, 
                                    lapply(counts_files, read_counts_table))) %>%
    inner_join(meta_soma, by="sample_id")
  

  save(counts_soma, file=file.path(file_path_ref$project_directory, file_path_ref$results$soma_counts_data))
  
  
}

get_soma_counts_data()

get_axon_counts_data <- function() {
  
  meta_axon <- subset(metadata, grepl("Axon", location))
  
  counts_files <- file.path(file_path_ref$project_directory, 
                            file_path_ref$results$conversion_counts_dir,
                            sprintf("%s.conversions.txt", meta_axon$sample_id))
  
  counts_axon <- data.frame(do.call(rbind, 
                                    lapply(counts_files, read_counts_table))) %>%
    inner_join(meta_axon, by="sample_id")
  

  save(counts_axon, file=file.path(file_path_ref$project_directory, file_path_ref$results$axon_counts_data))
  
}

get_axon_counts_data()

get_gc_counts_data <- function() {
  
  meta_gc <- subset(metadata, location == "GC")
  
  counts_files <- file.path(file_path_ref$project_directory, 
                            file_path_ref$results$conversion_counts_dir,
                            sprintf("%s.conversions.txt", meta_gc$sample_id))
  
  counts_gc <- data.frame(do.call(rbind, 
                                    lapply(counts_files, read_counts_table))) %>%
    inner_join(meta_gc, by="sample_id")
  

  save(counts_gc, file=file.path(file_path_ref$project_directory, file_path_ref$results$gc_counts_data))
  
}

get_gc_counts_data()


## ---------------------------------------------------------------------------------------------------------------
compare_error_rates_samples <- function() {
  
   meta_soma <- subset(metadata, location == "Soma" & IUE_plasmid == "pPV25") 

  valid_genes_soma <- (counts_soma %>% 
                         subset(sample_id %in% meta_soma$sample_id) %>%
    group_by(Geneid) %>%
    summarize(meets_thresh=sum(reads_total > 100)/length(reads_total) == 1 & sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(meets_thresh))$Geneid
  
  res_long <- counts_soma %>%
                      subset(Geneid %in% valid_genes_soma & sample_id %in% meta_soma$sample_id) %>%
    rowwise() %>%
    mutate(error_rate=(Conversions_total-TC_count_total)/(Coverage_total))
  
  # plot the median log error rates for each gene
  plt <- res_long %>% 
    group_by(Geneid) %>% 
    summarize(median_error_rate=median(error_rate)) %>%
    ggplot(aes(median_error_rate)) + 
    geom_histogram(bins=100) +
    theme_classic() + 
    scale_x_log10()
  
  plot(plt)
  
  res <- data.frame(res_long %>%
                      rowwise() %>%
                      mutate(log_error_rate=log(((Conversions_total-TC_count_total)+1)/(Coverage_total))) %>%
      dplyr::select(Geneid, sample_id, log_error_rate) %>%
      pivot_wider(id_cols=Geneid,
                  names_from=sample_id,
                  values_from=log_error_rate))
  
  rownames(res) <- res$Geneid
  res$Geneid <- NULL
  res <- na.omit(res)
  error_rate_cor <- cor(res)
  
  annotation_row <- data.frame(sample_time_min=meta_soma[,c("sample_time_min")],
                               row.names=meta_soma$sample_id)
  annotation_row$sample_time_min <- sprintf("%s_min", annotation_row$sample_time_min)
  annotation_colors <- accent_colors[c(10, 8, 6, 5, 3, 1)]
  names(annotation_colors) <- sprintf("%s_min", c(0, 120, 240, 360, 480, 720))
 
  pheatmap(error_rate_cor, annotation=annotation_row)
  
  
}

compare_error_rates_samples()


## ---------------------------------------------------------------------------------------------------------------

get_t_counts_mat <- function(counts_df) {
  res <- data.frame(counts_df %>% 
                             dplyr::select(Geneid, sample_id, T_count_total) %>%
                             pivot_wider(id_cols=Geneid,
                                         names_from=sample_id,
                                         values_from=T_count_total))
  
  rownames(res) <- res$Geneid
  res$Geneid <- NULL
  
  return(res)
  
}

get_tc_counts_mat <- function(counts_df) {
  res <- data.frame(counts_df %>% 
                             dplyr::select(Geneid, sample_id, TC_count_total) %>%
                             pivot_wider(id_cols=Geneid,
                                         names_from=sample_id,
                                         values_from=TC_count_total))
  
  rownames(res) <- res$Geneid
  res$Geneid <- NULL
  
  return(res)
  
}


get_tc_conversions_mat <- function(counts_df) {
  res <- data.frame(counts_df %>% 
                      rowwise() %>%
                      mutate(tc_conversions=TC_count_total/T_count_total) %>%
                             dplyr::select(Geneid, sample_id, tc_conversions) %>%
                             pivot_wider(id_cols=Geneid,
                                         names_from=sample_id,
                                         values_from=tc_conversions,
                                         values_fill = list(tc_conversions=0)))
  res[is.na(res)] <- 0
  
  rownames(res) <- res$Geneid
  res$Geneid <- NULL
  return(res)
}

get_rangenorm_conversions_mat <- function(counts_df, zero_cols) {
  tc_conversions <- get_tc_conversions_mat(counts_df)

  row_max <- apply(tc_conversions, 1, max)
  zero_convrate <- rowMeans(tc_conversions[,zero_cols])
  
  tc_conversions_minus_zero <- apply(sweep(tc_conversions, 1, zero_convrate, "-"),
                                     c(1, 2), function(x) max(0, x))
  
  tc_conversions_pct_max <- sweep(tc_conversions_minus_zero, 1, row_max, "/")
  
  tc_conversions_pct_max[is.na(tc_conversions_pct_max)] <- 0
  
  return(tc_conversions_pct_max)
}

get_normalized_conversions <- function() {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$soma_counts_data))
  
  meta_soma <- subset(metadata, location == "Soma" & IUE_plasmid == "pPV25") # sample_time_min < 720)
  meta_soma$log_sample_time_h <- sapply(meta_soma$sample_time_min, function(t_min) ifelse(t_min == 0, 0, log2(t_min/60)))

  valid_genes_soma <- (counts_soma %>% 
    group_by(Geneid) %>%
    summarize(meets_thresh=sum(reads_total > 100)/length(reads_total) > 2/3 & sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(meets_thresh))$Geneid
  
  counts_soma_filt <- counts_soma %>% subset(Geneid %in% valid_genes_soma & sample_id %in% meta_soma$sample_id)
  
  tc_conversion_df <- get_tc_conversions_mat(counts_soma_filt)
  
  bkg_samples <- subset(meta_soma, sample_time_min == 0)$sample_id
  tc_conversion_bkg <- tc_conversion_df[,bkg_samples]
  tc_conversion_bkg_avg <- rowMeans(tc_conversion_bkg)
  
  signal_samples <- subset(meta_soma, sample_time_min > 0)$sample_id
  tc_conversion_signal <- tc_conversion_df[,signal_samples]
  
  tc_conversion_df_minus_bkg <- sweep(tc_conversion_signal, 1, tc_conversion_bkg_avg, "-")
  
  tc_conversion_df_minus_bkg_long <- tc_conversion_df_minus_bkg
  tc_conversion_df_minus_bkg_long$Geneid <- rownames(tc_conversion_df_minus_bkg_long)
  tc_conversion_df_minus_bkg_long <- tc_conversion_df_minus_bkg_long %>%
    pivot_longer(contains("SLAM"), names_to="sample_id", values_to="conversion_rate") %>%
    inner_join(meta_soma, by="sample_id")
  
 
  
  return( tc_conversion_df_minus_bkg_long)
}


get_utrs <- function(peak_ids) {
  load(file.path(file_path_ref$project_directory, 
                 file_path_ref$results$quantseq_3end_gene_annotated))
  
  ensembl_transcript_ids <- subset(quantseq_utr_annotated_filt, 
                                  peak_name %in% peak_ids)$utr_name
  
  utr3_seqs <- biomaRt::getBM(c("ensembl_transcript_id", "3utr"),
                              filters="ensembl_transcript_id",
                              values=ensembl_transcript_ids,
                              mart=grcm38_101) %>%
    inner_join(quantseq_utr_annotated_filt, by=c("ensembl_transcript_id" = "utr_name")) %>%
    subset(!grepl("N|S", `3utr`) & nchar(`3utr` > 10)) %>%
    rowwise() %>%
    mutate(utr3_truncated=substring(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end))
  
  utr3_seqs_lst <- utr3_seqs$utr3_truncated
  names(utr3_seqs_lst) <- utr3_seqs$peak_name
  
  return(utr3_seqs_lst)
}

get_genes <- function(peak_ids) {
  load(file.path(file_path_ref$project_directory, 
                 file_path_ref$results$quantseq_3end_gene_annotated))
  
  ensembl_transcript_ids <- subset(quantseq_utr_annotated_filt, 
                                  peak_name %in% peak_ids)$utr_name
  
  ensembl_gene <- biomaRt::getBM(c("ensembl_transcript_id", "ensembl_gene_id"),
                                 filters="ensembl_transcript_id", 
                                 values=ensembl_transcript_ids,
                                 mart=grcm38_101) %>%
    inner_join(quantseq_utr_annotated_filt, by=c("ensembl_transcript_id" = "utr_name"))
  
  ensembl_gene_lst <- ensembl_gene$ensembl_gene_id
  names(ensembl_gene_lst) <- ensembl_gene$peak_name
  
  return(ensembl_gene_lst)
}

run_pca <- function(norm_counts_df, metadata, ntop=10000) {
  # norm_counts_df looks like DESeq2's assay(vst(dds))
  
   # calculate the variance for each gene
  rv <- rowVars(norm_counts_df)

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(norm_counts_df[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  pc_df <- data.frame(pca$x)
  pc_df$sample_id <- rownames(pc_df)
  pc_df <- pc_df %>% inner_join(metadata, by="sample_id")
  return(list(pc_df=pc_df, percent_var=percentVar, pca=pca))
  
  
}



## ---------------------------------------------------------------------------------------------------------------

plot_pca_sizefactor_conversions <- function() {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$soma_counts_data))
  
  meta_soma <- subset(metadata, location == "Soma" & IUE_plasmid == "pPV25" & sample_id != "PV_SLAM_L1") # sample_time_min < 720)
  meta_soma$log_sample_time_h <- sapply(meta_soma$sample_time_min, function(t_min) ifelse(t_min == 0, 0, log2(t_min/60)))

  valid_genes_soma <- (counts_soma %>% 
    group_by(Geneid) %>%
    summarize(meets_thresh=sum(reads_total > 100)/length(reads_total) > 2/3 & sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(meets_thresh))$Geneid
  
  counts_soma_filt <- counts_soma %>% subset(Geneid %in% valid_genes_soma)
  
  t_counts_mat <- get_t_counts_mat(counts_soma_filt)
  
  dds_t_counts <- DESeqDataSetFromMatrix(t_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  dds_t_counts <- DESeq(dds_t_counts)
                             
  tc_counts_mat <- get_tc_counts_mat(counts_soma_filt)
  
  dds_tc_counts <- DESeqDataSetFromMatrix(tc_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  sizeFactors(dds_tc_counts) <- sizeFactors(dds_t_counts)
  
  
  dds_tc_counts <- DESeq(dds_tc_counts)
  
  norm_tc_counts <- assay(varianceStabilizingTransformation(dds_tc_counts, blind=TRUE)) #log2(counts(dds_tc_counts, normalize=TRUE)+1)
  
  pc_res <- run_pca(norm_tc_counts, meta_soma)
  
  ggplot(pc_res$pc_df, aes(PC1, PC2, label=sample_id, color=factor(sample_time_min))) + geom_point(size=2) + scale_color_manual(values=accent_colors[c(10, 8, 6, 5, 3, 1)]) +
    theme_classic()
  
  vst_cor <- cor(norm_tc_counts)
  library(pheatmap)
  annotation_row <- data.frame(sample_time_min=meta_soma[,c("sample_time_min")],
                               row.names=meta_soma$sample_id)
  annotation_row$sample_time_min <- sprintf("%s_min", annotation_row$sample_time_min)
  annotation_colors <- accent_colors[c(10, 8, 6, 5, 3, 1)]
  names(annotation_colors) <- sprintf("%s_min", c(0, 120, 240, 360, 480, 720))
  
  pheatmap(vst_cor, annotation_row=annotation_row, annotation_colors=list(sample_time_min=annotation_colors))
}




## ---------------------------------------------------------------------------------------------------------------
# load soma counts data, get conversions, and subtract 0h conversions
tc_convs_minus_bkg <- get_normalized_conversions()

# get averages
avg_tc_convs_minus_bkg <- tc_convs_minus_bkg %>% 
  group_by(Geneid, sample_type, sample_time_min) %>%
  summarize(mean_conversion_rate=mean(conversion_rate))

# get maximum timepoint
max_tc_convs_minus_bkg <- avg_tc_convs_minus_bkg %>%
  group_by(Geneid) %>%
  summarize(max_timepoint=sample_time_min[which.max(mean_conversion_rate)])

# get counts
max_tc_convs_minus_bkg %>%
  group_by(max_timepoint) %>%
  summarize(number=n(),
            proportion=n()/nrow(max_tc_convs_minus_bkg))

early_geneset <- get_genes(subset(max_tc_convs_minus_bkg, max_timepoint <= 240)$Geneid)
late_geneset <- get_genes(subset(max_tc_convs_minus_bkg, max_timepoint == 720)$Geneid)

# early genes are enriched for transcription factors
enrich_early_geneset <- enrichGO(early_geneset, keyType="ENSEMBL", ont="ALL", OrgDb=org.Mm.eg.db, universe=all_geneset)
# late genes enriched for "extracellular region/space" and "calcium ion binding"
enrich_late_geneset <- enrichGO(late_geneset, keyType="ENSEMBL", ont="ALL", OrgDb=org.Mm.eg.db, universe=all_geneset)

all_geneset <- get_genes(max_tc_convs_minus_bkg$Geneid)
enrich_early_geneset <- enrichGO(early_geneset, keyType="ENSEMBL", ont="ALL", OrgDb=org.Mm.eg.db, universe=all_geneset)



# now get the associated utrs for <= 240, 360, 480, and 720. Will run streme
time1_utrs <- get_utrs(subset(max_tc_convs_minus_bkg, max_timepoint <= 240)$Geneid)
time2_utrs <- get_utrs(subset(max_tc_convs_minus_bkg, max_timepoint == 360)$Geneid)
time3_utrs <- get_utrs(subset(max_tc_convs_minus_bkg, max_timepoint == 480)$Geneid)
time4_utrs <- get_utrs(subset(max_tc_convs_minus_bkg, max_timepoint == 720)$Geneid)

all_utrs <- c(time1_utrs, time2_utrs, time3_utrs, time4_utrs)

write.fasta(lapply(time1_utrs, function(x) strsplit(x, "")[[1]]), 
            names(time1_utrs), 
            file.path(file_path_ref$project_directory, "data/motifs/max_label_leq_4h.fa"),
            open="w", nbchar=10000000)

write.fasta(lapply(time2_utrs, function(x) strsplit(x, "")[[1]]), 
            names(time2_utrs), 
            file.path(file_path_ref$project_directory, "data/motifs/max_label_6h.fa"),
            open="w", nbchar=10000000)


write.fasta(lapply(time3_utrs, function(x) strsplit(x, "")[[1]]), 
            names(time3_utrs), 
            file.path(file_path_ref$project_directory, "data/motifs/max_label_8h.fa"),
            open="w", nbchar=10000000)


write.fasta(lapply(time4_utrs, function(x) strsplit(x, "")[[1]]), 
            names(time4_utrs), 
            file.path(file_path_ref$project_directory, "data/motifs/max_label_12h.fa"),
            open="w", nbchar=10000000)


write.fasta(lapply(all_utrs, function(x) strsplit(x, "")[[1]]), 
            names(all_utrs), 
            file.path(file_path_ref$project_directory, "data/motifs/all_utrs.fa"),
            open="w", nbchar=10000000)


## ---------------------------------------------------------------------------------------------------------------
signal_metadata <- subset(metadata, IUE_plasmid == "pPV25" & sample_time_min != 360)
counts_soma_signal <- subset(counts_soma, sample_id %in% signal_metadata$sample_id)

tc_conversions_pct_max <- get_rangenorm_conversions_mat(counts_soma_signal, sprintf("PV_SLAM_L%s", c(9, 18, 23)))

pca_rangenorm_conv <- run_pca(tc_conversions_pct_max, signal_metadata, ntop=500) 

ggplot(pca_rangenorm_conv$pc_df, aes(PC1, PC2, color=factor(sample_time_min))) + geom_point() + xlab(sprintf("PC1: %s Percent Variance Explained", round(100*pca_rangenorm_conv$percent_var[1]))) +
  ylab(sprintf("PC1: %s Percent Variance Explained", round(100*pca_rangenorm_conv$percent_var[2])))

annotations <- data.frame(sample_time_min=signal_metadata$sample_time_min)
annotations$sample_time_min <- factor(annotations$sample_time_min)
rownames(annotations) <- signal_metadata$sample_id
pheatmap(cor(tc_conversions_pct_max), annotation_row = annotations)
library(bios2mds)
library(NbClust)
library(factoextra)
library(cluster)

gap_stat <- clusGap(tc_conversions_pct_max, FUN = kmeans, nstart = 30, K.max = 20, B = 50)

fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")

fviz_nbclust(mammals_scaled, kmeans, method = "silhouette", k.max = 20) + theme_minimal() + ggtitle("The Silhouette Plot")



# Elbow method
run_elbow_method <- function(conversion_df, maxn=15) {
# Decide how many clusters to look at
n <- maxn

# Initialize total within sum of squares error: wss
wss <- numeric(n)

# Look over 1 to n possible clusters
for (i in 1:n) {
  # Fit the model: km.out
  km.out <- kmeans(conversion_df, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n, wss = wss)
 
scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    xlab('Number of clusters')
return(list(wss_df=wss_df, scree_plot=scree_plot))
}
tc_conversions_pct_max
clusters <- hclust(dist(iris[, 3:4]))
plot(clusters)

run_silhouette_method <- function(conversion_df, maxn=15) {
  res <- sil.score(conversion_df, nb.clus = c(2:maxn), nb.run = 100, iter.max = 1000,
method = "euclidean")
  
  return(list(silhouette_score=res))
}

sil_score <- run_silhouette_method(tc_conversions_pct_max)


## ---------------------------------------------------------------------------------------------------------------
fit_linear_exp_model_norm_conv <- function(tc_conv_df, plot=FALSE) {
  avg_conv_rate_df<- tc_conv_df %>% group_by(sample_time_min) %>%
    summarize(avg_conv_rate=mean(conversion_rate)) %>%
    rowwise() %>%
    mutate(x=sample_time_min/60)
  max_conv_rate <- max(avg_conv_rate_df$avg_conv_rate)

  tc_conv_df <- tc_conv_df %>%
    subset(conversion_rate >= 0) %>%
    rowwise() %>%
    mutate(x=sample_time_min/60,
           y=log(conversion_rate/max_conv_rate+1e-4))
   
   
  
   # remove the increasing values, focus only on decay
   t0 <- avg_conv_rate_df$x[which.max(avg_conv_rate_df$avg_conv_rate)][1]
   tc_conv_df <- subset(tc_conv_df, x >= t0) %>%
     rowwise() %>%
     mutate(x=x-t0)
   
   if(length(unique(tc_conv_df$x)) == 1) {
     return(data.frame(intercept=NA,
                       slope=NA,
                       rsquared=NA,
                       adj_rsquared=NA,
                       t0=NA))
   }
   
   fit <- lm(y~x, data=tc_conv_df)
   #print(summary(fit))
   
   #plot(ggplot(tc_conv_df, aes(x, y)) + 
  #        geom_point() +
  #        geom_abline(intercept=summary(fit)$coefficients[1,1],
  #                    slope=summary(fit)$coefficients[2, 1]))
   
   fit_summary <- summary(fit)
   
   
   if (plot == TRUE) {
     plot(ggplot(tc_conv_df, aes(x, y)) + geom_point() + geom_abline(intercept=fit_summary$coefficients[1,1],
                                                                     slope=fit_summary$coefficients[2,1]))
   }
   return(data.frame(intercept=fit_summary$coefficients[1,1],
                     slope=fit_summary$coefficients[2,1], 
                     rsquared=fit_summary$r.squared,
                     adj_rsquared=fit_summary$adj.r.squared,
                     t0=t0))
}

loglinear_model_fit <- data.frame(do.call(rbind, lapply(
  unique(tc_convs_minus_bkg$Geneid),
  function(peak) fit_linear_exp_model_norm_conv(
    subset(tc_convs_minus_bkg, Geneid == peak)) %>%
    rowwise() %>%
    mutate(peak_name=peak))
)
)

slope_ordered <- loglinear_model_fit %>% 
  subset(!is.na(slope)) %>%
  arrange(slope)

ordered_3utrs <- get_utrs(slope_ordered$peak_name)
slope_ordered_3utrs <- na.omit(ordered_3utrs[slope_ordered$peak_name])

kmer_spma <- run_kmer_spma(slope_ordered_3utrs)
kmer_spma_10bin <- run_kmer_spma(slope_ordered_3utrs, n_bins = 10)
kmer_spma_20bin <- run_kmer_spma(slope_ordered_3utrs, n_bins = 20)
matrix_spma_20bin <- run_matrix_spma(slope_ordered_3utrs, n_bins = 20)

get_lm_res <- function(df) { 
 mod <- lm(enrichment~bin, data=df)
 slope <- summary(mod)$coefficients[2,1]
 rsquared <- summary(mod)$r.squared
 return(data.frame(slope=slope, rsquared=rsquared))
}
all_kmer_enrich <- data.frame(do.call(rbind, lapply(1:20, function(i) kmer_spma_20bin$foreground_scores[[i]]$enrichment_df %>% rowwise() %>% mutate(bin=i) )))
all_kmer_lms <- all_kmer_enrich %>% group_by(kmer) %>% do(get_lm_res(.))

plot_kmer_bins <- function(k, all_kmer_enrich) {
  ggplot(subset(all_kmer_enrich, kmer==k),
         aes(bin, enrichment)) + geom_point()
}

save(kmer_spma, kmer_spma_10bin, kmer_spma_20bin, matrix_spma, matrix_spma_10bin, matrix_spma_20bin, file=file.path(file_path_ref$project_directory, "data/motifs/kmer_matrix_spma_based_on_soma_degradation_rates.RData")) 

# cannot use gene set enrichment just to group genes by stability
ordered_genes <- data.frame(ens_gene=get_genes(slope_ordered$peak_name)[slope_ordered$peak_name],
                            peak_name=slope_ordered$peak_name) %>%
  inner_join(slope_ordered, by="peak_name") %>%
  arrange(desc(slope)) %>%
  subset(!is.na(slope) & !is.na(ens_gene))

ordered_gene_slopes <- ordered_genes$slope
names(ordered_gene_slopes) <- ordered_genes$ens_gene

gse_res <- gseGO(ordered_gene_slopes, ont="ALL", OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")

early_up <- (dge_timepoints[["time_4h_8h"]] %>% slice_min(log2FoldChange, n=500))$peak_name
ggplot(loglinear_model_fit %>% rowwise() %>% mutate(is_early=peak_name %in% early_up), aes(slope, color=is_early)) + geom_histogram(bins=50, position='identity', alpha=0.5) 

fit_exponential_model_norm_conv <- function(tc_conv_df) {
  
  tc_conv_df <- tc_conv_df %>%
    rowwise() %>%
    mutate(x=sample_time_min/60,
           y=conversion_rate)
  
  Asym <- max(tc_conv_df$y)
  R0 <- 
  fit <- nls(y ~)
  
  tc_convs_minus_bkg <- get_normalized_conversions()
  
  fit <- nls(y ~ SSasymp(x, Asym, R0, lrc), data = dt)

}



## ---------------------------------------------------------------------------------------------------------------
 pc_res <- run_pca(tc_conversion_df_minus_bkg, meta_soma)
  ggplot(pc_res$pc_df, aes(PC1, PC2, label=sample_id, color=factor(sample_time_min))) + geom_point(size=2) + scale_color_manual(values=accent_colors[c(10, 8, 6, 5, 3, 1)]) +
    theme_classic()
  
  conv_cor <- cor(tc_conversion_df_minus_bkg)
  annotation_row <- data.frame(sample_time_min=meta_soma[,c("sample_time_min")],
                               row.names=meta_soma$sample_id)
  annotation_row$sample_time_min <- sprintf("%s_min", annotation_row$sample_time_min)
  annotation_colors <- accent_colors[c(10, 8, 6, 5, 3, 1)]
  names(annotation_colors) <- sprintf("%s_min", c(0, 120, 240, 360, 480, 720))
  
  pheatmap(conv_cor, annotation_row=annotation_row, annotation_colors=list(sample_time_min=annotation_colors))


## ---------------------------------------------------------------------------------------------------------------
get_relative_decay_480_720 <- function() {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$soma_counts_data))
  
  meta_soma <- subset(metadata, location == "Soma" & IUE_plasmid == "pPV25" & sample_time_min >= 480)
  meta_soma$log_sample_time_h <- sapply(meta_soma$sample_time_min, function(t_min) ifelse(t_min == 0, 0, log2(t_min/60)))

  valid_genes_soma <- (counts_soma %>% 
    group_by(Geneid) %>%
    summarize(meets_thresh=sum(reads_total > 100)/length(reads_total) > 2/3 & sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(meets_thresh))$Geneid
  
  counts_soma_filt <- counts_soma %>% subset(Geneid %in% valid_genes_soma)
  
  t_counts_mat <- get_t_counts_mat(counts_soma_filt)
  
  dds_t_counts <- DESeqDataSetFromMatrix(t_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  dds_t_counts <- DESeq(dds_t_counts)
                             
  tc_counts_mat <- get_tc_counts_mat(counts_soma_filt)
  
  dds_tc_counts <- DESeqDataSetFromMatrix(tc_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  sizeFactors(dds_tc_counts) <- sizeFactors(dds_t_counts)
  
  
  dds_tc_counts <- DESeq(dds_tc_counts)
  
  res <- results(dds_tc_counts, name="log_sample_time_h")
  
  load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))

 res_annot <- data.frame(res) 
 res_annot$peak_name <- rownames(res_annot)
 res_annot <- res_annot %>% 
   inner_join(quantseq_utr_annotated_filt)
 
 # get the UTRs
 library(biomaRt)
 grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")
 
 annotated_utrs <- biomaRt::getBM(c("ensembl_transcript_id", "3utr"), filters="ensembl_transcript_id", values=unique(res_annot$utr_name), mart=grcm38_101)
        
 res_annot_utr <- res_annot %>% 
   inner_join(annotated_utrs, by=c("utr_name" = "ensembl_transcript_id")) %>%
   subset(!grepl("N|S|s", `3utr`) & nchar(`3utr`) > 10) %>%
   rowwise() %>%
   mutate(utr_seq_truncated=substring(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end))
 
 ordered_res <- res_annot_utr %>%
   arrange(desc(log2FoldChange))
 
 utr_lst_ordered <- ordered_res$utr_seq_truncated
 names(utr_lst_ordered) <- ordered_res$peak_name
 
 library(transite)
 utr_20bin_spma <- run_matrix_spma(utr_lst_ordered, n_bins=20)
 utr_10bin_spma <- run_matrix_spma(utr_lst_ordered, n_bins=10)

 save(res_annot_utr, utr_20bin_spma, utr_10bin_spma, file=file.path(file_path_ref$project_directory,
                                     file_path_ref$results$soma_spma))
}


## ---------------------------------------------------------------------------------------------------------------
get_valid_genes_utrs <- function() {
  # load soma counts data, only indexed by "Geneid" which is the same as "peak_name"
  load(file.path(file_path_ref$project_directory, 
            file_path_ref$results$soma_counts_data))
  
  # load quantseq utr location definition, includes mapping of quantseq peaks to 
  # the respective 3'UTRs, with annotation for internal peaks of where the cutoff is
  load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))
  
  
valid_genes_soma <- (counts_soma %>% 
    group_by(Geneid) %>%
    summarize(meets_thresh=sum(reads_total > 100)/length(reads_total) > 2/3 & sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(meets_thresh))$Geneid
  
  counts_soma_filt <- counts_soma %>% 
    subset(Geneid %in% valid_genes_soma) %>%
     inner_join(quantseq_utr_annotated_filt,
                by=c("Geneid" = "peak_name"))
  
 # get the UTRs
 library(biomaRt)
 grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")
 
 annotated_utrs <- biomaRt::getBM(c("ensembl_transcript_id", "3utr"), filters="ensembl_transcript_id", values=unique(counts_soma_filt$utr_name), mart=grcm38_101) 
 
 annotated_utrs <- annotated_utrs %>%
   inner_join(quantseq_utr_annotated_filt, 
              by=c("ensembl_transcript_id" = "utr_name")) %>%
   subset(!grepl("N|S|s", `3utr`) & nchar(`3utr`) > 10) %>%
   rowwise() %>%
   mutate(utr_seq_truncated=substring(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end))
 
 return(list(valid_genes_soma=valid_genes_soma, 
             counts_soma_filt=counts_soma_filt,
             annotated_utrs=annotated_utrs))
}

valid_genes_utrs <- get_valid_genes_utrs()
save(valid_genes_utrs, file=file.path(file_path_ref$project_directory, 
                                      "data/sequencing/soma_utr_info.RData"))
                                      #file_path_ref$results$soma_utr_info))

get_dge_timepoints <- function(time0, time1) {
   
  # load valid_genes_utrs which contains
  # counts_soma_filt
  load(file.path(file_path_ref$project_directory, 
                                      "data/sequencing/soma_utr_info.RData"))
                                      #file_path_ref$results$soma_utr_info)))
  counts_soma_filt <- valid_genes_utrs$counts_soma_filt
  
  meta_soma <- subset(metadata, location == "Soma" & IUE_plasmid == "pPV25" & sample_time_min %in% c(time0, time1))
  meta_soma$log_sample_time_h <- sapply(meta_soma$sample_time_min, function(t_min) ifelse(t_min == time0, 0, 1))


  
  t_counts_mat <- get_t_counts_mat(counts_soma_filt)
  
  dds_t_counts <- DESeqDataSetFromMatrix(t_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  dds_t_counts <- DESeq(dds_t_counts)
                             
  tc_counts_mat <- get_tc_counts_mat(counts_soma_filt)
  
  dds_tc_counts <- DESeqDataSetFromMatrix(tc_counts_mat[,meta_soma$sample_id], 
                                colData=meta_soma, 
                                design=~log_sample_time_h)
  
  sizeFactors(dds_tc_counts) <- sizeFactors(dds_t_counts)
  
  
  dds_tc_counts <- DESeq(dds_tc_counts)
  
  res <- results(dds_tc_counts, name="log_sample_time_h")
  
 res_annot <- data.frame(res) 
 res_annot$peak_name <- rownames(res_annot)
 res_annot <- res_annot %>% 
   inner_join(valid_genes_utrs$annotated_utr, by="peak_name")
 
 return(res_annot)
}

get_relative_motif_spectrum <- function(time0, time1, 
                                        counts_soma_filt, annotated_utrs) {
  
 res_annot <- get_dge_timepoints(time0, time1, counts_soma_filt)
        
 res_annot_utr <- res_annot %>% 
   inner_join(annotated_utrs, by=c("utr_name" = "ensembl_transcript_id")) %>%
   subset(!grepl("N|S|s", `3utr`) & nchar(`3utr`) > 10) %>%
   rowwise() %>%
   mutate(utr_seq_truncated=substring(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end))
 
 ordered_res <- res_annot_utr %>%
   arrange(desc(log2FoldChange))
 
 utr_lst_ordered <- ordered_res$utr_seq_truncated
 names(utr_lst_ordered) <- ordered_res$peak_name
 
 library(transite)
 #utr_20bin_spma <- run_matrix_spma(utr_lst_ordered, n_bins=20)
 utr_10bin_spma <- run_matrix_spma(utr_lst_ordered, n_bins=10)
 
 return(list(utr_10bin_spma=utr_10bin_spma, dge_res=res_annot))

 }
  
time_comparisons <- list(time_2h_4h=c(120, 240),
                         time_2h_6h=c(120, 360),
                         time_4h_6h=c(240, 360),
                         time_4h_8h=c(240, 480),
                         time_6h_8h=c(360, 480),
                         time_8h_12h=c(480, 720))

dge_timepoints <- lapply(time_comparisons,
                         function(comp) get_dge_timepoints(comp[1], comp[2]))
names(dge_timepoints) <- names(time_comparisons)
save(dge_timepoints, file=file.path(file_path_ref$project_directory, 
                                      "data/sequencing/soma_dge.RData"))
                                      #file_path_ref$results$soma_dge))


utr_enrichment_deltas <- lapply(time_comparisons, 
                                function(comp) suppressWarnings(get_relative_motif_spectrum(
                                  comp[1], comp[2], 
                                  valid_genes_utrs$counts_soma_filt,
                                  valid_genes_utrs$annotated_utrs)))
names(utr_enrichment_deltas) <- names(time_comparisons)
                                      

save(utr_enrichment_deltas, file=file.path(file_path_ref$project_directory, "data/motifs/soma_spma_all_deltas.RData"))



## ---------------------------------------------------------------------------------------------------------------
load(file.path(file_path_ref$project_directory, 
                                      "data/sequencing/soma_dge.RData"))
                                      #file_path_ref$results$soma_dge)))

for (comparison in names(dge_timepoints)) {
  comparison_dge <- dge_timepoints[[comparison]]
  all_seqs <- comparison_dge$utr_seq_truncated
  t0_seqs <- comparison_dge %>% slice_min(log2FoldChange, n=500)
  
  
}


## ---------------------------------------------------------------------------------------------------------------
utr_scores <- utr_enrichment_deltas$time_8h_12h$utr_10bin_spma$background_scores

tc_conv_long_avg <- get_normalized_conversions() %>%
  group_by(Geneid) %>%
  summarize(maxval=max(conversion_rate)) %>% 
  inner_join(tc_conversion_df_minus_bkg_long, by="Geneid") %>%
  rowwise() %>%
  mutate(rangenorm_conv_rate=ifelse(maxval <= 0, 0, max(0, conversion_rate/maxval))) %>%
  group_by(Geneid, sample_type, IUE_plasmid, sample_time_min) %>%
  summarize(mean_rangenorm_conv_rate=mean(rangenorm_conv_rate))
  


## ---------------------------------------------------------------------------------------------------------------

get_motifs_per_gene <- function() {
  motif_names <- utr_enrichment_deltas$time_8h_12h$utr_10bin_spma$spectrum_info_df$motif_id
  utr_scores <- utr_enrichment_deltas$time_8h_12h$utr_10bin_spma$background_scores
  
  peak_names <- sapply(names(utr_scores$total_sites[[1]]),
                       function(peak_info) strsplit(as.character(peak_info), "\\|")[[1]][1])
  
  all_motif_hits <- data.frame(do.call(rbind,
          lapply(1:length(motif_names),
                 function(idx)
                 data.frame(peak_name=peak_names,
                            motif_counts=utr_scores$absolute_hits[[idx]],
                            motif_id=motif_names[[idx]]))))
  
  return(all_motif_hits)
                            
}

all_motif_hits <- get_motifs_per_gene()

plot_motif_txpt_genes <- function(motif_char) {
  
  idx <- which(utr_enrichment_deltas$time_8h_12h$utr_10bin_spma$spectrum_info_df$motif_id == motif_char)
  
motif_res <- data.frame(peak_info=names(utr_scores$total_sites[[idx]]),
                        motif_counts=utr_scores$absolute_hits[[idx]]) %>%
  rowwise() %>%
  mutate(peak_name=strsplit(as.character(peak_info), "\\|")[[1]][1])

motif_annot_conversions <- tc_conv_long_avg %>% 
  inner_join(motif_res, by=c("Geneid" = "peak_name"))

motif_annot_conversions_grouped <-
  motif_annot_conversions %>%
  rowwise() %>%
  mutate(num_motif=motif_counts) %>%
  group_by(sample_time_min, IUE_plasmid, sample_type, num_motif) %>%
  summarize(rate=mean(mean_rangenorm_conv_rate),
            upper_ci=quantile(mean_rangenorm_conv_rate, 0.975),
            lower_ci=quantile(mean_rangenorm_conv_rate, 0.025))

motif_annot_conversions_grouped$num_motif <- factor(motif_annot_conversions_grouped$num_motif, 
         levels=0:5)

plt <- ggplot(motif_annot_conversions_grouped, aes(sample_time_min, rate, group=num_motif, color=num_motif)) +
  geom_line() +
  scale_color_manual(values=accent_colors[c(10, 9, 7, 5, 3, 1)]) +
  ylab("Proportion of Maximum") +
  theme_classic()

return(plt)
}


## ---------------------------------------------------------------------------------------------------------------
#hnrnpr, syncrip
# CAAAAAG
plot_motif_txpt_genes("M236_0.6") + ggtitle("M236_0.6: Hnrnpr, Syncrip")
# AAAAAAG
plot_motif_txpt_genes("M243_0.6") + ggtitle("M243_0.6: Hnrnpr, Syncrip")
# CCAAAUU <- least clear
plot_motif_txpt_genes("M242_0.6") + ggtitle("M242_0.6: Hnrnpr, Syncrip")


# rbms3
# aUAUAua <- look at more AU repeats??
plot_motif_txpt_genes("M055_0.6") + ggtitle("M055_0.6: Rbms")
# uAUAUA
plot_motif_txpt_genes("M164_0.6") + ggtitle("M164_0.6: Rbms")
# uAUAu[a/c]c <- less pronounced
plot_motif_txpt_genes("M143_0.6") + ggtitle("M143_0.6: Rbms")


# sf3b4
plot_motif_txpt_genes("M195_0.6") + ggtitle("M195_0.6: Sf3b4")

# pcbp4
plot_motif_txpt_genes("M043_0.6") + ggtitle("M043_0.6: Pcbp4")

# ybx1
plot_motif_txpt_genes("1177_19561594") + ggtitle("1177_19561594: Ybx1")

# hnrnpk
plot_motif_txpt_genes("M026_0.6") + ggtitle("M026_0.6: Hnrnpk")

# srsf4
plot_motif_txpt_genes("M126_0.6") + ggtitle("M126_0.6: Srsf4")

plot_motif_txpt_genes("M141_0.6") + ggtitle("M141_0.6: ESRP2, ESRP1")

plot_motif_txpt_genes("M178_0.6") + ggtitle("M178_0.6: Celf3, Celf6")

plot_motif_txpt_genes("M105_0.6") + ggtitle("M105_0.6: Srsf1")

# khdrbs2
plot_motif_txpt_genes("M176_0.6") + ggtitle("M176_0.6: Khdrbs2")

# GC rich
plot_motif_txpt_genes("M050_0.6") + ggtitle("M050_0.6: Rbm4")
plot_motif_txpt_genes("M044_0.6") + ggtitle("M044_0.6")
plot_motif_txpt_genes("M054_0.6")+ ggtitle("M054_0.6")



## ---------------------------------------------------------------------------------------------------------------
get_axon_dge <- function(valid_genes) {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$axon_counts_data))
  
  
  meta_axon <- subset(metadata, location == "Axon")

  valid_axon_genes <- (counts_axon %>%
    group_by(Geneid) %>%
    summarize(is_valid=sum(tc_geq1_reads_total > 10) >= 3) %>%
    subset(is_valid))$Geneid
  
  counts_axon_filt <- counts_axon %>% subset(Geneid %in% intersect(valid_genes, valid_axon_genes))
    

  t_counts_mat <- get_t_counts_mat(counts_axon_filt)
  print(dim(t_counts_mat))
  
  dds_t_counts <- DESeqDataSetFromMatrix(t_counts_mat[,meta_axon$sample_id], 
                                colData=meta_axon, 
                                design=~sample_type)
  
  dds_t_counts <- DESeq(dds_t_counts)
                             
  tc_counts_mat <- get_tc_counts_mat(counts_axon_filt)
  
  dds_tc_counts <- DESeqDataSetFromMatrix(tc_counts_mat[,meta_axon$sample_id], 
                                colData=meta_axon, 
                                design=~sample_type)
  
  sizeFactors(dds_tc_counts) <- sizeFactors(dds_t_counts)
  
  
  dds_tc_counts <- DESeq(dds_tc_counts)
  print(resultsNames(dds_tc_counts))
  
  res <- results(dds_tc_counts, name="sample_type_UPRT_Axon_vs_Neg_Ctrl_Axon")
  
  load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))

 res_annot <- data.frame(res) 
 res_annot$peak_name <- rownames(res_annot)
 res_annot <- res_annot %>% 
   inner_join(quantseq_utr_annotated_filt)
 
 return(list(dds=dds_tc_counts, dge_res=res_annot))
  
}
dge_axon <- get_axon_dge(valid_genes_utrs$valid_genes_soma)
get_axon_betabinom <- function(loc, time, valid_genes_soma=valid_genes_utrs$valid_genes_soma) {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$axon_counts_data))
  
  
  meta_axon <- subset(metadata, location == loc & sample_time_min == time & !grepl("R1", sample_id))
  
  print(meta_axon)

  valid_axon_genes <- intersect((counts_axon %>%
                                   subset(sample_id %in% meta_axon$sample_id) %>% 
    group_by(Geneid) %>%
    summarize(is_valid=sum(reads_total > 100)/length(reads_total) > 2/3 & sum(tc_geq1_reads_total > 10)/length(reads_total) >= 1/2) %>%
    subset(is_valid))$Geneid, valid_genes_soma)
  
  print(length(valid_axon_genes))
  
  all_counts_axon <- counts_axon %>% subset(sample_id %in% meta_axon$sample_id)

  betabinom_res <- data.frame(peak_name=valid_axon_genes,
                              pvalue=sapply(valid_axon_genes, function(x) get_ttest_peak(x, all_counts_axon)))
  betabinom_res$padj <- p.adjust(betabinom_res$pvalue, method='fdr')
  
   load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))

  betabinom_res <- betabinom_res %>%
   inner_join(quantseq_utr_annotated_filt)
 
  return(betabinom_res)
  
}
axon_betabinom_12h_dist <- suppressMessages(get_axon_betabinom("Distal_Axon", 720))
axon_betabinom_18h_dist <- suppressMessages(get_axon_betabinom("Distal_Axon", 1080))

axon_betabinom_pvals <- get_axon_betabinom()


## ---------------------------------------------------------------------------------------------------------------
get_gc_dge <- function(valid_genes) {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$gc_counts_data))
  
  
  meta_gc <- subset(metadata, location == "GC")

  valid_gc_genes <- (counts_gc %>%
    group_by(Geneid) %>%
    summarize(is_valid=sum(tc_geq1_reads_total > 2) >= 3) %>%
    subset(is_valid))$Geneid
  
  counts_gc_filt <- counts_gc %>% subset(Geneid %in% intersect(valid_genes, valid_gc_genes))
    

  t_counts_mat <- get_t_counts_mat(counts_gc_filt)
  print(dim(t_counts_mat))
  
  dds_t_counts <- DESeqDataSetFromMatrix(t_counts_mat[,meta_gc$sample_id], 
                                colData=meta_gc, 
                                design=~sample_type)
  
  dds_t_counts <- DESeq(dds_t_counts)
                             
  tc_counts_mat <- get_tc_counts_mat(counts_gc_filt)
  
  dds_tc_counts <- DESeqDataSetFromMatrix(tc_counts_mat[,meta_gc$sample_id], 
                                colData=meta_gc, 
                                design=~sample_type)
  
  sizeFactors(dds_tc_counts) <- sizeFactors(dds_t_counts)
  
  
  dds_tc_counts <- DESeq(dds_tc_counts)
  print(resultsNames(dds_tc_counts))
  
  res <- results(dds_tc_counts, name="sample_type_UPRT_GC_vs_Neg_Ctrl_GC")
  
  load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))

 res_annot <- data.frame(res) 
 res_annot$peak_name <- rownames(res_annot)
 res_annot <- res_annot %>% 
   inner_join(quantseq_utr_annotated_filt)
 
 return(list(dds=dds_tc_counts, dge_res=res_annot))
  
}
library(countdata)
get_ttest_peak <- function(peak_id, all_counts_df) {
  counts_df <- subset(all_counts_df, Geneid == peak_id)
  positive_tc <- subset(counts_df, IUE_plasmid == "pPV25")$TC_count_total
  positive_t <- subset(counts_df, IUE_plasmid == "pPV25")$T_count_total
  negative_tc <- subset(counts_df, IUE_plasmid == "dollar_51")$TC_count_total
  negative_t <- subset(counts_df, IUE_plasmid == "dollar_51")$T_count_total

  if (sum(positive_t > 0) == length(positive_t) & sum(negative_t > 0) == length(negative_t)) {
  betabinom <- bb.test(c(positive_tc, negative_tc), 
          c(positive_t, negative_t),
          c(rep("UPRT", length(positive_tc)), rep("Neg", length(negative_tc))),
          alternative="greater")
  
  return(betabinom$p.value) 
  } else {
    return(NA)
  }

}
get_gc_betabinom <- function() {
  load(file.path(file_path_ref$project_directory, file_path_ref$results$gc_counts_data))
  
  
  meta_gc <- subset(metadata, location == "GC")

  valid_gc_genes <- (counts_gc %>%
    group_by(Geneid) %>%
    summarize(is_valid=sum(tc_geq1_reads_total > 5) >= 3) %>%
    subset(is_valid))$Geneid
  
  betabinom_res <- data.frame(peak_name=valid_gc_genes,
                              pvalue=sapply(valid_gc_genes, get_ttest_peak))
  betabinom_res$padj <- p.adjust(betabinom_res$pvalue, method='fdr')
  
   load(file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))

  betabinom_res <- betabinom_res %>%
   inner_join(quantseq_utr_annotated_filt)
 
  return(betabinom_res)
  
}
gc_betabinom_pvals <- get_gc_betabinom()

plot_gene_soma <- function(peak_name) {
  ggplot(subset(counts_soma, Geneid == peak_name & sample_type=="UPRT_Soma"),
         aes(sample_time_min, TC_count_total/T_count_total)) +
    geom_point()
}
plot_gene_axon <- function(peak_name) {
  ggplot(subset(counts_axon, Geneid == peak_name),
         aes(sample_type, TC_count_total/T_count_total, color=factor(sample_time_min))) +
    geom_boxplot() +
    geom_point(position=position_dodge(width=1)) + 
    facet_wrap(vars(location))
}
plot_gene_gc <- function(peak_name) {
  ggplot(subset(counts_gc, Geneid == peak_name),
         aes(sample_type, TC_count_total/T_count_total)) +
    geom_boxplot() +
    geom_point()
}

# look for motif overrepresentation
motifs_gc <- subset(all_motif_hits, peak_name %in% subset(gc_betabinom_pvals, pvalue < 0.1)$peak_name) %>% subset(motif_counts > 0) %>% group_by(motif_id) %>% summarize(num_total_gc=sum(motif_counts), num_contain_motif_gc=n())

motifs_all <- all_motif_hits %>% subset(motif_counts > 0) %>% group_by(motif_id) %>% summarize(num_total_soma=sum(motif_counts), num_contain_motif_soma=n())

motifs_combo <- motifs_gc %>% inner_join(motifs_all, by="motif_id")
motifs_combo$total_genes_gc <- length(subset(gc_betabinom_pvals, pvalue < 0.1)$peak_name)
motifs_combo$total_genes_soma <- length(unique(all_motif_hits$peak_name))

motifs_combo <- motifs_combo %>%
  rowwise() %>%
  mutate(fisher_res=list(fisher.test(matrix(c(num_contain_motif_gc, num_contain_motif_soma, total_genes_gc-num_contain_motif_gc, total_genes_soma-num_contain_motif_soma), nrow=2)))) %>%
  rowwise() %>%
  mutate(fisher_pval=fisher_res$p.value,
         fisher_odds_ratio=fisher_res$estimate,
         fisher_ci_lower=fisher_res$conf.int[1],
         fisher_ci_upper=fisher_res$conf.int[2])

plot_gene_soma("peak_12470")

plot_gene_gc("peak_12470")

gc_dge <- get_gc_dge(valid_genes_utrs$valid_genes_soma)

