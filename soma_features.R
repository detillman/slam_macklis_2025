## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(tidyr)
library(seqinr)
library(pheatmap)
library(transite)

source("Y:/users/priyav/slam_final/scripts/helper_script.R")

# load the AREscore functionailty from lianoglou
source(file.path(file_path_ref$project_directory, 
                 file_path_ref$scripts$rna_features_scripts$arescore))


plot_subdir <- file.path(file_path_ref$project_directory, "plots/SLAM_Soma_RNA_Features")
data_subdir <- file.path(file_path_ref$project_directory, "data/SLAM_Soma_RNA_Features")

overwrite <- FALSE


## ---------------------------------------------------------------------------------------------------------------
bin_x_lab <- "Shorter < -- Halflife -- > Longer"
bin_colors <- accent_colors[c(2, 4, 5, 8, 10)]


## ---------------------------------------------------------------------------------------------------------------
load(file.path(file_path_ref$project_directory, 
               file_path_ref$results$quantseq_3end_gene_annotated))

t2g <- getBM(c("ensembl_gene_id", "ensembl_transcript_id"), 
             filters="ensembl_transcript_id", 
             values=quantseq_utr_annotated_filt$utr_name, mart=grcm38_101)


## ---------------------------------------------------------------------------------------------------------------
overwrite <- FALSE
utr_seqs_annot_file <- file.path(data_dir, "202308_soma_halflife_utr_seqs_annot.RData")

if (!(file.exists(utr_seqs_annot_file) | overwrite)) {
  
load(file.path(file_path_ref$project_directory, 
               file_path_ref$results$soma_halflife_data)) 

# get utr sequences
utr_seqs <- getBM(c("ensembl_transcript_id", "3utr"), 
                          filters="ensembl_transcript_id",
                          values=unique(cpn_halflife_filt$ensembl_transcript_id),
                   mart=grcm38_101) 

# remove sequences that are too short
utr_seqs_filt <- utr_seqs %>%
  subset(!grepl("S|N", `3utr`) & nchar(`3utr`) > 10)


#Next load the quantseq data to define 3'ends and trim from the 3'UTR, as well as intron information
load(file.path(file_path_ref$project_directory, 
               file_path_ref$results$quantseq_3end_gene_annotated))


utr_introns <- getBM(c("chromosome_name", "strand", "ensembl_transcript_id", "3_utr_start", "3_utr_end", "transcript_biotype"),
                     filters="ensembl_transcript_id",
                     values=unique(cpn_halflife_filt$ensembl_transcript_id),
                     mart=grcm38_101) %>%
  tidyr::drop_na() %>% 
  inner_join(cpn_halflife_filt, by=c("ensembl_transcript_id")) 


#For FIMO, STREME, and other bed-based utilities, write bed-files of full utrs 
#and windows within utrs.

## make bed file of utr introns
get_coords_trunc <- function(peak_info_df, maxintrons=10, window_size=NA, window_step=NA) {
  region_name <- peak_info_df$peak_name[1]
    #print(region_name)

  chromosome_name <- peak_info_df$chromosome_name
  utr_starts <- as.numeric(peak_info_df$`3_utr_start`)
  utr_ends <- as.numeric(peak_info_df$`3_utr_end`)
  strand <-peak_info_df$strand 
  quantseq_start <- as.numeric(peak_info_df$quantseq_start[1])
  quantseq_end <- as.numeric(peak_info_df$quantseq_end[1])
  
  if (length(utr_starts) > maxintrons) {
    return(data.frame(seqnames=character(0),
                      start=integer(0),
                      end=integer(0), 
                      name=character(0),
                      width=integer(0),
                      strand=character(0)))
  }
  utr_ranges=GRanges(seqnames=Rle(as.character(chromosome_name)), 
                                      strand=Rle(strand), 
                                      ranges=IRanges(start=utr_starts, end=utr_ends))
  
  if (strand[1] == "+") {
    
      quantseq_range=GRanges(seqnames=as.character(chromosome_name[1]),
                         strand=Rle(strand[1]),
                         ranges=IRanges(start=min(utr_starts), end=quantseq_end+100))
      
  } else {
    quantseq_range=GRanges(seqnames=as.character(chromosome_name[1]),
                         strand=Rle(strand[1]),
                         ranges=IRanges(start=quantseq_start-100, end=max(utr_ends)))
  }
  
   if (is.na(window_size)) {
    
    overlap <- data.frame(intersect(utr_ranges, quantseq_range))
    overlap$name <- region_name
  
    return(overlap[c("seqnames", "start", "end", "name", "width", "strand")])
  } else {
    overlap <- intersect(utr_ranges, quantseq_range)
    
    windows <- slidingWindows(overlap, window_size, step=window_step)
    windows_df <- data.frame(windows)
    windows_df$name <- sprintf("%s:%s", region_name, 1:nrow(windows_df))
    
    return(windows_df[c("seqnames", "start", "end", "name", "width", "strand")])

  }
 }

all_coord_info <- utr_introns %>% dplyr::select(-strand) %>% inner_join(quantseq_utr_annotated_filt, 
                                             by=c("chromosome_name", "peak_name"    ))
coords_bed_df <- data.frame(do.call(rbind,
                                    lapply(unique(all_coord_info$peak_name),
                                           function(name) get_coords_trunc(subset(all_coord_info, peak_name == name)))))
# the bigwigSummary utility requires unique region names, so arbitrarily append index
coords_bed_df$name <- sprintf("%s:%s", coords_bed_df$name, 1:nrow(coords_bed_df))
write.table(coords_bed_df, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE,
            file=file.path(data_subdir, "regions/full_utr_coords.bed"))

coords_bed_df$seqnames <- sprintf("chr%s", coords_bed_df$seqnames)
write.table(coords_bed_df, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE,
            file=file.path(fdata_subdir, "regions/full_utr_coords.ucsc.bed"))

coords_bed_window_df <- data.frame(do.call(rbind,
                                    lapply(unique(all_coord_info$peak_name),
                                           function(name) get_coords_trunc(subset(all_coord_info, peak_name == name), window_size=50, window_step=10))))
coords_bed_window_df$seqnames <- sprintf("chr%s", coords_bed_window_df$seqnames)
write.table(coords_bed_window_df, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE,
            file=file.path(data_subdir, "regions/full_utr_coords.window50bp.ucsc.bed"))


utr_introns_agg <- utr_introns %>%
  group_by(peak_name, ensembl_transcript_id, chromosome_name, transcript_biotype) %>%
  summarize(utr_coords=paste0(sprintf("%s-%s", `3_utr_start`, `3_utr_end`), collapse=";"),
            num_exons=length(`3_utr_start`))
            

utr_seqs_annot <- utr_introns_agg %>%
  inner_join(cpn_halflife_filt, by=c("peak_name", "ensembl_transcript_id")) %>%
  inner_join(quantseq_utr_annotated_filt, by=c("peak_name", "ensembl_transcript_id" = "utr_name")) %>%
  inner_join(utr_seqs_filt, by=c( "ensembl_transcript_id")) %>% 
  rowwise() %>% 
  mutate(utr3_trunc=substr(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end))

# bin the halflives for future analysis
halflife_quantiles <- quantile(utr_seqs_annot$halflife_Veeraraghavan_CPN, 
                               seq(0.1, 0.9, 0.1), na.rm=TRUE)
utr_seqs_annot$halflife_bin=sprintf("Bin_%s", 1:5)[.bincode(
  utr_seqs_annot$halflife_Veeraraghavan_CPN, 
  c(0, halflife_quantiles[c(2,4,6,8)], Inf))]


utr_seqs_annot$halflife_bin <- factor(utr_seqs_annot$halflife_bin, 
                                      levels=sprintf("Bin_%s", 1:5))
save(utr_seqs_annot, 
     file=utr_seqs_annot_file)

} else {
  load(utr_seqs_annot_file)
}



## ---------------------------------------------------------------------------------------------------------------
utrs_to_write <- subset(utr_seqs_annot, nchar(utr3_trunc) < 10000 & nchar(utr3_trunc) > 10)
utr_names <- utrs_to_write$peak_name
utrs <- as.list(utrs_to_write$utr3_trunc)

write.fasta(utrs, utr_names, 
            file.path(data_subdir, "202308_utr_sequences_cpn_valid_labeled_genes.fa"), 
                      open = "w",
                      nbchar = 10000, 
                      as.string = TRUE)



## ---------------------------------------------------------------------------------------------------------------
get_coords_trunc <- function(peak_info_df, maxintrons=10, window_size=NA, window_step=NA) {
  region_name <- peak_info_df$peak_name[1]
  chromosome_name <- peak_info_df$chromosome_name.x[1]
  utr_starts <- peak_info_df$`3_utr_start`
  utr_ends <- peak_info_df$`3_utr_end`
  strand <- peak_info_df$strand[1]
  quantseq_start <- peak_info_df$quantseq_start[1]
  quantseq_end <- peak_info_df$quantseq_end[1]
  
  if (length(utr_starts) > maxintrons) {
    return(data.frame(seqnames=character(0),
                      start=integer(0),
                      end=integer(0), 
                      name=character(0),
                      width=integer(0),
                      strand=character(0)))
  }
  
  utr_ranges=GRanges(seqnames=Rle(as.character(chromosome_name)), 
                                      strand=Rle(strand), 
                                      ranges=IRanges(utr_starts, utr_ends))
  

  if (strand == "+") {
    
      quantseq_range=GRanges(seqnames=as.character(chromosome_name),
                         strand=Rle(strand),
                         ranges=IRanges(min(utr_starts), quantseq_end+100))
      
  } else {
    quantseq_range=GRanges(seqnames=as.character(chromosome_name),
                         strand=Rle(strand),
                         ranges=IRanges(quantseq_start-100, max(utr_ends)))
  }
  
  if (is.na(window_size)) {
    overlap <- data.frame(intersect(utr_ranges, quantseq_range))
    overlap$name <- region_name
  
    return(overlap[c("seqnames", "start", "end", "name", "width", "strand")])
  } else {
    overlap <- intersect(utr_ranges, quantseq_range)
    windows <- slidingWindows(overlap, window_size, step=window_step)
    windows_df <- data.frame(windows)
    windows_df$name <- region_name
    
    return(windows_df[c("seqnames", "start", "end", "name", "width", "strand")])

  }
}

coords_bed_df <- data.frame(do.call(rbind,
                                    lapply(unique(utr_seqs_annot$peak_name),
                                           function(name) get_coords_trunc(subset(utr_seqs_annot, peak_name == name)))))
# the bigwigSummary utility requires unique region names, so arbitrarily append index
coords_bed_df$name <- sprintf("%s:%s", coords_bed_df$name, 1:nrow(coords_bed_df))
write.table(coords_bed_df, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE,
            file=file.path(data_subdir, "regions/full_utr_coords.bed"))

coords_bed_df$seqnames <- sprintf("chr%s", coords_bed_df$seqnames)
write.table(coords_bed_df, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE,
            file=file.path(data_subdir, "regions/full_utr_coords.ucsc.bed"))



## ---------------------------------------------------------------------------------------------------------------
load(file.path(file_path_ref$project_directory, "external_data/veeraraghavan_engmann_all_dge_results.RData"))

gc_txpts <- subset(dge_results$transcript_level$gc_vs_gcf, padj.CPN_GC_DPC23.vs.Forebrain_Background < 0.05 & log2FoldChange.CPN_GC_DPC23.vs.Forebrain_Background > 0.5)$longest_txpt

total_utrs <- unique(utr_seqs_annot$ensembl_transcript_id)

gc_txpt_annot <- utr_seqs_annot %>% rowwise() %>% 
  mutate(is_gc_enrich=ensembl_transcript_id %in% gc_txpts) %>% 
  subset(!is.na(halflife_bin))

gc_uniq <- unique(subset(gc_txpt_annot, 
                         is_gc_enrich)$ensembl_transcript_id)
total_uniq <- unique(gc_txpt_annot$ensembl_transcript_id)

gc_txpt_per_bin <- gc_txpt_annot %>% 
  group_by(halflife_bin) %>% 
  summarize(num_gc_bin=sum(is_gc_enrich),
            num_soma_bin=sum(!is_gc_enrich),
            num_gc_notbin=length(gc_uniq)-sum(is_gc_enrich),
            num_soma_notbin=length(setdiff(setdiff(total_uniq, gc_uniq), ensembl_transcript_id[!is_gc_enrich]))) %>%
  rowwise() %>%
  mutate(fisher_res=list(fisher.test(matrix(c(num_gc_bin, num_soma_bin, num_gc_notbin, num_soma_notbin), nrow=2)))) %>%
  rowwise()%>%
  mutate(pval=fisher_res$p.value,
         odds_ratio=fisher_res$estimate,
         lower_ci=fisher_res$conf.int[1],
         upper_ci=fisher_res$conf.int[2])

gc_txpt_per_bin$halflife_bin <- factor(gc_txpt_per_bin$halflife_bin, 
                                         levels=sprintf("Bin_%s", 1:5))


ggplot(gc_txpt_per_bin, aes(halflife_bin, odds_ratio, color=pval < 0.05)) + 
  geom_point(size=2)  + 
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,  width=0.1)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.2) + 
  theme_classic() + 
  ggtitle("Enrichment for GC-localized transcripts\nAcross Distinct Halflife Bins")   +
  scale_y_continuous(trans="log2", breaks=c(0.5, 0.75, 1.5, 2)) +
  scale_color_manual(values=c(background_colors[4], accent_colors[2])) +
  xlab(bin_x_lab) +
  theme(text=element_text(size=5))
ggsave(file.path(plot_subdir, "figure3_gc_txpt_enrichment_halflife_bins.pdf"),
       height=60, width=60, units="mm")



## ---------------------------------------------------------------------------------------------------------------
# can plot as violin plot with anova comparing means
utr_seqs_annot %>%
  rowwise() %>%
  mutate(length_utr=nchar(utr3_trunc)) %>%
  subset(!is.na(length_utr) & !is.na(halflife_bin)) %>%
  ggplot(aes(halflife_bin, length_utr, fill=halflife_bin)) +
  geom_violin(draw_quantiles=0.5, linewidth=0.1) +
  scale_y_log10() +
  theme_classic() +
  geom_hline(yintercept=median(nchar(utr_seqs_annot$utr3_trunc)),
             linetype="dashed", linewidth=0.1) +
  ylab("3'UTR Length") +
  stat_compare_means(method="anova") +
  theme(legend.position="None", text = element_text(size=7)) +
  scale_fill_manual(values=bin_colors) +
  xlab(bin_x_lab)
ggsave(file.path(plot_subdir, "figure2_utr_length_halflife_bins.pdf"),
       height=30, width=30, units="mm")



## ---------------------------------------------------------------------------------------------------------------
to_plot <- utr_seqs_annot %>%
  subset(transcript_biotype != "nonsense_mediated_decay") %>%
  rowwise() %>%
  mutate(contains_intron=num_exons > 1) 

pval <- wilcox.test(subset(to_plot, contains_intron)$halflife_Veeraraghavan_CPN, 
                    subset(to_plot, !contains_intron)$halflife_Veeraraghavan_CPN)$p.value

ggplot(to_plot, aes(halflife_Veeraraghavan_CPN, color=contains_intron))+
  stat_ecdf(linewidth=0.1) +
  theme_classic() +
  scale_x_log10() +
  ggtitle(sprintf("Compare 3'UTRs with and without Intron\nP-value: %s", pval)) +
  scale_color_manual(values=c(background_colors[4], accent_colors[2])) +
  xlab("Relative Halflife in Soma (h)") +
  theme(text=element_text(size=5))
ggsave(file.path(plot_subdir, "figure2_supplement_presence_intron_utr_halflife.pdf"),
       height=40, width=60, units="mm")


## ---------------------------------------------------------------------------------------------------------------
#read in results from vienna rnafold
folding_energy <- read.csv(file.path(data_subdir, "structure/utr_sequences_vienna_mfe.csv"), header = FALSE)

# change colnames
colnames(folding_energy) <- c("utr3_trunc_rna", "vienna_mfe")

folding_energy$utr3_trunc <- gsub("U", "T", folding_energy$utr3_trunc_rna)

# need to normalize the MFE to account for length
# use the approach used in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4236180/ (also used in Tushev et al.)
# formula: mfe_den=100*(mfe-mfe_ref_len)/(L-L0)

run_simulation <- FALSE

if (run_simulation) {
sim_dir <- file.path(data_subdir, "structure/mfe_simulation_seqs")

dir.create(sim_dir)

# create examples of 20,30,40,50...70 %GC
# and lengths 50, 100, 150, 200, 250, 300, 350, 500, 1000, 
pct_gc <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
#lengths <- c(50, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
lengths <- seq(50, 950, 50)
simulate_seq <- function(l, gc) {
paste(sample(c("A", "U", "G", "C"), size=l, replace=TRUE, prob=c((1-gc)/2, (1-gc)/2, gc/2, gc/2)),
      collapse="")
}

for (gc in c(0.5))  {
  for (len in lengths) {
    write.csv(sapply(1:100, function(x) simulate_seq(len, gc)),
              file=file.path(sim_dir, sprintf("simulated_sequences.%s_pct_gc.%s_length.csv", gc*100, len)),
              quote=FALSE, row.names = FALSE)
  }
}
}


## ---------------------------------------------------------------------------------------------------------------
null_dist_dir <- file.path(data_subdir, "structure/mfe_simulation_output")
poss_files <- list.files(null_dist_dir)
null_dist_files <- file.path(null_dist_dir,
                             poss_files[grepl(".csv", poss_files)])

sim_res <- data.frame(
  do.call(rbind,
          lapply(null_dist_files,
                            function(f) read.csv(f, header = FALSE) %>%
                   dplyr::rename(utr_seq=V1, MFE=V2) %>%
                              rowwise() %>%
                              mutate(length=as.numeric(gsub("_length", "", strsplit(f, "\\.")[[1]][3])),
                                     gc_content=as.numeric(gsub("_pct_gc", "", strsplit(f, "\\.")[[1]][2]))/100,
                                     fname=f))
  )
)

null_distribution <- sim_res %>%
  group_by(length, gc_content) %>%
  summarize(num_seq=n(),
            mean_MFE=mean(MFE),
            median_MFE=median(MFE),
            min_MFE=min(MFE),
            max_MFE=max(MFE),
            upper_CI_MFE=quantile(MFE, 0.975),
            lower_CI_MFE=quantile(MFE, 0.025))


color_scale <- colorRampPalette(c("white", accent_colors[8], "black"))(10)
ggplot(null_distribution, aes(as.numeric(length), mean_MFE, fill=gc_content, group=gc_content, color=gc_content)) + 
  geom_point(size=0.2) +
  geom_line(linewidth=0.1) +
  geom_ribbon(aes(ymin=lower_CI_MFE, ymax=upper_CI_MFE, group=gc_content),
              alpha=0.5, linewidth=0.1) +
  theme_classic() +
  theme(text=element_text(size=5)) +
  scale_color_gradient(low=color_scale[3], high=color_scale[8]) +
  scale_fill_gradient(low=color_scale[3], high=color_scale[8]) +
  xlab("UTR Length (bp)")
ggsave(file.path(plot_subdir, "figure2_supplement_MFE_folding_null.pdf"),
       height=40, width=60, units="mm")

# create linear model to interpolate
mfe_model_gc <- lm(MFE~length+gc_content, data=sim_res)
mfe_model_equimolar <- lm(MFE~length, data=subset(sim_res, gc_content == 0.5))

# function to calculate reference MFE for a particular sequence
calc_mfe_ref_equimolar <- function(s) {
  len <- nchar(s)
  coef <- summary(mfe_model_equimolar)$coef
  return(coef[1,1]+coef[2,1]*len)
}

calc_mfe_ref_gc <- function(s) {
  len <- nchar(s)
  gc_content <- str_count(s, "G|C")/nchar(s)
  coef <- summary(mfe_model_gc)$coef
  return(coef[1,1]+coef[2,1]*len+coef[3,1]*gc_content)
}

# function to calculate the normalized mfe
# input is sequence s with numeric mfe
calc_mfe_den <- function(s, mfe, mfe_ref_func) {
  
  seqlen <- nchar(s)
  mfe_ref <- mfe_ref_func(s)
  # function is MFEden=100*(MFE-MFE_ref_Len)/(Len-L_0) 
  # where L_0 was empirically determined to be 8
  return(100*(mfe-mfe_ref)/(seqlen-8))
  
}




## ---------------------------------------------------------------------------------------------------------------
norm_mfe_df <- folding_energy %>% 
  inner_join(utr_seqs_annot, by="utr3_trunc") %>%
  subset(!is.na(halflife_bin)) %>%
  rowwise() %>% 
  mutate(mfe_norm_equimolar=calc_mfe_den(utr3_trunc, vienna_mfe, calc_mfe_ref_equimolar),
         mfe_norm_gc=calc_mfe_den(utr3_trunc, vienna_mfe, calc_mfe_ref_gc))
         
norm_mfe_df %>%
  ggplot(aes(halflife_bin, mfe_norm_equimolar, fill=halflife_bin)) + 
  geom_violin(draw_quantiles = 0.5) +
  geom_hline(yintercept=median(norm_mfe_df$mfe_norm_equimolar, na.rm=TRUE),
             linetype="dashed") +
  stat_compare_means(method="anova", label.y=30) +
  stat_compare_means(aes(halflife_bin, mfe_norm_anova),
                     method="wilcox.test",
                     comparisons=list(c("Bin_1", "Bin_2"), c("Bin_2", "Bin_3"), c("Bin_3", "Bin_4"), c("Bin_4", "Bin_5"))) +
  theme_classic() +
  ggtitle("Normalized MFE of Folding\n(Vienna RNA)\nNormalized for Length") +
  ylab("Length-normalized MFE") +
  theme(legend.position="None") +
  scale_fill_manual(values=bin_colors) +
  xlab(bin_x_lab)



norm_mfe_df %>%
  ggplot(aes(halflife_bin, mfe_norm_gc, fill=halflife_bin)) + 
  geom_violin(draw_quantiles = 0.5) +
  stat_compare_means(method="anova", label.y=700) +
  stat_compare_means(aes(halflife_bin, mfe_norm_gc),
                     comparisons=list(c("Bin_1", "Bin_2"), c("Bin_2", "Bin_3"), c("Bin_3", "Bin_4"), c("Bin_4", "Bin_5")),
                     method="wilcox.test") +
  theme_classic() +
  geom_hline(yintercept=median(norm_mfe_df$mfe_norm_gc, na.rm = TRUE), linetype="dashed") +
  ggtitle("Normalized MFE of Folding\n(Vienna RNA)\nNormalized for Length and %GC") +
  ylab("Length- and %GC-normalized MFE") +
  ylim(-10, 10) +
  theme(legend.position="None") +
  scale_fill_manual(values=bin_colors) +
  xlab(bin_x_lab)

#CDF accomplishes this much better
norm_mfe_df %>%
  ggplot(aes(mfe_norm_gc, color=halflife_bin)) + 
  geom_vline(xintercept=median(norm_mfe_df$mfe_norm_gc, na.rm = TRUE), linetype="dashed", linewidth=0.1) +
  stat_ecdf(linewidth=0.1) + 
  theme_classic() +
  theme(text=element_text(size=5)) +
  xlim(-100, 100) + 
  xlab("Length- and %GC-normalized MFE") +
  scale_color_manual(values=bin_colors)
ggsave(file.path(plot_subdir, "figure2_halflife_bins_gc_length_normalized_mfe.pdf"),
       height=50, width=75, units="mm")


# format for regression
peak_covar_mfe <- rbind(norm_mfe_df %>%
  dplyr::select(c(peak_name, mfe_norm_equimolar)) %>%
  rowwise() %>%
  mutate(covar_name="norm_length_MFE_vienna") %>%
  dplyr::rename(covar_value=mfe_norm_equimolar),
  norm_mfe_df %>%
  dplyr::select(c(peak_name, mfe_norm_gc)) %>%
  rowwise() %>%
  mutate(covar_name="norm_length_GC_MFE_vienna") %>%
  dplyr::rename(covar_value=mfe_norm_gc))
  
  



## ---------------------------------------------------------------------------------------------------------------
peak_covar_nt_content <- data.frame(
  do.call(rbind, 
          lapply(c("A", "T", "G", "C"),
                 function(nt) utr_seqs_annot %>%
                   rowwise() %>%
                   mutate(covar_value=str_count(utr3_trunc, nt)/nchar(utr3_trunc)*100,
                          covar_name=sprintf("percent_%s", nt)
                   )
          )
  )
) %>%
  dplyr::select(c(peak_name, covar_value, covar_name))

nt_content_cor <- data.frame(peak_covar_nt_content %>% 
  inner_join(utr_seqs_annot, by="peak_name") %>% 
  dplyr::select(c(covar_name, covar_value, halflife_Veeraraghavan_CPN)) %>% 
  tidyr::drop_na() %>% group_by(covar_name) %>% 
  summarize(correlation=cor(covar_value, halflife_Veeraraghavan_CPN, method = "spearman")))
rownames(nt_content_cor) <- gsub("percent_", "", nt_content_cor$covar_name)
nt_content_cor$covar_name <- NULL
nt_content_cor <- t(nt_content_cor[c("A", "T", "G", "C"),])
colnames(nt_content_cor) <- c("A", "T", "G", "C")

peak_covar_dint_content <- data.frame(
  do.call(rbind, 
          lapply(c("AT|TA", "AA", "AG|GA", "AC|CA", 
                   "TT", "TG|GT", "TC|CT", 
                   "CC", "GC|CG", "GG"),
                 function(pattern) utr_seqs_annot %>%
                   rowwise() %>%
                   mutate(covar_value=str_count(utr3_trunc, pattern)/(nchar(utr3_trunc)-1)*100,
                          covar_name=sprintf("percent_%s", pattern)
                   )
          )
  )
) %>%
  dplyr::select(c(peak_name, covar_value, covar_name))



dint_cor_halflife <- data.frame(peak_covar_dint_content %>% 
  inner_join(utr_seqs_annot, by="peak_name") %>% 
  dplyr::select(c(covar_name, covar_value, halflife_Veeraraghavan_CPN)) %>% 
  tidyr::drop_na() %>% group_by(covar_name) %>% 
  summarize(correlation=cor(covar_value, halflife_Veeraraghavan_CPN, method = "spearman")) %>%
  rowwise() %>%
  mutate(dinuc1=strsplit(gsub("percent_", "", covar_name), "\\|")[[1]][1]) %>%
  rowwise() %>%
  mutate(nt1=substr(dinuc1, 1, 1),
         nt2=substr(dinuc1, 2, 2)) %>%
  tidyr::pivot_wider(id_cols=nt1, names_from=nt2, values_from=correlation))
  
  rownames(dint_cor_halflife) <- dint_cor_halflife$nt1
  dint_cor_halflife$nt1 <- NULL
  
  dint_cor_halflife <- dint_cor_halflife[c("A", "T", "G", "C"), c("A", "T", "G", "C")]

pdf(file.path(plot_subdir, "figure2_nt_and_dint_content_cor_halflife.pdf"),
    height=6/2.54, width=6/2.54)
  pheatmap(rbind(nt_content_cor, dint_cor_halflife), cluster_rows = FALSE, cluster_cols=FALSE, na_col=background_colors[4],
         color=accent_colors) 
dev.off()


## ---------------------------------------------------------------------------------------------------------------
library(enrichMiR)
# load mouse conserved
data("consTSmm", package="enrichMiR")
consTSmm

df <- utr_seqs_annot# %>%
  #subset(!is.na(halflife_bin)) %>%
  #group_by(external_gene_name) %>%
  #summarize(halflife_bin=halflife_bin[which.max(mean_peak_tpm)])

test_bin_enrichment <- function(bin_id) {
  is_bin <- df$halflife_bin == bin_id
  names(is_bin) <- df$external_gene_name
  er <- testEnrichment(x=is_bin, sets=consTSmm)
  return(er)
  
}

enrichment_bin <- data.frame(
  do.call(rbind,
          lapply(1:5,
                 function(i) data.frame(test_bin_enrichment(
                   sprintf("Bin_%s", i))$siteoverlap.features) %>%
                     rowwise() %>%
                     mutate(halflife_bin=sprintf("Bin_%s", i))
                 )
                 
  )
)
enrichment_bin$halflife_bin <- factor(enrichment_bin$halflife_bin,
                                        levels=sprintf("Bin_%s", 1:5))

mirnas_enriched <- subset(enrichment_bin, FDR < 0.1)$members


mirna_sig_enrich_wide <- data.frame(enrichment_bin %>% 
                                      subset(members %in% mirnas_enriched) %>%
                                      tidyr::pivot_wider(id_cols=members, names_from=halflife_bin, values_from=enrichment))
rownames(mirna_sig_enrich_wide) <- sapply(mirna_sig_enrich_wide$members, 
                                          function(x) str_wrap(gsub(";", " ", x), width=15))
mirna_sig_enrich_wide$members <- NULL

# use predef accent colors
better_col_palette <- colorRampPalette(accent_colors)(30)[30:1]
#dev.off()
pdf(file.path(plot_subdir, "miRNA_halflife_bin_enrichment.pdf"),
    height=6/2.54, width=8/2.54)
pheatmap(mirna_sig_enrich_wide, cluster_cols = FALSE, color = better_col_palette, fontsize = 5 )
dev.off()

# get go-term enrichment for targeted genes
mirna_genes_long <- enrichment_bin %>% subset(members %in% mirnas_enriched) %>%
  tidyr::unnest(genes.features) 

mirna_genes_long <- getBM(c("ensembl_gene_id", "external_gene_name"),
                          filters="external_gene_name",
                          values=mirna_genes_long$genes.features,
                          mart=grcm38_101) %>%
  inner_join(mirna_genes_long, by=c("external_gene_name" = "genes.features" ))

mirna_genes_list <- lapply(mirnas_enriched,
                           function(m) subset(mirna_genes_long, members == m)$ensembl_gene_id)
names(mirna_genes_list) <- mirnas_enriched
background <- unique(subset(utr_seqs_annot, !is.na(halflife_bin))$ensembl_gene_id)

mirna_gene_go <- compareCluster(mirna_genes_list, fun="enrichGO", ont="ALL", keyType="ENSEMBL", OrgDb=org.Mm.eg.db, universe=background)
dotplot(mirna_gene_go, showCategory=8) + 
  scale_color_gradient(low=accent_colors[2], high=background_colors[4]) + 
  theme(text=element_text(size=5))
ggsave(file=file.path(plot_subdir, "miRNA_target_gene_GO_enrichment.pdf"),
       height=150, width=150, units="mm")


# now check with mendonsa et al miRNAs that were sequenced from primary coritcal neurons (somata and neurites)
mendonsa_mirnas <- na.omit(read.csv(file.path(file_path_ref$project_directory, "data/external_data/mendonsa_mirna_primary_cortical_pn.csv"), header=TRUE))

mirna_site_enrich_vs_expression <- 
  enrichment_bin %>% 
  tidyr::separate_rows(members, sep = ";") %>%
  subset(grepl("miR", members)) %>%
  rowwise() %>% mutate(miRNA.family=tolower(paste(strsplit(members, "-")[[1]][2:3], collapse="-"))) %>%
  inner_join(mendonsa_mirnas) %>%
  tidyr::pivot_wider(id_cols=c(members, average_norm_counts_soma), names_from=halflife_bin, values_from=enrichment) %>%
  arrange(desc(average_norm_counts_soma))

to_plot <- data.frame(mirna_site_enrich_vs_expression)
rownames(to_plot) <- to_plot$members
to_plot$average_norm_counts_soma <- NULL
to_plot$members <- NULL
pheatmap(to_plot,  cluster_rows=FALSE, cluster_cols=FALSE, color=better_col_palette)



## ---------------------------------------------------------------------------------------------------------------
soma_spma_savefile <- file.path(data_subdir, "sequence/soma_spma.RData")

if (!file.exits(soma_spma_savefile) | overwrite) {
  # utr_seqs_trunc  
  utr_seqs_ordered <- utr_seqs_annot %>%
    subset(!is.na(halflife_Veeraraghavan_CPN)) %>%
    arrange(halflife_Veeraraghavan_CPN)
  
  utr_seqs_lst <- utr_seqs_ordered$utr3_trunc
  names(utr_seqs_lst) <- utr_seqs_ordered$peak_name
  
  halflife_spma_matrix20 <- suppressWarnings(run_matrix_spma(utr_seqs_lst, n_bins=20))
  save(halflife_spma_matrix20, file=soma_spma_savefile)
} else {
  load(soma_spma_savefile)
}

# plot
sig_motifs <- subset(halflife_spma_matrix20$spectrum_info_df, aggregate_classifier_score > 0)$motif_id
enrichment_scores <- data.frame(do.call(rbind, lapply(1:length(halflife_spma_matrix20$enrichment_dfs),
                                 function(i) subset(halflife_spma_matrix20$enrichment_dfs[[i]], motif_id %in% sig_motifs) %>%
                                   rowwise() %>% mutate(bin=i))) %>%
                                  pivot_wider(id_cols=motif_id, names_from=bin, values_from=enrichment))
rownames(enrichment_scores) <- enrichment_scores$motif_id
enrichment_scores$motif_id <- NULL
pdf(file.path(plot_subdir, "figure2_binned_spma_heatmap.pdf"),
    height = 8/2.54, width=12/2.54)
pheatmap(log(enrichment_scores), cluster_cols=FALSE,          
         color =colorRampPalette(background_colors[c(1, 2, 3, 7:10)])(50),
         cutree_rows = 2, fontsize = 5)
dev.off()





## ---------------------------------------------------------------------------------------------------------------
go_enrich_motif_spma <- function(spma_result, rbp_regex) {
  
  idx <- which(grepl(rbp_regex, spma_result$spectrum_info_df$motif_rbps))
  print(spma_result$spectrum_info_df$motif_rbps[idx])
  sum_motif_hits <- colSums(data.frame(do.call(rbind, 
                                               lapply(idx, function(i) spma_result$background_scores$absolute_hits[[i]] > 0))))
  
  motif_gene_info <- data.frame(is_motif_hit=sum_motif_hits > 0,
                                peak_name=sapply(names(spma_result$background_scores$total_sites[[1]]), function(x) strsplit(x, "\\|")[[1]][1])) %>%
    inner_join(utr_seqs_annot, by="peak_name") 
  
  print(motif_gene_info %>% group_by(is_motif_hit) %>% summarize(num_genes=length(unique(ensembl_gene_id))))
  
  return(enrichGO(unique(subset(motif_gene_info, is_motif_hit)$ensembl_gene_id),
                  OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont="ALL", 
                  universe=unique(motif_gene_info$ensembl_gene_id)))
  
}

rbfox_go_enrich <- go_enrich_motif(halflife_spma_matrix20, "A2BP1|RBFOX")
au_go_enrich <- go_enrich_motif(halflife_spma_matrix20, "RBMS")


## ---------------------------------------------------------------------------------------------------------------
m1 <- transite::get_motif_by_id("M055_0.6")[[1]]
m2 <- transite::get_motif_by_id("M164_0.6")[[1]]

pwm_to_pfm <- function(pwm) {
  
  return(2^pwm*0.25)
}

meme_header <- "MEME version 4\nALPHABET = ACGT\nstrands: +\n"
get_motif_string <- function(transite_motif) {
  m_pfm <- pwm_to_pfm(transite_motif@matrix)
  motif_name <- sprintf("MOTIF %s %s\n", 
                        transite_motif@id,
                        paste(transite_motif@rbps, collapse=";"))
  
  header_line <- sprintf("letter-probability matrix:\n")
  matrix_str <- paste0(sapply(1:nrow(m_pfm),
                             function(i) paste0(m_pfm[i,], collapse=" ")),
                      collapse="\n")
  
  return(paste0(motif_name, header_line, matrix_str, "", collapse="\n"))
}
 
str1 <- get_motif_string(m1) 
str2 <- get_motif_string(m2)

meme_file <- paste0(list(meme_header, str1, str2), collapse = "\n\n")  

fileConn<-file(file.path(data_subdir, "au_repeat_motifs.meme"))
writeLines(meme_file, fileConn)
close(fileConn)



## ---------------------------------------------------------------------------------------------------------------
motif_locs <- read.table(file.path(data_subdir, 
                                   "sequence/202308_fimo_au_repeat_utr_sequences_cpn_valid_labeled_genes.txt"), 
                         sep="\t", header = TRUE) 

collapse_motif_locs <- function(peak_name) {
  
  peak_motifs <- subset(motif_locs, sequence_name == peak_name)
  peak_motif_range <- data.frame(reduce(IRanges(start=peak_motifs$start, end=peak_motifs$stop)))
  utr_seq_full <- utr_seqs_annot$utr3_trunc[utr_seqs_annot$peak_name == peak_name]
  ensembl_gene <- utr_seqs_annot$ensembl_gene_id[utr_seqs_annot$peak_name == peak_name]
  if (length(peak_motif_range) > 0) {
    return(data.frame(max_width=max(peak_motif_range$width),
                      num_motifs=nrow(peak_motif_range),
                      total_coverage=sum(peak_motif_range$width),
                      starts=paste0(peak_motif_range$start, collapse=";"),
                      ends=paste0(peak_motif_range$end, collapse=";"),
                      max_start=peak_motif_range$start[which.max(peak_motif_range$width)],
                      max_end=peak_motif_range$end[which.max(peak_motif_range$width)],
                      utr_seqs=paste0(sapply(1:nrow(peak_motif_range), function(i) {
                        middle <- (peak_motif_range$start[i]+peak_motif_range$end[i])/2 
                        return(substr(utr_seq_full, max(0, middle-25), min(nchar(utr_seq_full), middle+25)))
                        }), collapse=";"),
                      peak_name=peak_name,
                      ensembl_gene_id=ensembl_gene))
  } else {
    return(data.frame(max_width=0,
                      num_motifs=0,
                      total_coverage=0,
                      max_start=NA,
                      max_end=NA,
                      utr_seqs="",
                      starts="",
                      ends="",
                      peak_name=peak_name,
                      ensembl_gene_id=ensembl_gene))
  }
  
}

motif_locs_info <- data.frame(do.call(rbind,
                                      lapply(unique(motif_locs$sequence_name), 
                                             collapse_motif_locs)))

# write motif_locs to a bed file as well
motif_loc_coords <- motif_locs_info %>% 
  dplyr::select(peak_name, starts, ends) %>%
  separate_longer_delim(c(starts, ends), delim=";") %>%
  inner_join(utr_seqs_annot, by=c( "peak_name")) %>%
  rowwise() %>%
  mutate(five_prime=ifelse(strand == "+", 
                              as.numeric(strsplit(utr_coords, "-")[[1]][1])+as.numeric(starts),
                              as.numeric(strsplit(utr_coords, "-")[[1]][2])-as.numeric(starts)),
         three_prime=ifelse(strand == "+", 
                              as.numeric(strsplit(utr_coords, "-")[[1]][1])+as.numeric(ends),
                              as.numeric(strsplit(utr_coords, "-")[[1]][2])-as.numeric(ends)),
         score=1,
         strand=strand,
         chromosome_name=sprintf("chr%s", chromosome_name.x)) %>%
  rowwise() %>%
  mutate(genomic_start=min(five_prime, three_prime),
         genomic_stop=max(five_prime, three_prime))

motif_loc_coords$sequence_name <- sprintf("%s_%s", motif_loc_coords$peak_name, 1:nrow(motif_loc_coords))
         
write.table(na.omit(unique(motif_loc_coords[,c("chromosome_name", "genomic_start", "genomic_stop", "sequence_name", "score", "strand")])),
            file=file.path(data_subdir, "sequence/202308_fimo_AU_repeats_locations.bed"),
            sep="\t", row.names = FALSE, col.names=FALSE, quote = FALSE)


half_window <- 25 # use 50bp windows centered at middle of motif
motif_loc_coords_window <- motif_locs %>% 
  dplyr::select(sequence_name, start, stop) %>%
  rowwise() %>%
  mutate(middle=(as.numeric(start)+as.numeric(stop))/2) %>% 
  inner_join(utr_seqs_annot, by=c( "sequence_name" = "peak_name")) %>%
  dplyr::rename(peak_name=sequence_name) %>% 
  rowwise() %>%
  mutate(five_prime=ifelse(strand == "+", 
                              as.numeric(strsplit(utr_coords, "-")[[1]][1])+(middle-half_window),
                              as.numeric(strsplit(utr_coords, "-")[[1]][2])-(middle-half_window)),
         three_prime=ifelse(strand == "+", 
                              as.numeric(strsplit(utr_coords, "-")[[1]][1])+(middle+half_window),
                              as.numeric(strsplit(utr_coords, "-")[[1]][2])-(middle+half_window)),
         score=1,
         strand=strand,
         chromosome_name=sprintf("chr%s", chromosome_name.x)) %>%
  rowwise() %>%
  mutate(genomic_start=min(five_prime, three_prime),
         genomic_stop=max(five_prime, three_prime))

motif_loc_coords_window$sequence_name <- sprintf("%s_%s", motif_loc_coords_window$peak_name, 1:nrow(motif_loc_coords_window))
         
write.table(na.omit(unique(motif_loc_coords_window[,c("chromosome_name", "genomic_start", "genomic_stop", "sequence_name", "score", "strand")])),
            file=file.path(data_subdir, "sequence/202308_fimo_AU_repeats_locations.50bp_window.bed"),
            sep="\t", row.names = FALSE, col.names=FALSE, quote = FALSE)

# now work backwards to get motif_loc_coords_window$sequence_name to width of motif


join_starts <- motif_locs_info[,c("max_start", "peak_name", "max_width")] %>% inner_join(motif_loc_coords_window, by=c("max_start" = "start", "peak_name"))
join_stops <- motif_locs_info[,c("max_end", "peak_name", "max_width")]%>% inner_join(motif_loc_coords_window, by=c("max_end" = "stop", "peak_name"))

peak_idx_to_width <- data.frame(do.call(rbind, lapply(unique(join_starts$peak_name),
                                 function(p) 
                                   { 
                                   max_start <- subset(motif_locs_info, peak_name == p)$max_start
                                   max_end <- subset(motif_locs_info, peak_name == p)$max_end
                                   motif_width <- subset(motif_locs_info, peak_name == p)$max_width

                                   valid_seqs <- subset(motif_loc_coords_window, start >= max_start & stop <= max_end)$sequence_name
                                   return(data.frame(peak_name=p,
                                                     motif_width=motif_width,
                                                     sequence_name=valid_seqs))
                                 }
)))


## ---------------------------------------------------------------------------------------------------------------

au_repeats_halflife <- utr_seqs_annot %>% 
  full_join(motif_locs_info, by=c("ensembl_gene_id", "peak_name")) %>%
  rowwise() %>%
  mutate(max_width=ifelse(is.na(max_width), 0, max_width),
         num_motifs=ifelse(is.na(num_motifs), 0, num_motifs),
         total_coverage=ifelse(is.na(total_coverage), 0, total_coverage),
         max_motif_dist_3end=nchar(utr3_trunc)-max_end,
         max_motif_dist_5end=max_start,
         rel_loc_motif=(max_end+max_start)/(2*nchar(utr3_trunc)))

# find the quantiles of the repeat lengths
rep_length_quant <- quantile(au_repeats_halflife$max_width, probs=seq(0, 1, 0.1))

au_bins <- c("0-6nt", "7-8nt", "9-10nt", ">10nt")
au_bin_int <- c(-1, 6, 8, 10, Inf)
au_repeats_halflife$au_repeat_bin <- factor(au_bins[
  .bincode(au_repeats_halflife$max_width,
           au_bin_int)], levels=au_bins) 


au_repeats_halflife %>%
  ggplot(aes(rel_loc_motif, color=au_repeat_bin)) +
  geom_density() + 
  scale_color_manual(values=bin_colors[c(4, 2, 1)]) +
  theme_classic() 

plt <- au_repeats_halflife %>%
  ggplot(aes(max_width)) +
  geom_histogram(bins=100, fill=background_colors[5]) +
  scale_y_log10() +
  theme_classic()
for (xval in au_bin_int[2:(length(au_bin_int)-1)]) {
  plt <- plt + geom_vline(xintercept=xval, linewidth=0.5)
}
  
au_repeats_halflife %>%
    ggplot(aes(halflife_Veeraraghavan_CPN, color=factor(au_repeat_bin))) +
    stat_ecdf(linewidth=0.1) +
    theme_classic() +
    xlab("Relative Halflife in Soma (h)") + 
  theme(text=element_text(size=5)) +
    scale_x_log10() + 
  scale_color_manual(values=bin_colors[c(5, 4, 2, 1)]) + 
    geom_vline(xintercept=median(utr_seqs_annot$halflife_Veeraraghavan_CPN, na.rm=TRUE), 
               linetype='dashed', linewidth=0.1, color=background_colors[5])

ggsave(file.path(plot_subdir, "figure3_au_repeats_cdf_plot.pdf"),
       height=50, width=75, units="mm")

# get go enrichment for longest categories
# 9-10 and 10+ both had the same trend for faster turnover
au_repeat_genes <-unique((au_repeats_halflife %>% subset(max_width > 8))$ensembl_gene_id)
au_repeat_enrich_go <- enrichGO(au_repeat_genes, OrgDb = org.Mm.eg.db, keyType="ENSEMBL", ont="CC", universe=unique(utr_seqs_annot$ensembl_gene_id))

dotplot(simplify(au_repeat_enrich_go)) + scale_color_gradient(low=accent_colors[2], high=background_colors[5]) + theme(text=element_text(size=5))
ggsave(file.path(plot_subdir, "figure3_au_repeats_go_terms.pdf"),
       height=100, width=150, units="mm")

# now look at number of motifs
ggplot(au_repeats_halflife %>% rowwise() %>% mutate(long=max_width > 8), aes(factor(num_motifs), halflife_Veeraraghavan_CPN)) + geom_violin(draw_quantiles = 0.5) + scale_y_log10() + geom_hline(yintercept=median(utr_seqs_annot$halflife_Veeraraghavan_CPN))

au_multiple_repeats <- enrichGO(subset(au_repeats_halflife, num_motifs > 3)$ensembl_gene_id,
                                OrgDb = org.Mm.eg.db, keyType="ENSEMBL", ont="CC", universe=unique(utr_seqs_annot$ensembl_gene_id))

au_multiple_repeat_compare <- compareCluster(list(Greater3=subset(au_repeats_halflife, num_motifs > 3)$ensembl_gene_id,Less3=subset(au_repeats_halflife, num_motifs <= 3)$ensembl_gene_id),
                                OrgDb = org.Mm.eg.db, keyType="ENSEMBL", ont="CC", universe=unique(utr_seqs_annot$ensembl_gene_id))


## ---------------------------------------------------------------------------------------------------------------
bwavgbed_colnames <- c("sequence_name", "size", "covered", "sum_phyloP", "mean0_phyloP", "mean_phyloP", "min_phyloP", "max_phyloP")

get_conservation_scores <- function(alignment_name) {
utr_conservation <- read.table(file.path(data_subdir, 
                                         sprintf("conservation_scores/full_utr_coords.window50bp.ucsc.60way_phyloP_%s.txt", alignment_name)), header=FALSE, col.names=bwavgbed_colnames, sep="\t") %>%
  rowwise() %>%
  mutate(peak_name=strsplit(sequence_name, ":")[[1]][1]) %>%
  group_by(peak_name) %>%
  summarize(mean_phyloP_bkg=mean(mean_phyloP),
            median_phyloP_bkg=median(mean_phyloP),
            min_phyloP_bkg=min(mean_phyloP),
            max_phyloP_bkg=max(mean_phyloP),
            quantile_75_bkg=quantile(mean_phyloP, 0.75),
            quantile_25_bkg=quantile(mean_phyloP, 0.25)) %>%
  rowwise() %>%
  mutate(au_repeat_bin="Background")


au_conservation <- read.table(file.path(data_subdir, 
                                         sprintf("conservation_scores/202308_fimo_AU_repeats_locations.50bp_window.60way_phyloP_%s.txt", alignment_name)), header=FALSE, col.names = bwavgbed_colnames, sep="\t") %>%
  inner_join(peak_idx_to_width, by="sequence_name") %>%
  group_by(peak_name, motif_width) %>%
  summarize(
            num_ranges=length(mean_phyloP),
            maxwindow_phyloP_longest_motif=max(mean_phyloP),             # get the maximal mean phyloP for the longest motif as determined by IDX (could be over multiple windows)
            meanwindow_phyloP_longest_motif=mean(mean_phyloP), # get the mean phyloP over the longest motif
            
            medianwindow_phyloP_longest_motif=median(mean_phyloP)
            ) %>%
  inner_join(utr_conservation, by="peak_name", suffix=c(".motif", ".bkg")) %>%
  rowwise() %>%
  mutate(obs_cons_vs_exp_cons=meanwindow_phyloP_longest_motif-mean_phyloP_bkg)

au_conservation$au_repeat_bin <- factor(au_bins[.bincode(au_conservation$motif_width, au_bin_int)], levels=au_bins)

return(list(utr_conservation=utr_conservation, au_conservation=au_conservation))
}

placental_conservation <- get_conservation_scores("placental")

vertebrate_conservation <- get_conservation_scores("vertebrate")

# plot comparative phyloP
df_motifs_au <- placental_conservation$au_conservation %>% 
        dplyr::rename(phyloP=medianwindow_phyloP_longest_motif) %>%
        ungroup() %>%
        dplyr::select(phyloP, au_repeat_bin)

df_utrs_au <- subset(placental_conservation$utr_conservation, 
             peak_name %in% placental_conservation$au_conservation$peak_name) %>% 
        dplyr::rename(phyloP=median_phyloP_bkg) %>%
        ungroup() %>% 
        dplyr::select(phyloP, au_repeat_bin) %>%
  rowwise() %>% mutate(au_repeat_bin="Background_AU_UTRs")


df_utrs_all <- subset(placental_conservation$utr_conservation) %>% 
        dplyr::rename(phyloP=median_phyloP_bkg) %>%
        ungroup() %>% 
        dplyr::select(phyloP, au_repeat_bin) %>%
  rowwise() %>% mutate(au_repeat_bin="Background_All_UTRs")

rbind(rbind(df_motifs_au,
      df_utrs_au),
df_utrs_all) %>%
  ggplot(aes(phyloP, color=au_repeat_bin)) + stat_ecdf() +
  scale_color_manual(values=c(bin_colors[c(4,2,1)], background_colors[3], background_colors[8])) +
  theme_classic() +
  xlab("Median(phyloP)")
ggsave(file.path(plot_subdir, "figure2_median_phyloP_motif_vs_AU_utr_vs_all_utr.50bp_windows.pdf"),
       height=50, width=130, units="mm")





## ---------------------------------------------------------------------------------------------------------------
utr_seqs_annot$are_score <- sapply(utr_seqs_annot$utr3_trunc, AREscore)
# use the are score algo for AUs instead
utr_seqs_annot$au_score <- sapply(utr_seqs_annot$utr3_trunc, function(seq) AREscore(seq, pentamer = "ATATA", overmer="ATATATATA"))

utr_seqs_annot$are_score_bin <- .bincode(utr_seqs_annot$are_score, c(-1, 0, 1.5, 3,  10, Inf), include.lowest=TRUE)

utr_seqs_annot %>%
  ggplot(aes(halflife_Veeraraghavan_CPN, color=factor(are_score_bin))) +
    stat_ecdf(linewidth=0.1) +
    theme_classic() +
    xlab("Relative Halflife in Soma (h)") + 
  theme(text=element_text(size=5)) +
    scale_x_log10() + 
  scale_color_manual(values=bin_colors[c(5, 4, 3, 2, 1)]) + 
    geom_vline(xintercept=median(utr_seqs_annot$halflife_Veeraraghavan_CPN, na.rm=TRUE), 
               linetype='dashed', linewidth=0.1, color=background_colors[5])
ggsave(file.path(plot_subdir, "figure3_are_score_cdf_plot.pdf"),
       height=50, width=75, units="mm")

utr_seqs_annot %>%
  ggplot(aes(are_score+0.1)) + 
  geom_histogram(bins=100) +
  scale_x_log10() +
  theme_classic() +
  xlab("ARE Score + 0.1") +
    geom_vline(xintercept=0.2) +
  geom_vline(xintercept=1.6) + 
  geom_vline(xintercept=3.1) +
  geom_vline(xintercept=10.1)
ggsave(file.path(plot_subdir, "figure3_are_score_distribution_bins.pdf"),
       height=50, width=75, units="mm")


## ---------------------------------------------------------------------------------------------------------------
are_score_savefile <- file.path(data_subdir, "sequence/ARE_scores_50bp_windows.RData")

if (!file.exists(are_score_savefile) | overwrite) {
  # grab the 50bp regions surrounding motifs
  au_regions <- unique(separate_rows(motif_locs_info, utr_seqs, sep=";")$utr_seqs)
  au_are_scores <- sapply(au_regions, AREscore)
  
  # get background ARE propensity in 3'UTRs (given 50bp regions)
  utr_seqs_50bp_window <- 
    utr_seqs_annot %>% 
    rowwise() %>%
    mutate(utr_seqs=paste0(sapply(0:(nchar(utr3_trunc)/10-1),
                                  function(i) substr(utr3_trunc, i*10, i*10+50)), collapse=";"))
  
  all_bp_regions <- unique(separate_rows(utr_seqs_50bp_window, utr_seqs, sep=";")$utr_seqs)
  # remove small regions
  all_bp_regions <- all_bp_regions[nchar(all_bp_regions) > 30]
  all_bp_are_scores <- sapply(all_bp_regions, AREscore)
  
  save(au_regions, au_are_scores, all_bp_regions, all_bp_are_scores, 
       file=are_score_savefile)
} else {
  load(are_score_savefile)
}

are_score_compare <- data.frame(are_score=c(all_bp_are_scores, au_are_scores),
           region=c(rep("All UTRs", length(all_bp_are_scores)), rep("AU(n) Elements", length(au_are_scores)))) 

are_score_compare$are_score_bin <- .bincode(are_score_compare$are_score, c(-1, 0, 1.5, 3,  10, Inf), include.lowest=TRUE)

total_sites_region <- are_score_compare %>%
  ungroup() %>%
  group_by(region) %>%
  summarize(total_sites=n())

total_sites_region %>%
  inner_join(are_score_compare, by="region") %>%
  group_by(are_score_bin, region) %>%
  summarize(prop_scores=n()/total_sites[1]) %>%
  ggplot(aes(factor(are_score_bin), prop_scores, fill=region)) + 
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  scale_fill_manual(values=c(background_colors[5], accent_colors[2])) +
  theme_classic() 
ggsave(file.path(plot_subdir, "figure2supp_are_score_distribution_compare_AU_windows_vs_3UTR.pdf"),
       height=50, width=75, units="mm")


## ---------------------------------------------------------------------------------------------------------------
are_au_compare <- au_repeats_halflife %>% inner_join(utr_seqs_annot[,c("peak_name", "are_score", "are_score_bin")]) 

cor(are_au_compare$are_score, are_au_compare$max_width, method = "spearman")

au_not_are <- unique((au_repeats_halflife %>% inner_join(utr_seqs_annot[,c("peak_name", "are_score", "are_score_bin")]) %>% subset(are_score_bin <= 2 & !(au_repeat_bin %in% c("0-6nt"))))$ensembl_gene_id.x)

goau_not_are <- enrichGO(au_not_are, OrgDb = org.Mm.eg.db, keyType="ENSEMBL", ont="ALL", universe=unique(t2g$ensembl_gene_id))

color_gradient <- colorRampPalette(c(accent_colors[2], background_colors[5]))(10)
dotplot(goau_not_are) + scale_color_gradient(low=color_gradient[1], high=color_gradient[8]) + theme(text=element_text(size=5))
ggsave(file.path(plot_subdir, "figure3_au_repeats_not_are_go_terms.pdf"),
       height=100, width=150, units="mm")

au_and_are <- unique((au_repeats_halflife %>% inner_join(utr_seqs_annot[,c("peak_name", "are_score", "are_score_bin")]) %>% subset(are_score_bin >=3 & !(au_repeat_bin %in% c("0-6nt"))))$ensembl_gene_id.x)
are_not_au <- unique((au_repeats_halflife %>% inner_join(utr_seqs_annot[,c("peak_name", "are_score", "are_score_bin")]) %>% subset(are_score_bin >=3 & au_repeat_bin == "0-6nt"))$ensembl_gene_id.x)
# look at AU, ARE, and 
au_are_go <- compareCluster(list(AU_and_ARE=au_and_are, AU_not_ARE=au_not_are, ARE_not_AU=are_not_au),
               OrgDb = org.Mm.eg.db, keyType="ENSEMBL", ont="ALL", universe=unique(t2g$ensembl_gene_id))

