## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------
library(rjson)


## ---------------------------------------------------------------------------------------------------------------
input <- "C:/Users/prv892/Documents/GitHub/slam/metadata.json"

file_path_ref <- fromJSON(file = input)
#source(file.path(file_path_ref$project_directory, file_path_ref$sequencing$setup_script))
#source(file.path(file_path_ref$project_directory, file_path_ref$plots$setup_script))

background_colors <- colorRampPalette(c("white", "black"))(10)
bkg_sample_colors <- background_colors[c(3,5,7)]


## ---------------------------------------------------------------------------------------------------------------
get_quantseq_result <- function(sample_id) {
quantseq_gz=gzfile(
  file.path(file_path_ref$project_directory,
            file_path_ref$results$quantseq_3end_definition,
            sprintf("%s_intersect_3utr.remove_lowmap.tsv.gz", sample_id)),'rt') 

quantseq_utr <- read.table(quantseq_gz, sep="\t")
# note that V1 is chromsome, V2 is the read 5' start, and V4 is the number of reads that have that 5' end
quantseq_utr_closest <- quantseq_utr %>% 
                group_by(V1, V2, V4) %>% 
  summarize(utr_start=V8[which.min(V14)], 
            utr_end=V9[which.min(V14)], 
            strand=V11[which.min(V14)],
            attributes=V13[which.min(V14)],
            fragment_length=V14[which.min(V14)]) %>%
  dplyr::rename(chromosome_name=V1,
                read5prime=V2,
                num_reads=V4)

quantseq_utr_closest$sample_id <- sample_id

return(quantseq_utr_closest)
}

quantseq_results <- data.frame(do.call(rbind, 
                                       lapply(c("PV_QuantSeq_L23", "PV_QuantSeq_L18", "PV_QuantSeq_L9"), 
                                              get_quantseq_result)))

save(quantseq_results, file=file.path(file_path_ref$project_directory,
                                 file_path_ref$results$quantseq_3end_definition))

plt <- ggplot(quantseq_results, aes(fragment_length, color=sample_id)) + 
    geom_density() + 
  scale_x_log10() + 
  geom_vline(xintercept=500, linetype="dashed") +
  geom_vline(xintercept=median(quantseq_results$fragment_length), linetype="dotted") + 
  theme_classic() + 
  ggtitle(sprintf("QuantSeq Estimated Fragment Size Distribution\nBased on Annotated 3'ends\nMedian: %s", median(quantseq_results$fragment_length))) +
  scale_color_manual(values=bkg_sample_colors)


ggsave(file.path(file_path_ref$project_directory, file_path_ref$plots$plot_directory, "quantseq_supplement_size_distribution.pdf"), plot=plt)


## ---------------------------------------------------------------------------------------------------------------

grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")

all_utrs <- getBM(c("ensembl_transcript_id", 
                         "ensembl_gene_id",
                         "chromosome_name",
                         "strand",
                         "3_utr_start", 
                         "3_utr_end"),
                       mart=grcm38_101) 

extended_utrs <-all_utrs %>% 
  subset(!is.na(`3_utr_start`) & !is.na(`3_utr_end`)) %>%
  group_by(ensembl_gene_id, chromosome_name, strand) %>% 
  summarize(end=ifelse(strand[1] > 0, max(`3_utr_end`)+1000, min(`3_utr_start`)),
            start=ifelse(strand[1] > 0, max(`3_utr_end`), min(`3_utr_start`)-1000)) %>%
  rowwise() %>%
  mutate(name=sprintf("%s_utr_extension", ensembl_gene_id)) %>%
  ungroup() %>% 
  dplyr::select(c(chromosome_name, start, end, name, strand)) 
extended_utrs$score <- "."
extended_utrs$strand <- ifelse(all_utrs_extend_bed$strand > 0, "+", "-")

# write the extension
write.table(extended_utrs[,c("chromosome_name", "start", "end", "name", "score", "strand")], sep="\t", quote=FALSE, col.names=FALSE, row.names = FALSE, file=file.path(file_path_ref$project_directory, "data/sequencing/grcm38_101_3utr_extended.bed"))

# write the normal utr bed
all_utrs <- all_utrs %>% subset(!is.na(`3_utr_start`) & !is.na(`3_utr_end`))
all_utrs$score <- "."
all_utrs$strand <- ifelse(all_utrs$strand > 0, "+", "-")
write.table(all_utrs[,c("chromosome_name", "3_utr_start", "3_utr_end", "ensembl_transcript_id", "score", "strand")],
            sep="\t", quote=FALSE, col.names=FALSE, row.names = FALSE,
            file=file.path(file_path_ref$project_directory, 
                           "data/sequencing/grcm38_101_3utr.bed"))


## ---------------------------------------------------------------------------------------------------------------
quantseq_utr_overlap <- read.table(file.path(file_path_ref$project_directory, "data/sequencing/PV_QuantSeq_peaks.intersect_3utr.bed"), sep="\t", header=FALSE)

quantseq_peak_named <- read.table(file.path(file_path_ref$project_directory, "data/sequencing/PV_QuantSeq_peaks.intersect_3utr.unique.bed"), sep="\t", header=FALSE)


colnames(quantseq_utr_overlap) <- c("chromosome_name", 
                                    "quantseq_start", "quantseq_end", 
                                    "num_reads", "chromosome_name2", 
                                    "utr_start", "utr_end", "utr_name", 
                                    "empty", "strand", "length_overlap")

colnames(quantseq_peak_named) <- c("chromosome_name", 
                                   "quantseq_start", "quantseq_end",
                                   "peak_name", "num_reads", "strand")

quantseq_utr_annotated <- quantseq_peak_named %>% 
  inner_join(quantseq_utr_overlap %>% 
               dplyr::select(c(chromosome_name, quantseq_start, quantseq_end, utr_start, utr_end, utr_name, length_overlap)))

# filter now to get only the relevant associated utr
quantseq_utr_annotated_filt <- quantseq_utr_annotated %>%
  rowwise() %>%
  mutate(dist_to_3end=ifelse(strand == "+", quantseq_end-utr_end, utr_start-quantseq_start)) %>%
  subset(utr_start > 0) %>% 
  group_by(chromosome_name, quantseq_start, quantseq_end, peak_name, num_reads, strand) %>%
  summarize(utr_start=utr_start[which.min(abs(dist_to_3end))],
            utr_end=utr_end[which.min(abs(dist_to_3end))],
            utr_name=utr_name[which.min(abs(dist_to_3end))],
            dist_to_3end=dist_to_3end[which.min(abs(dist_to_3end))],
            bp_to_remove_from_3end=ifelse(strand[1] == "+", max(0, utr_end-(quantseq_end+100)), max(0, (quantseq_start-100-utr_start))))

save(quantseq_utr_annotated, quantseq_utr_annotated_filt, 
     file=file.path(file_path_ref$project_directory,
                    file_path_ref$results$quantseq_3end_gene_annotated))


## ---------------------------------------------------------------------------------------------------------------
reads_per_utr_bysample <- quantseq_results %>% 
  rowwise() %>%
  mutate(is_fragment_geq500=ifelse(fragment_length > 500, "Fragment Length >500bp", "Fragment Length <500bp")) %>%
  group_by(chromosome_name, utr_start, utr_end, strand, attributes, sample_id, is_fragment_geq500) %>%
  summarize(num_reads=sum(num_reads))

ggplot(reads_per_utr_bysample, aes(num_reads, fill=sample_id)) + 
  geom_histogram() + 
  scale_x_log10() + 
  theme_classic() + 
  geom_vline(xintercept=10, linetype="dashed") +
  scale_fill_manual(values=bkg_sample_colors) + 
  facet_wrap(vars(is_fragment_geq500), scales="free_y") +
  xlab("Number of Reads Aligning to 3'End") + 
  ylab("Number of UTR 3'Ends")


## ---------------------------------------------------------------------------------------------------------------
unannotated_utrs <- subset(quantseq_results, fragment_length > 500) %>%
  group_by(attributes, sample_id) %>%
  summarize(num_reads=sum(num_reads)) %>%
  subset(num_reads > 10) %>%
  group_by(attributes) %>%
  summarize(n_samples_valid=n()) %>%
  subset(n_samples_valid >= 2)



unannotated_txptid <- sapply(unannotated_utrs$attributes,
                           function(attrib) {
                             res <- str_extract(as.character(attrib), 
                                       "transcript_id ENSMUST[:digit:]+;")
                             res <- gsub(";", "", gsub("transcript_id ", "", res))
                             return(res)
                           })
names(unannotated_txptid) <- NULL


grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")

unannot_utr_seqs <- getBM(c("ensembl_transcript_id", "3utr"), 
                          filters="ensembl_transcript_id",
                          values="")


## ---------------------------------------------------------------------------------------------------------------
utr3_coordinate_file <- file.path(file_path_ref$project_dir, file_path_ref$external_data$utr3_coordinates)

if (file.exists(utr3_coordinate_file)) {
  load(utr3_coordinate_file)
  
} else {
  # create ensembl biomart instance 
 grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")
 coords <- getBM(c("ensembl_gene_id", "ensembl_transcript_id",
                   "chromosome_name", "strand",
                   "3_utr_start", "3_utr_end"),
                 filters="ensembl_gene_id",
                 values=valid_genes,
                 mart=grcm38_101)
 
 coords_nona <- coords %>%
   subset(!is.na(`3_utr_start`) & !is.na(`3_utr_end`)) 
 

 # tested on get_last250(coords_nona, "ENSMUST00000000572") long last exon
 # can't use granges resize because this just auto adds more transcript real estate even if outside utr
 get_last250 <- function(coords_df, transcript_id) {
   
   txpt_df <- subset(coords_df, ensembl_transcript_id == transcript_id) 
   
   num_exons <- nrow(txpt_df)
   direction <- txpt_df$strand[1]

   
   if (num_exons == 1) {
     utr3_start <- txpt_df$`3_utr_start`[1]
     utr3_end <- txpt_df$`3_utr_end`[1]

     
     if (direction == 1) {
       start_end <- data.frame(start=max(utr3_start, utr3_end-250),
                               end=utr3_end,
                               rank=1,
                               strand=direction) %>%
         rowwise() %>%
         mutate(exon_length=abs(start-end))
       
     } else {
       start_end <- data.frame(start=utr3_start, 
                               end=min(utr3_end, utr3_start+250),
                               rank=1,
                               strand=direction
                               ) %>%
         rowwise() %>%
         mutate(exon_length=abs(start-end))
     }
   } else {
      
     exon_df <- txpt_df %>% dplyr::rename(start=`3_utr_start`,
                                          end=`3_utr_end`) %>%
       rowwise() %>%
       mutate(exon_length=abs(start-end)) %>%
       dplyr::select(-c(chromosome_name, ensembl_gene_id, ensembl_transcript_id))

     if (direction ==1) {
       exon_df_ordered <- exon_df %>% arrange(desc(start)) 
     } else {
       
       exon_df_ordered <- exon_df %>% arrange(start)
     }
       exon_df_ordered$rank <- 1:nrow(exon_df_ordered)

       
       total_coverage <- sapply(1:nrow(exon_df_ordered), 
                                function(i) sum(exon_df_ordered$exon_length[1:i]))

       whole_exons_include <- which(total_coverage < 250)
       if (length(whole_exons_include) == 0) {
         # case when last exon is geq 250
         last_exon <- 1
         start_end <- data.frame(start=integer(0),
                                 end=integer(0),
                                 rank=integer(0),
                                 strand=integer(0),
                                 exon_length=integer(0))
         bp_include_last_exon <- 250
         
       } else {
         last_exon <- max(whole_exons_include) + 1
         start_end <- exon_df_ordered[1:max(whole_exons_include),]
         bp_include_last_exon <- 250-total_coverage[max(whole_exons_include)]


       }
       
       
       if (nrow(exon_df_ordered) >= last_exon) {
       
         if (direction == 1) {
           last_exon_start= max(exon_df_ordered$start[last_exon], 
                              exon_df_ordered$end[last_exon]-bp_include_last_exon)
           last_exon_end = exon_df_ordered$end[last_exon]
         } else {
           last_exon_start = exon_df_ordered$start[last_exon]
           last_exon_end = min(exon_df_ordered$end[last_exon],
                             exon_df_ordered$start[last_exon]+bp_include_last_exon)
        
         }
         
         start_end <- data.frame(rbind(start_end, 
                                       data.frame(start=last_exon_start,
                                                end=last_exon_end,
                                                rank=last_exon,
                                                strand=direction) %>%
                                       rowwise() %>%
                                       mutate(exon_length=abs(start-end))))
       }
      

   }
   
   start_end$ensembl_gene_id <- txpt_df$ensembl_gene_id[1]
   start_end$ensembl_transcript_id <- txpt_df$ensembl_transcript_id[1]
   start_end$chromosome_name <- txpt_df$chromosome_name[1]
   
   return(start_end)
 }
 
 coords_last250 <- lapply(unique(coords_nona$ensembl_transcript_id),
                                             function(txid) get_last250(coords_nona, txid))
 
 coords_last250_df <- data.frame(do.call(rbind, coords_last250))

 coords_last250_granges <- coords_last250_df %>%
   group_by(ensembl_gene_id, ensembl_transcript_id) %>%
   summarize(txpt_ranges=list(GRanges(seqnames=Rle(as.character(chromosome_name)), 
                                      strand=Rle(strand(ifelse(strand > 0, "+", "-"))), 
                                      ranges=IRanges(start, end))))
 
 save(coords_last250_granges, file=file.path(file_path_ref$project_directory, "data/external_data/ensembl_GRCm38_101_coords_last250.RData"))
 
 # collapse overlapping transcripts
 collapse_gene <- function(gene_utrs_df) {

   all_combos <- expand.grid(gene_utrs_df$ensembl_transcript_id, 
                             gene_utrs_df$ensembl_transcript_id) %>%
   dplyr::rename(from=Var1, 
                to=Var2) %>%
     rowwise() %>%
     mutate(is_overlap=length(findOverlaps(subset(gene_utrs_df, ensembl_transcript_id == from)$txpt_ranges[[1]],
                                    subset(gene_utrs_df, ensembl_transcript_id == to)$txpt_ranges[[1]])) > 0) 
   
   overlap_combos <- data.frame(all_combos %>%
                                  subset(is_overlap))
   
   # get connected components
   comps <- components(graph_from_data_frame(overlap_combos))
   comp_groups <- groups(comps)
   
   # use connected components named by transcript id to grab the iranges
   # Genomic ranges for some reason does not support lapply etc
   group_names <- c()
   group_ranges <- c()
   for (i in comp_groups) {
     group_names <-  c(group_names, paste(i, collapse=";"))
     comp_ranges <- NULL
     # iter through txpts
     for (j in i) {
       if (is.null(comp_ranges)) {
         comp_ranges <- 
           subset(gene_utrs_df, ensembl_transcript_id == j)$txpt_ranges[[1]]
       } else {
       comp_ranges <- GenomicRanges::reduce(c(comp_ranges,
                        subset(gene_utrs_df, ensembl_transcript_id == j)$txpt_ranges[[1]]))
     }
     }
      group_ranges <- c(group_ranges, comp_ranges)

   }

   group_ranges_collapsed_df <- data.frame(
     do.call(rbind,
             lapply(1:length(group_ranges), 
                    function(x) data.frame(group_ranges[x]) %>%
                      rowwise() %>% mutate(ensembl_transcript_id=group_names[x])))
   )

   # add chromosome and gene info
   group_ranges_collapsed_df$chromosome_name <- group_ranges_collapsed_df$seqnames
   group_ranges_collapsed_df$seqnames <- NULL
   group_ranges_collapsed_df$ensembl_gene_id <- gene_utrs_df$ensembl_gene_id[1]
   #group_ranges_collapsed_df$strand <- gene_utrs_df$strand[1]


   
   return(group_ranges_collapsed_df)
 
     
 }
 

 
 coords_last250_df_collapse <- lapply(unique(coords_last250_granges$ensembl_gene_id),
                                                         function(g) collapse_gene(subset(coords_last250_granges, ensembl_gene_id == g)))
   
 coords_last250_collapsed_df <- data.frame(do.call(rbind, coords_last250_df_collapse))
  
 save(coords_last250_df, coords_last250_collapsed_df, file=file.path(file_path_ref$project_directory, "data/external_data/ensembl_GRCm38_101_coords_last250.RData"))
 
 coords_last250_collapsed_gtf <- coords_last250_collapsed_df %>%
   rowwise() %>%
   mutate(source="GRCm38_101_last250_utrs", 
          feature="utr",
          frame=".",
          attribute=paste0(c(sprintf('gene_id "%s"', ensembl_gene_id),
                          sprintf(' transcript_id "%s"', gsub(";", ":",ensembl_transcript_id))), collapse=";"))
 
 write.table(coords_last250_collapsed_gtf[,c("chromosome_name", "source", "feature", "start", "end", "width", "strand", "frame", "attribute")], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, file=file.path(file_path_ref$project_directory, "data/external_data/GRCm38_101_last250bp_collapse.gtf"))
 }



## ---------------------------------------------------------------------------------------------------------------
read_counts_last250 <- read.table(file.path(file_path_ref$project_directory,
                                            file_path_ref$sequencing$results$reads_last250_3utr), 
                                  sep="\t",
                                  header=TRUE)

colnames(read_counts_last250) <- gsub("X.n.holyscratch01.macklis_lab.priyav.lab_notebook.cpn_cth_paper_final.star_align.", "", gsub(".Aligned.sorted.bam|.noUMI.noAdapter.noPolyA.fastq_slamdunk_mapped_filtered.UMIDedup.bam", "", colnames(read_counts_last250)))

# need to fix a few names that were truncated
load(file.path(file_path_ref$project_directory, "data/external_data/ensembl_GRCm38_101_coords_last250.RData"))

read_counts_last250$Geneid_first <- sapply(read_counts_last250$Geneid, 
                                           function(x) strsplit(as.character(x), ":")[[1]][1])

coords_last250_collapsed_df$Geneid_first <- sapply(coords_last250_collapsed_df$ensembl_transcript_id,
                                                  function(x) strsplit(as.character(x), ";")[[1]][1])


read_counts_last250_fix <- read_counts_last250 %>%
    inner_join(unique(coords_last250_collapsed_df[,c("Geneid_first", "ensembl_transcript_id", "ensembl_gene_id")]), by="Geneid_first")

read_counts_meta_cols <- c("ensembl_transcript_id", "ensembl_gene_id", "Geneid_first", "Chr",
                           "Start", "End", "Strand", 
                           "Length")



## ---------------------------------------------------------------------------------------------------------------
counts_cols <- colnames(read_counts_last250_fix)
quantseq_names <- counts_cols[grepl("QuantSeq",counts_cols)]

quantseq_counts <- read_counts_last250_fix[,quantseq_names]

# divide counts by sum of total counts per sample
quantseq_counts_norm <- sweep(quantseq_counts, 2, colSums(quantseq_counts), "/")*1e6

# get the average counts per transcript
quantseq_counts_norm$avg_counts_norm <- rowMeans(quantseq_counts_norm)

# add the metadata
quantseq_counts_norm <- data.frame(cbind(read_counts_last250_fix[,read_counts_meta_cols], 
                                         quantseq_counts_norm))

valid_txpt_df <- subset(quantseq_counts_norm, avg_counts_norm > 2)



# plot to determine expression cutoff
plt <- ggplot(quantseq_counts_norm, aes(avg_counts_norm)) + 
  geom_density() +
  scale_x_log10() +
  theme_classic() + 
  geom_vline(xintercept=2, linetype="dashed") +
  ggtitle(sprintf("QuantSeq Counts Per Million\nFor the Last 250bp of a 3'UTR in CPN\n%s 3'UTRs with CPM > 2", nrow(valid_txpt_df)))

plot(plt)

ggsave(file.path(file_path_ref$project_directory, file_path_ref$plots$plot_directory, "quantseq_supplement_defining_expressed_3utrs.pdf"), plot=plt, height=100, width=100, units="mm")

# now look at number of transcripts per each gene


txpt_per_gene <- quantseq_counts_norm %>% 
  subset(avg_counts_norm > 2) %>%
  group_by(ensembl_gene_id) %>%
  summarize(num_txpt=n()) 

plt <- ggplot(txpt_per_gene, aes(num_txpt)) + 
  geom_bar() +
  theme_classic() +
  ggtitle(sprintf("QuantSeq Number Transcripts per Gene\n%s Genes Total", length(unique(quantseq_counts_norm$ensembl_gene_id))))
  
ggsave(file.path(file_path_ref$project_directory, file_path_ref$plots$plot_directory, "quantseq_supplement_defining_num3prime_ends_per_gene.pdf"), plot=plt, height=100, width=100, units="mm")




## ---------------------------------------------------------------------------------------------------------------
valid_txpt_df$txpt_group <- 1:nrow(valid_txpt_df)

txpts_expanded <- valid_txpt_df %>% 
                              separate_rows(ensembl_transcript_id, sep=";")

txpt_utr_lengths <- getBM(c("ensembl_transcript_id", "chromosome_name", "strand", "3_utr_start", "3_utr_end"),
                          filters="ensembl_transcript_id",
                          values=unique(txpts_expanded$ensembl_transcript_id),
                          mart=grcm38_101) %>%
  drop_na() %>% # required because they always for some reason return NA
  group_by(ensembl_transcript_id) %>%
  summarize(utr_length=sum(abs(`3_utr_end`-`3_utr_start`)),
            utr_coords=paste(sprintf("%s:%s-%s:%s", chromosome_name, `3_utr_start`, `3_utr_end`, ifelse(strand < 0, "-", "+")), collapse=";"))

# retain the longest or one of the longest if multiple are same length
txpt_longest <- txpts_expanded %>% 
  inner_join(txpt_utr_lengths, by="ensembl_transcript_id") %>%
  group_by(txpt_group) %>%
  summarize(longest_txpt=ensembl_transcript_id[which.max(utr_length)[1]],
            utr_length_longest=utr_length[which.max(utr_length)[1]],
            utr_coords_longest=utr_coords[which.max(utr_length)[1]])
  
  
valid_txpt_df <- valid_txpt_df %>% inner_join(txpt_longest, by="txpt_group")

save(valid_txpt_df, file=file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$valid_transcripts_quantseq))

read_counts_last250_utr_annot <- read_counts_last250_fix %>% inner_join(unique(valid_txpt_df[,c("Geneid_first", "txpt_group", "longest_txpt", "utr_length_longest", "utr_coords_longest", "ensembl_gene_id", "ensembl_transcript_id")]))

save(read_counts_last250_utr_annot, file=file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$counts_last250_valid_transcripts_quantseq))


## ---------------------------------------------------------------------------------------------------------------
# must load
# file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$valid_transcripts_quantseq)
# file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$counts_last250_valid_transcripts_quantseq))
counts_data <- data.frame(read_counts_last250_utr_annot %>% 
  group_by(ensembl_transcript_id) %>%
  summarize(across(all_of(metadata$sample_id), sum)))

rownames(counts_data) <- counts_data$ensembl_transcript_id

calc_gc_vs_bkg_last250 <- function(dev_stage, stype, bkg) {
   meta_subs <- subset(metadata, 
                      prep_type == "GC" & 
                        ((stage ==  dev_stage & subtype == stype) | subtype == bkg))
  
  print(meta_subs)
  
  
  # run deseq
  dds <- DESeqDataSetFromMatrix(countData=counts_data[,meta_subs$sample_id],
                                colData=meta_subs,
                                design=~enrichment_type)
  
  dds <- DESeq(dds)
  
  print(resultsNames(dds))
  
  # get the result
  res <- results(dds, name="enrichment_type_signal_vs_background",
                    altHypothesis = "greater")
  
  print(summary(res))
  
  res_df <- data.frame(res)
  res_df$ensembl_transcript_id <- rownames(res_df)
  
  res_df <- res_df %>% 
    rowwise() %>%
    mutate(is_sig=padj < 0.05) %>%
    inner_join(unique(valid_txpt_df[,c("ensembl_transcript_id", "ensembl_gene_id", "utr_length_longest", "longest_txpt")]))
  
  return(res_df)
  
}


calc_gc_vs_gc_last250 <- function(stype, geneids_above_bkg) {
   meta_subs <- subset(metadata, 
                      prep_type == "GC" & 
                        ( subtype == stype))
  
  print(meta_subs)
  
  
  # run deseq
  dds <- DESeqDataSetFromMatrix(countData=counts_data[,meta_subs$sample_id],
                                colData=meta_subs,
                                design=~stage)
  
  dds <- DESeq(dds)
  
  print(resultsNames(dds))
  
  # get the result
  res <- results(dds, name="stage_DPC23_vs_DPC21")
  
  print(summary(res))
  
  res_df <- data.frame(res)
  res_df$ensembl_transcript_id <- rownames(res_df)
  
  res_df <- res_df %>% 
    rowwise() %>%
    mutate(is_sig=padj < 0.05) %>%
    inner_join(unique(valid_txpt_df[,c("ensembl_transcript_id", "ensembl_gene_id", "utr_length_longest", "longest_txpt")])) %>%
    rowwise() %>%
    mutate(is_valid_gc=longest_txpt %in% geneids_above_bkg)
  
  return(res_df)
  
}




## ---------------------------------------------------------------------------------------------------------------
dpc23_cpn_res <- calc_gc_vs_bkg_last250("DPC23", "CPN", "Forebrain_Background")
dpc21_cpn_res <- calc_gc_vs_bkg_last250("DPC21", "CPN", "Forebrain_Background")
dpc23_cpn_res$stype <- "P3_CPN_GC"
dpc21_cpn_res$stype <- "P1_CPN_GC"

stage_cpn_res <- calc_gc_vs_gc_last250("CPN", union(subset(dpc23_cpn_res, is_sig)$longest_txpt,
                                                    subset(dpc21_cpn_res, is_sig)$longest_txpt))

# look at thalamus
dpc23_cth_res1 <- calc_gc_vs_bkg_last250("DPC23", "CThPN", "Forebrain_Background")
dpc23_cth_res2 <- calc_gc_vs_bkg_last250("DPC23", "CThPN", "Thalamus_Background")
dpc23_cth_res <- rbind(dpc23_cth_res1,
                       dpc23_cth_res2) %>%
  drop_na() %>%
  group_by(ensembl_transcript_id, longest_txpt, utr_length_longest) %>%
  summarize(is_sig=sum(is_sig) == 2)
dpc23_cth_res$stype <- "P3_CTh_GC"
dpc23_gc_res <- rbind(dpc23_cth_res, dpc23_cpn_res[,c("ensembl_transcript_id", "utr_length_longest", "is_sig", "stype")]) 
save(dpc23_gc_res, file=file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$dpc23_subtype_utr_analysis))

dpc23_gc_res_grouped <- dpc23_gc_res %>%
  drop_na() %>%
  group_by(ensembl_transcript_id, utr_length_longest) %>%
  summarize(expression_group=paste(stype[is_sig], collapse=";")) %>%
  rowwise() %>%
  mutate(expression_group_readable=ifelse(nchar(expression_group) == 0, "Soma",
                                          ifelse(grepl(";", expression_group), "Shared_GC",
                                          expression_group)))
dpc23_gc_res_grouped$expression_group_readable <- factor(dpc23_gc_res_grouped$expression_group_readable,
                                                    levels=c("Soma", "P3_CPN_GC", "Shared_GC", "P3_CTh_GC"))

ggplot(dpc23_gc_res_grouped, aes(expression_group_readable, utr_length_longest, 
                            fill=expression_group_readable)) + 
  geom_violin(draw_quantiles = 0.5) +
  scale_y_log10()


cpn_gc_res <- rbind(dpc21_cpn_res, dpc23_cpn_res) 
save(cpn_gc_res, file=file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$cpn_stage_utr_analysis))

cpn_res_grouped <- cpn_res %>%
  drop_na() %>%
  group_by(ensembl_transcript_id, longest_txpt, utr_length_longest) %>%
  summarize(expression_group=paste(stype[is_sig], collapse=";")) %>%
  rowwise() %>%
  mutate(expression_group_readable=ifelse(nchar(expression_group) == 0, "Soma",
                                          ifelse(grepl(";", expression_group), "Shared_GC",
                                          expression_group)))
cpn_res_grouped$expression_group_readable <- factor(cpn_res_grouped$expression_group_readable,
                                                    levels=c("Soma", "P1_CPN_GC", "Shared_GC", "P3_CPN_GC"))

mean_comparisons <- list(c("Soma", "P1_CPN_GC"), c("P1_CPN_GC", "Shared_GC"), c("Shared_GC", "P3_CPN_GC"))

plt <- ggplot(cpn_res_grouped, aes(expression_group_readable, utr_length_longest, 
                            fill=expression_group_readable)) + 
  geom_violin(draw_quantiles = 0.5) + 
  theme_classic() + 
  xlab("") +
  ylab("3'UTR Length") +
  ggtitle("UTR Lengths Across Time\nCPN GCs") +
  scale_fill_manual(values=c(background_color,
                             colorRampPalette(c(cpn_color_light, cpn_color_dark))(3))) +
  theme(legend.position="None", text=element_text(size=7)) +
  stat_compare_means(method = "anova", label.y=5.5, size=3)+
  stat_compare_means(aes(label=after_stat(p.signif)), comparisons=mean_comparisons, method="wilcox.test") +
  scale_y_log10()

ggsave(file.path(file_path_ref$project_directory, file_path_ref$plots$plot_directory, "figure5_cpn_utr_lengths_time.pdf"),
       plot=plt, 
       device="pdf",
       height=100, width=80, units="mm")


## ---------------------------------------------------------------------------------------------------------------
load(file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$valid_transcripts_quantseq))
utr_sequences <- getBM(c("3utr", "ensembl_transcript_id"),
                       filters="ensembl_transcript_id",
                       values=unique(valid_txpt_df$longest_txpt),
                       mart=grcm38_101)

utr_sequences_valid <- subset(utr_sequences, !grepl("N", `3utr`))
utr_lst <- as.list(utr_sequences_valid$`3utr`)
names(utr_lst) <- utr_sequences_valid$ensembl_transcript_id

# get ordered cpn
ordered_cpn_gc <- stage_cpn_res %>% subset(is_valid_gc & utr_length_longest > 20) %>% arrange(desc(log2FoldChange))
ordered_cpn_gc_utrs <- unlist(utr_lst[ordered_cpn_gc$longest_txpt])
spma_cpn_gc <- run_matrix_spma(ordered_cpn_gc_utrs, n_bins=7, cache=FALSE)

ordered_cpn_gc_lst <- ordered_cpn_gc$log2FoldChange
names(ordered_cpn_gc_lst) <- ordered_cpn_gc$ensembl_gene_id
cpn_gc_gsea <- gseGO(ordered_cpn_gc_lst, ont="ALL", OrgDb = org.Mm.eg.db, keyType="ENSEMBL")

save(spma_cpn_gc, cpn_gc_gsea, file=file.path(file_path_ref$project_directory, "data/sequencing/differential_expression/cpn_stage_isoform_level_motif_and_geneset_enrichment.RData"))

plt <- ridgeplot(simplify(cpn_gc_gsea), showCategory=20, label_format = 40)
ggsave(file.path(file_path_ref$project_directory, file_path_ref$plots$plot_directory, "cpn_stage_isoform_level_geneset_enrichment.pdf"), plot=plt, height=200, width=200, units="mm")

cpn_soma_utr_scores <- transite::score_transcripts(utr_lst, cache=FALSE)

cpn_gc_utr_scores <- transite::score_transcripts(unlist(utr_lst[subset(cpn_res_grouped, grepl("GC", expression_group_readable))$longest_txpt]), cache=FALSE)


cpn_gc_dpc23_txpts <- subset(dpc23_cpn_res, is_sig)$longest_txpt #subset(cpn_res_grouped, expression_group_readable == "P3_CPN_GC")$longest_txpt
cpn_gc_dpc23_scores <- transite::score_transcripts(unlist(utr_lst[cpn_gc_dpc23_txpts]), cache=FALSE)

cpn_gc_dpc21_txpts <- subset(dpc21_cpn_res, is_sig)$longest_txpt #subset(cpn_res_grouped, expression_group_readabe == "P1_CPN_GC")$longest_txpt
cpn_gc_dpc21_scores <- transite::score_transcripts(unlist(utr_lst[cpn_gc_dpc21_txpts]), cache=FALSE)

#cpn_gc_shared_txpts <- subset(cpn_res_grouped, expression_group_readable == "Shared_GC")$longest_txpt
#cpn_gc_shared_scores <- transite::score_transcripts(unlist(utr_lst[cpn_gc_shared_txpts]), cache=FALSE)

cpn_gc_dpc23_enrich <- transite::calculate_motif_enrichment(cpn_gc_dpc23_scores$df, 
                                                            cpn_gc_utr_scores$df,
                                                            cpn_gc_utr_scores$total_sites,
                                                            cpn_gc_utr_scores$absolute_hits,
                                                            length(cpn_gc_dpc23_txpts))


cpn_gc_dpc21_enrich <- transite::calculate_motif_enrichment(cpn_gc_dpc21_scores$df, 
                                                            cpn_gc_utr_scores$df,
                                                            cpn_gc_utr_scores$total_sites,
                                                            cpn_gc_utr_scores$absolute_hits,
                                                            length(cpn_gc_dpc21_txpts))


cpn_gc_shared_enrich <- transite::calculate_motif_enrichment(cpn_gc_shared_scores$df, 
                                                            cpn_gc_utr_scores$df,
                                                            cpn_gc_utr_scores$total_sites,
                                                            cpn_gc_utr_scores$absolute_hits,
                                                            length(cpn_gc_shared_txpts))



## ---------------------------------------------------------------------------------------------------------------
load(file.path(file_path_ref$project_directory, file_path_ref$sequencing$results$cpn_soma_comparison))

rbps  <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0003729", columns="ENSEMBL")

dge_res_rbp <- dge_res %>% rowwise() %>% mutate(is_mrna_rbp=ensembl_gene_id %in% unique(rbps$ENSEMBL))

ggplot(dge_res_rbp %>% subset(is_mrna_rbp), aes(baseMean, log2FoldChange, color=padj < 0.05)) + 
  geom_point() + 
  geom_text_repel(data=dge_res_rbp %>% subset(is_mrna_rbp & abs(log2FoldChange) > 0.8), 
                   mapping=aes(baseMean, log2FoldChange, label=external_gene_name)) +
  scale_x_log10() + scale_color_manual(values=c(background_color, cpn_color)) +
  theme_classic() +
  geom_hline(yintercept=0) +
  ylab("P1 CPN <- Log2 FoldChange -> P3 CPN") +
  theme(legend.position="None")

