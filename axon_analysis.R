## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------
library(MASS)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggupset)
library(ggrepel)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(rjson)
library(pheatmap)

source("Y:/users/priyav/slam_final/scripts/helper_script.R")

file_path_ref <- fromJSON(file = input)

data_subdir <- file.path(file_path_ref$project_directory, "data/SLAM_Axon_Analysis")
plot_subdir <- file.path(file_path_ref$project_directory, "plots/SLAM_Axon_Analysis")
dir.create(data_subdir, recursive = TRUE)
dir.create(plot_subdir, recursive = TRUE)

# filter out the problematic samples R5 12h Dist --> seem to be the same?
metadata <- metadata %>% subset(!grepl("R5_UPRT_12h_Dist|R5_Neg_12h_Dist", sample_id) & location %in% c("Distal_Axon", "Proximal_Axon")) 

save_pdf <- FALSE # whether to save plots as pdf
overwrite <- FALSE # whether to re-run analysis or load saved analysis


## ---------------------------------------------------------------------------------------------------------------

# fix a gene identifier from HTseq in wihc some are __ambiguous[ensembl1+ensembl2]
fix_gene_id <- function(gene_id) {
  return(strsplit(gsub("__ambiguous|\\[|\\]", "", gene_id), "\\+")[[1]])
}




## ---------------------------------------------------------------------------------------------------------------

get_all_gene_conversions <- function(axon_slam_dir=file.path(data_subdir, 
                                                             "axon_conversion_counts"),
                                     overwrite=FALSE) {
  
  savefile <- file.path(axon_slam_dir, "axon_conversion_counts_agg.RData")

  if (!file.exists(savefile) | overwrite) {
   meta_subset <- metadata
  axon_slam_files <- file.path(axon_slam_dir, sprintf("%s_.no_intron.tc_leq2.counts.csv", meta_subset$sample_id ))
  
  axon_slam_files <- axon_slam_files[file.exists(axon_slam_files)]
  print(length(axon_slam_files))
  conversion_df <- data.frame(do.call(rbind, lapply(axon_slam_files,
                 function(s) { 
                    df <- read.csv(s, header=FALSE)
                    colnames(df) <- c("GF", "nT", "nTC", "reads", "TC_reads")
                    df$sample_id <- gsub("_.no_intron.tc_leq2.counts.csv", "", basename(s))

                   return(df)
}
                )))
  
  save(conversion_df, file=savefile)
  
  # also write to csv
  write.csv(conversion_df, file=file.path(counts_dir, "axon_conversion_counts_agg.csv"))
  
  return(conversion_df)

  } else {
  load(savefile)
    return(conversion_df)
}
}

conversion_df <- get_all_gene_conversions()



total_counts <- conversion_df %>% group_by(sample_id) %>% summarize(total_reads=sum(reads)) 

total_counts %>% ggplot(aes(sample_id, total_reads)) + geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept=30e6)




## ---------------------------------------------------------------------------------------------------------------
read_save_intron_df <- function(loc, sample_time, 
                                axon_slam_dir=file.path(data_subdir, 
                                                             "axon_conversion_counts",
                                                        overwrite=FALSE)) {
  
  savefile <- file.path(axon_slam_dir, sprintf("%s_%s_intron_reads_aggregated.RData", loc, sample_time))

  if (!file.exists(savefile) | overwrite) {
  meta_subset <- subset(metadata, location == loc & sample_time_min == sample_time)
  axon_slam_files <- file.path(axon_slam_dir, sprintf("%s_counts.csv.gz", meta_subset$sample_id ))
  axon_slam_files <- axon_slam_files[file.exists(axon_slam_files)]

intron_slam_df <- data.frame(do.call(rbind, lapply(axon_slam_files, 
                 function(s) { df <- data.frame(fread(s))
                   
                   df_subs <- df %>% 
                   subset(ai == 1) %>% 
                   group_by(GF) %>% 
                   summarize(nT=sum(as.numeric(nref)),
                             nTC=sum(as.numeric(TC)),
                             nreads=n())
                 
                 df_subs$sample_id <-gsub("_counts.csv.gz", "", basename(s))
                 return(df_subs)
}
                )))
print(head(intron_slam_df))

save(intron_slam_df, file=savefile)

return(intron_slam_df)
  } else {
  load(savefile)
    return(intron_slam_df)
  }
}

get_all_intron_conversions <- function(overwrite=FALSE) {
    savefile <- file.path(data_subdir, 
                          "axon_conversion_counts", 
                          "all_genes_intron_reads_aggregated.RData")
    
    if (!file.exists(savefile) | overwrite) {
      
intron_df1 <- read_save_intron_df("Distal_Axon", 720, overwrite=TRUE) 

intron_df2 <- read_save_intron_df("Distal_Axon", 1080) 

intron_df3 <- read_save_intron_df("Proximal_Axon", 720)

intron_df4 <- read_save_intron_df("Proximal_Axon", 1080)

intron_df <- data.frame(do.call(rbind, list(intron_df1, intron_df2, intron_df3, intron_df4)))
save(intron_df, file=savefile)
return(intron_df)
      
    } else {
      load(savefile)
      return(intron_df)
    }
    

}

intron_conversion_df <- get_all_intron_conversions()


intron_median_rates <- intron_conversion_df %>% 
  inner_join(metadata, by="sample_id") %>%
  subset(nreads > 10) %>%
  group_by(sample_id) %>%
  summarize(median_intron_convrate=median(nTC/nT, na.rm=TRUE))

if (overwrite) {
  write.csv(intron_median_rates, file=file.path(data_subdir, 
                                              "axon_conversion_counts", 
                                              "axon_intron_median_rates.csv"))
}


# confirm that median intron rates are the same across biological samples
intron_median_rates %>% inner_join(metadata, by="sample_id") %>% rowwise() %>% mutate(bio_sample=strsplit(sample_id, "_")[[1]][1]) %>% pivot_wider(id_cols=c(bio_sample, IUE_plasmid, sample_time_min), names_from=location, values_from=median_intron_convrate) %>% ggplot(aes(Proximal_Axon, Distal_Axon, color=IUE_plasmid)) + geom_point() + stat_cor() + theme_classic() + scale_color_manual(values=c(background_colors[4], accent_colors[2])) + theme(legend.position = "None") + ggtitle("Median Intron Conversion Rates\nPer Biological Sample")



### join all data 
conversion_data <- conversion_df %>% inner_join(metadata, by="sample_id") %>% inner_join(intron_median_rates, by="sample_id") 


# calculate sample conversion rates and compare to intron convrate
sample_conv_rates <- conversion_data %>%
  subset(!grepl("R5_UPRT_12h_Dist|R5_Neg_12h_Dist", sample_id)) %>%
  group_by(sample_id, IUE_plasmid) %>%
  summarize(sample_convrate=sum(nTC, na.rm = TRUE)/sum(nT, na.rm = TRUE),
            median_sample_convrate=median(nTC/nT, na.rm=TRUE))

plt <- sample_conv_rates %>%
  inner_join(intron_median_rates, by="sample_id") %>%
  ggplot(aes(median_sample_convrate, median_intron_convrate, color=IUE_plasmid)) +
  geom_point(size=0.5) + 
  stat_cor(method="spearman") +
  theme_classic() +
  scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
  theme(text=element_text(size=7)) 
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_supplement_intron_convrate_vs_sample_convrate.pdf"),
       plot=plt,
       height=40, width=70, units="mm")
}


## ---------------------------------------------------------------------------------------------------------------
# find counts per sample
total_sample_counts <- conversion_data %>%
  group_by(sample_id) %>%
  summarize(read_total=sum(reads),
            T_total=sum(nT))

# get reads CPM and TC reads CPM
all_convs <- conversion_data %>%
  inner_join(total_sample_counts, by="sample_id") %>%
  rowwise() %>%
  mutate(reads_CPM=reads/read_total*1e6,
         TC_reads_CPM=TC_reads/read_total*1e6)

# find number of samples per group for thresholding proportion of samples 
n_samples_per_group <- all_convs %>% 
  dplyr::select(sample_id, IUE_plasmid, location, sample_time_min) %>% 
  unique() %>% 
  group_by(IUE_plasmid, location, sample_time_min) %>% 
  summarize(nsamples=n())  
  

get_valid_gene_df <- function(TC_reads_CPM_thresh, 
                              convs_df=all_convs) {
       
       n_samples_per_group <- all_convs %>% 
          dplyr::select(sample_id, IUE_plasmid, location, sample_time_min) %>% 
          unique() %>% 
          group_by(IUE_plasmid, location, sample_time_min) %>% summarize(nsamples=n())    

        
    valid_gene_df <- convs_df %>% subset(!grepl("__alignment_not_unique|__no_feature", GF)) %>% 
      inner_join(n_samples_per_group, by=c("location", "sample_time_min", "IUE_plasmid")) %>%
      rowwise() %>% mutate(is_tc_valid_sample=TC_reads_CPM > TC_reads_CPM_thresh & TC_reads > 10,
                          is_t_valid_sample=reads > 50) %>%
      group_by(GF, location, sample_time_min, IUE_plasmid) %>%
      summarize(is_t_valid_group=sum(is_t_valid_sample) >= nsamples[1]*3/4,
                is_tc_valid_group=sum(is_tc_valid_sample) >= nsamples[1]*3/4) %>%
      group_by(GF, location, sample_time_min) %>%
      summarize(is_valid=(sum(is_t_valid_group) == 2),
                group=paste0(IUE_plasmid[is_tc_valid_group], collapse=";")) 
    
    return(valid_gene_df)
    }
    


    
    thresh_df <- data.frame(do.call(rbind, lapply(0:20, function(i) get_valid_gene_df(i) %>% subset(is_valid) %>% group_by(location, sample_time_min, group) %>% summarize(num_genes=n(), thresh=i))))
 
plt <- ggplot(thresh_df %>% subset(nchar(group) > 0), aes(factor(thresh), num_genes, color=group)) + 
  geom_point(size=0.5) +
  facet_grid(rows=vars(sample_time_min), cols=vars(location)) +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=7), legend.position="None") +
  scale_color_manual(values=c(accent_colors[10], background_colors[5], accent_colors[2]))
plot(plt)
if (save_pdf) {
ggsave(file.path(plot_subdir, "axon_supplement_tc_thresh_num_genes_pass_prox_dist.pdf"),
       plot=plt,
       height=60, width=60, units="mm")
}

# define the valid genes

# write genes that pass thresh 2 to file
valid_genes_savefile <- file.path(data_subdir, "axon_conversion_counts/axon_valid_genes.csv")
if (overwrite) {
  valid_gene_df <- get_valid_gene_df(2)
  write.csv(valid_gene_df, file=valid_genes_savefile)
} else {
  valid_gene_df <- read.csv(valid_genes_savefile)
}


## ---------------------------------------------------------------------------------------------------------------
# function to fit the GLM on given gene counts and design
# Arguments
#   gene_subs: data frame with columns IUE_plasmid (must be a character vector of "pPV25" and "dollar_51" ), tc_scaled, 
#   fit_name: character, what to name the model parameter estimates
#   design: character representing design; usually "tc_scaled~IUE_plasmid"
# Returns:
#    data frame with model results including coefficient estimates, coefficient z-values, coefficient p-values, 
#    AIC, 2*log likelihood, theta estimate for neg. binomial model, deviance and null deviance, 
#    mean squared error and mean absolute error of model prediction vs truth, r-squared goodness of fit
run_fit_glm_nb <- function(gene_subs, fit_name, design, return_model=FALSE) {
  tryCatch({
    fit3 <- glm.nb(as.formula(design), data=gene_subs)
    res_summary <- summary(fit3)
    
    coefs <- data.frame(res_summary$coefficients)
    colnames(coefs) <- c("estimate", "stderr", "z_value", "p_value")
    coefs$coef_name <- gsub("\\(|\\)", "", rownames(coefs))
    coefs_wide <- data.frame(coefs %>% pivot_wider(names_from=coef_name, values_from=c(estimate, stderr, z_value, p_value)))
    fit_info <- data.frame(deviance=res_summary$deviance,
                           theta=res_summary$theta,
                           theta_se=res_summary$SE.theta,
                           aic=res_summary$aic,
                           twologlik=res_summary$twologlik,
                           null_deviance=res_summary$null.deviance
    )
    
    # calculate error metric
    abs_err <- abs(predict(fit3, newdata=gene_subs)-log(gene_subs$tc_scaled))
    sq_err <- (predict(fit3, newdata=gene_subs)-log(gene_subs$tc_scaled))^2
    
    # for rsquared
    ss_tot <- sum((log(gene_subs$tc_scaled)-mean(log(gene_subs$tc_scaled)))^2)
    
    fit_info$train_mean_abs_err <- mean(abs_err)
    fit_info$train_mean_sq_err <- mean(sq_err)
    fit_info$train_rsquared <- 1-sum(sq_err)/ss_tot
    
    return_df <- data.frame(cbind(coefs_wide, fit_info))
    return_df$glm_nb_fit_name <- fit_name
    
    if (return_model) {
      return(list(res_df=return_df, model=fit3))
    }
    return(return_df)
  }, error=function(e){ 
    
    return_df <- data.frame( estimate_Intercept=NA,
                             estimate_IUE_plasmidpPV25=NA,
                             stderr_Intercept=NA,
                             stderr_IUE_plasmidpPV25=NA,
                             z_value_Intercept=NA,
                             z_value_IUE_plasmidpPV25=NA,
                             p_value_Intercept=NA,
                             p_value_IUE_plasmidpPV25=NA,
                             deviance=NA,
                             theta=NA,
                             theta_se=NA,
                             aic=NA,
                             twologlik=NA,
                             null_deviance=NA,
                             train_mean_abs_err=NA,
                             train_mean_sq_err=NA,
                             train_rsquared=NA,
                             glm_nb_fit_name=fit_name)
    
    if (return_model) {
      return(list(res_df=return_df, model=NULL))
    }
        return(return_df)
  })
}

# function that wraps run_fit_glm_nb to run different variants of the neg. binom glm
# Arguments:
#   geneid: character geneid must exist in input_data GF column
#   loc: location, one of "Proximal_Axon" or "Distal_Axon"
#   stime: integer, one of 720 or 1080
#   design: string design for neg binom model
#   input_data: input data frame including columns GF, location, sample_time_min, nT, nTC, and a column for intron scaling
#   intron_scale_column: character name of the column with intron scaling info
run_t_intron_scaled_glm <- function(geneid, loc, stime, 
                                  design="tc_scaled~IUE_plasmid", 
                                  input_data=conv_info, 
                                  intron_scale_column="median_intron_convrate",
                                  run_type="all") {
  
  gene_subs <- subset(input_data, GF == geneid & location == loc & sample_time_min == stime)
  t_conv_scale <- median(gene_subs$nT)/gene_subs$nT
  intron_scale <- median(unlist(gene_subs[intron_scale_column]))/unlist(gene_subs[intron_scale_column])
  
  gene_subs$tc_t_scaled=round(gene_subs$nTC*t_conv_scale)
  gene_subs$tc_scaled=round(gene_subs$tc_t_scaled*intron_scale)
  
  if (run_type == "all") {
    fit_lst[[1]] <- run_fit_glm_nb(gene_subs, "all", design)
    
    # leave one out cross valiation 
  } else if (run_type == "cv") {
    
    
    # now leave one out xval
    for (i in 1:nrow(gene_subs)) {
      keep_idx <- setdiff(1:nrow(gene_subs), i)
      out_sample <- gene_subs$sample_id[i]
      glm_res <- run_fit_glm_nb(gene_subs[keep_idx, ], sprintf("%s_excluded", out_sample), design, return_model=TRUE)
      glm_res_df <- glm_res$res_df
      
      
      if (!is.null(glm_res$model)) {
        prediction <- predict(glm_res$model, newdata=data.frame(IUE_plasmid=gene_subs$IUE_plasmid[i]))
        truth <- log(gene_subs$tc_intron_scaled[i])
        glm_res_df$prediction <- prediction
        glm_res_df$truth <- truth
      } else {
        glm_res_df$prediction <- NA
        glm_res_df$truth <- NA
      }
      
      fit_lst[[(i)]] <- glm_res_df
      
      
    }
    
    # bootstrapping
  } else if (run_type == "bootstrap") {
    
    set.seed(67) # important to set the seed so that all genes have same random ordering
    lst_idx <- 1 #nrow(gene_subs) + 1
    for (j in 1:100) {
      random_idx <- sample(1:nrow(gene_subs), size=nrow(gene_subs), replace=TRUE) # bootstrap
      #subs_test <- gene_subs
      subs_test$IUE_plasmid <- gene_subs$IUE_plasmid[random_idx]
      random_order <- gene_subs$IUE_plasmid[random_idx]
      #subs_test$tc_intron_scaled <- gene_subs$tc_intron_scaled[random_idx] # this was for randomize_values 
      if (paste0(random_order, collapse = ";") != paste0(gene_subs$IUE_plasmid, collapse=";")) {
        fit_lst[[lst_idx]] <- run_fit_glm_nb(subs_test,
                                             sprintf("random_order_%s", paste0(random_idx, collapse = ";")),
                                             design)
        lst_idx <- lst_idx + 1
      } 
    }
  }
  
  combined_fits <- data.frame(
    do.call(rbind, fit_lst))
  
  combined_fits$gene_id <- geneid
  combined_fits$sample_time_min <- stime
  combined_fits$location <- loc
  
  return(combined_fits)
}

# Function to adjust p-values for each condition in the data
# requres a df with sample time, and loc facets
# and also pvalue column on which to do adjustment
adjust_p_values_condition <- function(df){
  res_lst <- list()
  for (l in c("Distal_Axon", "Proximal_Axon")) {
    for (s in c(720, 1080)) {
      stype <- sprintf("%s_%s", l, s)
      stype_res <- subset(df, sample_time_min == s & location == l)
      stype_res$padj <- p.adjust(stype_res$pvalue, method="BH")
      stype_res$sample_type <- stype
      res_lst[[stype]] <- stype_res
    }
  }
  
  return(data.frame(do.call(rbind, res_lst)))
}


full_model_savefile <- file.path(data_subdir, "glm_results", "axon_neg_binomial_GLM_results.full_model.RData")

if (!file.exists(full_model_savefile) | overwrite) {
  
  # run the glm for all rows
  all_glm_res <- list()
  for (i in 1:nrow(valid_gene_df)) {
    if (i %% 50 == 0) {
      print(i)
    }
    g <- valid_gene_df$GF[i]
    l <- valid_gene_df$location[i]
    s <- valid_gene_df$sample_time_min[i]
    all_glm_res[[i]] <- run_t_intron_scaled_glm(g, l, s, "tc_intron_scaled~IUE_plasmid")
    
  }
  all_glm_res_df <- data.frame(do.call(rbind, all_glm_res))
  
  # define a one-sided (upper) p-value
  all_glm_res_df$pvalue <- 1-pnorm(all_glm_res_df$z_value_IUE_plasmidpPV25)
  
  # fdr correction
  all_glm_res_df_padj <- adjust_p_values_condition(all_glm_res_df)
  
  # save results
  save(all_glm_res_df, all_glm_res_df_padj, 
       file=full_model_savefile)
} else {
  load(full_model_savefile)
}


## ---------------------------------------------------------------------------------------------------------------
# reduced_model_savefile <- file.path(data_subdir, "glm_results", "axon_neg_binomial_GLM_results.reduced_model.RData")
# 
# if (!file.exists(reduced_model_savefile) | overwrite) {
#   all_glm_res_reduced <- list()
#   for (i in 1:nrow(valid_gene_df)) {
#     if (i %% 50 == 0) {
#       print(i)
#     }
#     g <- valid_gene_df$GF[i]
#     l <- valid_gene_df$location[i]
#     s <- valid_gene_df$sample_time_min[i]
#     all_glm_res_reduced[[i]] <- run_t_intron_scaled_glm(g, l, s, "tc_intron_scaled~1")
#     
#   }
#   all_glm_res_reduced_df <- data.frame(do.call(rbind, all_glm_res_reduced))
#   
#   save(all_glm_res_reduced_df, file=reduced_model_savefile)
# } else {
#   load(reduced_model_savefile)
# }
# 
# # combine and compare 
# # AIC and Likelihood ratio test
# all_glm_compare <- all_glm_res_df_padj %>% 
#   inner_join(all_glm_res_reduced_df, by=c("sample_time_min", "location", "gene_id"),
#              suffix=c(".full", ".reduced")) %>%
#   rowwise() %>%
#   mutate(lr_stat=-twologlik.reduced+twologlik.full) %>%
#   rowwise() %>%
#   mutate(lr_pvalue=pchisq(lr_stat, 1, lower.tail=FALSE))
# 
# # Look at AIC
# plt <- all_glm_compare %>%
#   arrange(desc(padj)) %>%
#   rowwise() %>%
#   mutate(abs_diff=log(abs(aic.full-aic.reduced)),
#          sign=ifelse(aic.full-aic.reduced > 0, 1, -1)) %>%
#   rowwise() %>%
#   mutate(xvalue=abs_diff*sign) %>%
#   ggplot(aes(xvalue, color=padj < 0.2)) +
#   stat_ecdf(linewidth=0.5) + 
#   theme_classic() +
#   scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
#   facet_grid(rows=vars(sample_time_min), cols=vars(location)) +
#   xlab("AIC Reduced - AIC Full\nLog Absolute Difference") + 
#   geom_vline(xintercept=0, linewidth=0.2, linetype='dashed', color=background_colors[8]) +
#   theme(text=element_text(size=5))
# plot(plt)
# if (save_pdf) {
# ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_analysis_aic_uprt_vs_intercept_only_model_cdf.pdf"),
#        plot=plt,
#        height=60, width=80, units="mm")
# }
# 
# # look at Likelihood ratio test
# all_glm_compare_lr_padj <- adjust_p_values_condition(
#   all_glm_compare %>% dplyr::select(c(location, sample_time_min, gene_id, lr_pvalue)) %>%
#     dplyr::rename(pvalue=lr_pvalue))
# 
# plt <- all_glm_compare %>%
#   inner_join(all_glm_compare_lr_padj, by=c("sample_time_min", "location", "gene_id"),
#              suffix=c(".wald", ".lrt")) %>%
#   ggplot(aes(padj.wald, padj.lrt, color=padj.wald < 0.2)) + 
#   geom_point(size=0.1) +
#   facet_grid(rows=vars(sample_time_min), cols=vars(location)) +
#   scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
#   xlab("FDR Adjusted p-value:\nWald Test (one-sided)") +
#   ylab("FDR Adjusted p-value:\nLikelihood Ratio Test") +
#   theme_classic() +
#   theme(text=element_text(size=5))
# 
# plot(plt)
# if (save_pdf) {
# ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_analysis_Wald_vs_LRT_padj.pdf"),
#        plot=plt,
#        height=60, width=80, units="mm")
# }




## ---------------------------------------------------------------------------------------------------------------
loo_cv_savefile <- file.path(data_subdir, "glm_results/axon_cross_val.RData")

if (!file.exists(loo_cv_savefile) | overwrite) {
  cv_glm_res <- list()
  for (i in 1:nrow(valid_gene_df)) {
    if (i %% 50 == 0) {
      print(i)
    }
    g <- valid_gene_df$GF[i]
    l <- valid_gene_df$location[i]
    s <- valid_gene_df$sample_time_min[i]
    cv_glm_res[[i]] <- run_t_intron_scaled_glm(g, l, s, "tc_intron_scaled~IUE_plasmid", run_type="cv")
  }
  cv_glm_res_df <- data.frame(do.call(rbind, cv_glm_res))
  cv_glm_res_df <- cv_glm_res_df %>%
    rowwise() %>% 
    mutate(sample_id=gsub("_excluded", "", glm_nb_fit_name)) %>%
    inner_join(metadata, by=c("sample_id", "sample_time_min", "location"))
  
  # also run CV on reduced model
  cv_glm_res_reduced <- list()
  for (i in 1:nrow(valid_gene_df)) {
    if (i %% 50 == 0) {
      print(i)
    }
    g <- valid_gene_df$GF[i]
    l <- valid_gene_df$location[i]
    s <- valid_gene_df$sample_time_min[i]
    cv_glm_res_reduced[[i]] <- run_t_intron_scaled_glm(g, l, s, "tc_intron_scaled~1", run_type="cv")
  }
  cv_glm_res_reduced_df <- data.frame(do.call(rbind, cv_glm_res_reduced)) 
  
  cv_glm_res_reduced_df <- cv_glm_res_reduced_df %>%
    rowwise() %>% 
    mutate(sample_id=gsub("_excluded", "", glm_nb_fit_name)) %>%
    inner_join(metadata, by=c("sample_id", "sample_time_min", "location"))
  
  save(cv_glm_res_df, cv_glm_res_reduced_df, file=loo_cv_savefile)
} else {
  load(loo_cv_savefile)
}

# aggregate the result per gene, per location, per timepoint
# obtain squared error for validation set
# obtain absolute error for validation set
agg_cv_res <- function(cv_res_df) {
  agg_res <- cv_res_df %>% 
    group_by(sample_time_min, location, gene_id) %>%
    summarize(ss_resid=sum((prediction-truth)^2),
              ss_tot=sum((truth-mean(truth))^2),
              mean_express=mean(truth)) %>%
    rowwise() %>%
    mutate(valid_rsquared=1-ss_resid/ss_tot)
  
  return(agg_res)
    
}


cv_glm_res_df_rsq <- agg_cv_res(cv_glm_res_df)



# plot the training set r-squared vs leave one out cross validation r-squared
plt <- cv_glm_res_df %>%
  group_by(sample_time_min, location, gene_id) %>%
  summarize(mean_train_rsquared=mean(train_rsquared)) %>%
  inner_join(cv_glm_res_df_rsq, by=c("sample_time_min", "location", "gene_id")) %>%
  inner_join(all_glm_res_df_padj, by=c("sample_time_min", "location", "gene_id")) %>%
  arrange(desc(padj))  %>%
  ggplot(aes(mean_train_rsquared, valid_rsquared, color=padj < 0.2)) + 
  geom_point(size=0.1) + 
  geom_density_2d(color=background_colors[1]) +
  scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
  facet_grid(rows=vars(sample_time_min), cols=vars(location)) +
  geom_abline(linetype='dashed', linewidth=0.2, color=background_colors[8]) +
  theme_classic() + 
  theme(text=element_text(size=5))

plot(plt)

if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_analysis_rsquared_train_vs_validation.pdf"),
       plot=plt,
       height=60, width=80, units="mm")
} 


## ---------------------------------------------------------------------------------------------------------------
bootstrap_savefile <- file.path(data_subdir, "glm_results/axon_bootstrap_results.RData")

if (!file.exists(bootstrap_savefile) | overwrite) {
  boot_glm_res <- list()
  for (i in 1:nrow(valid_gene_df)) {
    if (i %% 50 == 0) {
      print(i)
    }
    g <- valid_gene_df$GF[i]
    l <- valid_gene_df$location[i]
    s <- valid_gene_df$sample_time_min[i]
    boot_glm_res[[i]] <- run_t_intron_scaled_glm(g, l, s, "tc_intron_scaled~IUE_plasmid", run_type="bootstrap")
  }
  
  boot_glm_res_df <- data.frame(do.call(rbind, boot_glm_res))
  
  save(boot_glm_res_df,  file=bootstrap_savefile)

} else {
  load(bootstrap_savefile)
}

# get aggregated stats from the bootstrap df
# compare "all" (properly ordered) to "random" (bootstrapped)
valid_cols <- c("sample_time_min", "location", 
                "gene_id", "glm_nb_fit_name", "estimate_IUE_plasmidpPV25")
all_vs_boot_df <- rbind(boot_glm_res_df[, valid_cols], all_glm_res_df[, valid_cols])

# get bootstrap pvals and difference between means
bootstrap_res <- all_vs_boot_df %>%
  group_by(sample_time_min, location, gene_id) %>%
  summarize(pvalue=sum(estimate_IUE_plasmidpPV25[grepl("random", glm_nb_fit_name)] >
                                 estimate_IUE_plasmidpPV25[glm_nb_fit_name == "all"],
                               na.rm = TRUE)/length(na.omit(estimate_IUE_plasmidpPV25[grepl("random", glm_nb_fit_name)])),
            bootstrap_mean=mean(estimate_IUE_plasmidpPV25[grepl("random", glm_nb_fit_name)], na.rm=TRUE),
            bootstrap_sd=sd(estimate_IUE_plasmidpPV25[grepl("random", glm_nb_fit_name)], na.rm=TRUE),
            n_boot=sum(grepl("random", glm_nb_fit_name)),
            estimate_true=estimate_IUE_plasmidpPV25[glm_nb_fit_name == "all"]) %>%
  rowwise() %>%
  mutate(true_vs_boot=estimate_true-bootstrap_mean) 

# add in fdr information for coloring purposes
bootstrap_res <- bootstrap_res %>%
  inner_join(all_glm_res_df_padj[,c("padj", valid_cols)], by=c("sample_time_min", "location", "gene_id"), 
             suffix=c(".boot", ".true"))


plt <- ggplot(bootstrap_res %>% arrange(desc(padj)), aes(true_vs_boot, color=(padj < 0.2))) +
  geom_density() +
  theme_classic() +
  #facet_wrap(vars(sample_type)) +
  facet_grid(rows=vars(sample_time_min), cols=vars(location)) +
  theme(legend.position="none", text=element_text(size=5)) +
  scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
  xlab("UPRT regression coefficient:\nAll sample estimate minus mean estimate from 100 bootstraps") 
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/uprt_regression_coeff_true_minus_boot100_colored_by_is_sig_true_and_cv.pdf"),
       plot=plt,
       height=50, width=50, units="mm")
}



## ---------------------------------------------------------------------------------------------------------------
plot_gene <- function(gene_id, input_df) {
  gene_subs <- subset(input_df, gene_id_fix == gene_id) %>%
    rowwise() %>%
    mutate(sample_group_type=sprintf("%s_%s", location, sample_time_min))
  print(gene_subs$sample_group_type)

  
  median_vals <- gene_subs %>%
    group_by(sample_group_type) %>%
    summarize(median_T=median(nT),
              median_intron=median(median_intron_convrate))
  
  gene_subs <- gene_subs %>%
    inner_join(median_vals, by="sample_group_type") %>%
    rowwise() %>%
    mutate(tc_intron_scaled=nTC*median_T/nT*median_intron/median_intron_convrate)
  
  gene_subs$sample_group_type <- factor(as.character(gene_subs$sample_group_type),
                                       levels=c("Proximal_Axon_720", "Distal_Axon_720", "Proximal_Axon_1080", "Distal_Axon_1080"))

  plt <- ggplot(gene_subs, aes(sample_group_type, tc_intron_scaled, fill=IUE_plasmid)) +
    geom_boxplot(linewidth=0.1, outlier.colour = NA) +
    geom_point(size=0.2, position=position_dodge(0.8)) +
    theme_classic() +
    theme(
          text=element_text(size=5), legend.position="None") +
    scale_fill_manual(values=c(background_colors[5], accent_colors[1], "green", "orange")) +
    facet_wrap(vars(GF), scales="free_y") +
    ylab("T->C Conversions Scaled") + 
    xlab("")
  
  return(plt)
}


input_df <- all_convs %>%
  rowwise() %>%
  mutate(gene_id_fix=paste0(fix_gene_id(GF), collapse=";")) %>%
  separate_rows(gene_id_fix, sep=";")


### Gene plotting

# cstb
plt <- plot_gene("ENSMUSG00000005054", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_cstb_boxplots.pdf"), 
       plot=plt,
       height=50, width=50, units="mm")
}

# rhoa
plt <- plot_gene("ENSMUSG00000007815", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_rhoa_boxplots.pdf"), 
       plot=plt,
       height=50, width=150, units="mm")
}

plt <- plot_gene("ENSMUSG00000018428", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_akap1_boxplots.pdf"), 
       plot=plt,
       height=50, width=100, units="mm")
}

# Cacna1e gene
plt <- plot_gene("ENSMUSG00000004110", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_cacna1e_boxplots.pdf"), 
       plot=plt,
       height=50, width=50, units="mm")
}

plt <- plot_gene("ENSMUSG00000021540", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_smad5_boxplots.pdf"), 
       plot=plt,
       height=50, width=150, units="mm")
}

plt <- plot_gene("ENSMUSG00000015002", input_df)
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/axon_figure_efr3a_boxplots.pdf"), 
       plot=plt,
       height=50, width=150, units="mm")
}


## ---------------------------------------------------------------------------------------------------------------
# upset plot
plt <- all_glm_res_df_padj %>% 
  subset(padj < 0.2) %>%
  group_by(gene_id) %>%
  summarize(sample_types=list(unique(sample_type))) %>%
  ggplot(aes(sample_types)) +
  geom_bar() +
  scale_x_upset() +
  theme_classic()

plot(plt)
  
if (save_pdf) {
pdf(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/figure_axon_gene_group_upset.pdf"),
      height=5/2.54, width=7.5/2.54)
  plt
dev.off()
}


## ---------------------------------------------------------------------------------------------------------------
expressed_genes <- unique(valid_gene_df$GF)
expressed_genes <- sapply(expressed_genes, function(gene_id) strsplit(gsub("__ambiguous|\\[|\\]", "", gene_id), "\\+")[[1]])
expressed_genes <- unique(unlist(expressed_genes))

all_glm_res_df_padj_geneid_fix <- all_glm_res_df_padj %>%
  rowwise() %>%
  mutate(gene_id_fix=paste0(fix_gene_id(gene_id), collapse=";")) %>%
  separate_rows(gene_id_fix, sep=";") 
  
sig_gene_lst_df <- all_glm_res_df_padj_geneid_fix %>% 
  subset(padj < 0.2) %>%
  group_by(sample_type) %>% 
  summarize(gene_ids=list(unique(intersect(expressed_genes, gene_id_fix))))
sig_gene_lst <- sig_gene_lst_df$gene_ids
names(sig_gene_lst) <- sig_gene_lst_df$sample_type
sapply(sig_gene_lst, length)

# go enrich
enriched_goterms_sig <- compareCluster(sig_gene_lst,
                                       ont="CC", keyType="ENSEMBL", 
                                       OrgDb=org.Mm.eg.db, 
                                       universe=expressed_genes)

plt <- dotplot(enriched_goterms_sig) + 
  scale_color_gradient(low=accent_colors[1], high=background_colors[4]) + 
  theme(text=element_text(size=5))

plot(plt)

if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, 
                 "plots/illustrator_pdfs/figure_axon_go_enrichment_genes_nb_binom_model.pdf"),
       plot=plt,
                 height=150, width=150, units="mm",
       useDingbats=FALSE)
}


# plot only the top 10 terms with most abundant counts for proimal 720
get_gene_prop <- function(gene_ratio) {
  vals <- strsplit(gene_ratio, "\\/")[[1]]
  return(as.numeric(vals[1])/as.numeric(vals[2]))
}
res <- data.frame(enriched_goterms_sig@compareClusterResult) %>% 
  rowwise() %>%
  mutate(gene_ratio=get_gene_prop(GeneRatio))

res_prox <- (subset(res, Cluster == "Proximal_Axon_720") %>%
  arrange(desc(Count)))[1:10,]


res_dist_18 <- (subset(res, Cluster == "Distal_Axon_1080") %>%
  arrange(desc(Count)))[1:5,]

res_dist_12 <- (subset(res, Cluster == "Distal_Axon_720") %>%
                  arrange(desc(Count)))[1:5,]
plt <- rbind(res_prox,
      res_dist_12,
      res_dist_18) %>%
  drop_na() %>%
  ggplot(aes(Cluster, Description, size=gene_ratio, color=p.adjust)) +
  geom_point() +
  theme_bw() +
  scale_size_continuous(breaks=c(0.05, 0.1)) +
  scale_color_gradient(low=accent_colors[1], high=background_colors[4], limits=c(NA, 0.05)) 
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, 
                 "plots/illustrator_pdfs/figure_axon_go_enrichment_genes_nb_binom_model_FOR_DISSERTATION_PRESENTATION.pdf"),
       plot=plt,
       height=150, width=150, units="mm",
       useDingbats=FALSE)
}


## ---------------------------------------------------------------------------------------------------------------
axon_motif_savefile <- file.path(file_path_ref$project_directory, "data/SLAM_Axon_Analysis/rna_features/202308_RNA_motifs_axon.RData")

if (!file.exists(axon_motif_savefile) | overwrite) {
  # now load soma data and look for motif enrichment
  load(file.path(file_path_ref$project_directory, "data/sequencing/quantseq_3end_gene_annotated.RData"))
  
  grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                    dataset="mmusculus_gene_ensembl")
  
  
  soma_genes <- unique(getBM(c("ensembl_gene_id", "ensembl_transcript_id"),
                             filters="ensembl_transcript_id", 
                             values=quantseq_utr_annotated_filt$utr_name,
                             mart=grcm38_101)$ensembl_gene_id)
  
  utrs <- getBM(c("3utr", "ensembl_transcript_id", "ensembl_gene_id"),
                filters="ensembl_transcript_id",
                values=unique(quantseq_utr_annotated_filt$utr_name),
                mart=grcm38_101)
  
  utrs_filt <-  utrs %>% 
    subset(!grepl("S|N", `3utr`) & nchar(`3utr`) > 10  & nchar(`3utr`) < 10000) %>%
    inner_join(unique(quantseq_utr_annotated_filt[,c("utr_name", "bp_to_remove_from_3end", "peak_name")]),
               by=c("ensembl_transcript_id" = "utr_name")) %>%
    rowwise() %>%
    mutate(utr3_trunc=substr(`3utr`, 1, nchar(`3utr`)-bp_to_remove_from_3end),
           utr_name=sprintf("%s_%s_%s", ensembl_gene_id, ensembl_transcript_id, peak_name))
  
  save(utrs_filt, file=file.path(file_path_ref$project_directory, "data/SLAM_Axon_Analysis/rna_features/expressed_utr_sequences_truncated.RData"))
  
  get_utr_vec <- function(ens_gene_vec) {
    utrs_filt_subs <- subset(utrs_filt, ensembl_gene_id %in% ens_gene_vec)
    utrs_vec <- utrs_filt_subs$utr3_trunc
    names(utrs_vec) <- utrs_filt_subs$utr_name
    return(utrs_vec)
  }
  
  sig_utr_lst <- lapply(sig_gene_lst, get_utr_vec)
  
  # write these utrs to file
  for (utr_lst_name in names(sig_utr_lst)) {
    utr_names <- names(sig_utr_lst[[utr_lst_name]])
    utrs <- as.list(sig_utr_lst[[utr_lst_name]])
    
    write.fasta(utrs, utr_names, 
                file.path(file_path_ref$project_directory,
                          sprintf("data/SLAM_Axon_Analysis/rna_features/202308_%s.utr_sequences_cpn_valid_labeled_genes.fa",
                                  utr_lst_name)), 
                open = "w",
                nbchar = 10000, 
                as.string = TRUE)
  }
  
  bkg_utr_lst <- utrs_filt$utr3_trunc
  names(bkg_utr_lst) <- utrs_filt$utr_name
  
  write.fasta(as.list(bkg_utr_lst), names(bkg_utr_lst), 
              file.path(file_path_ref$project_directory,
                        sprintf("data/SLAM_Axon_Analysis/rna_features/202308_%s.utr_sequences_cpn_valid_labeled_genes.fa",
                                "Axon_Background")), 
              open = "w",
              nbchar = 10000, 
              as.string = TRUE)
  
  motif_enrich_axon <- transite::run_matrix_tsma(sig_utr_lst, bkg_utr_lst, n_cores = 5)
  
  save(motif_enrich_axon, file=axon_motif_savefile)
} else {
  load(axon_motif_savefile)
}

# 20230831 plot a dot-plot version for the dissertation presentation
valid_motifs <- union(subset(motif_enrich_axon$enrichment_dfs$Proximal_Axon_720, 
                       adj_p_value < 0.1 & enrichment > 1)$motif_id,
                      subset(motif_enrich_axon$enrichment_dfs$Distal_Axon_720, 
                             adj_p_value < 0.1 & enrichment > 1)$motif_id)


dotplot_df <- rbind(motif_enrich_axon$enrichment_dfs$Proximal_Axon_720 %>% 
                      subset(motif_id %in% valid_motifs) %>% 
                      rowwise() %>% 
                      mutate(cluster="Proximal_Axon_720"), 
                    motif_enrich_axon$enrichment_dfs$Distal_Axon_720 %>% 
                      subset(motif_id %in% valid_motifs) %>% 
                      rowwise() %>% mutate(cluster="Distal_Axon_720")) %>% 
  rowwise() %>% mutate(motif_name=paste0(c(motif_rbps, motif_id), collapse=";")) %>%
  arrange(enrichment)
dotplot_df$cluster <- factor(dotplot_df$cluster, levels=c("Proximal_Axon_720", "Distal_Axon_720"))
#dotplot_df$motif_name <- factor(dotplot_df$motif_name, levels=dotplot_df$motif_name)

plt <- dotplot_df %>% 
  ggplot(aes(cluster, motif_name, size=enrichment, color=adj_p_value)) + 
  geom_point() +
  theme_bw() +
  scale_color_gradient(low=accent_colors[1], high=background_colors[4])

plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, 
                 "plots/illustrator_pdfs/figure_axon_enriched_motifs_12h_only_FOR_DISSERTATION_PRESENTATION.pdf"),
       plot=plt,
       height=150, width=150, units="mm",
       useDingbats=FALSE)
}

enrichment_res <- motif_enrich_axon$enrichment_dfs

enriched_motifs <- unlist(sapply(enrichment_res,
                                 function(df) subset(df, adj_p_value < 0.1 & enrichment > 1)$motif_id))

to_plot <- data.frame(do.call(rbind,
                              lapply(names(enrichment_res),
                                     function(i) enrichment_res[[i]] %>%
                                       subset(motif_id %in% enriched_motifs) %>%
                                       rowwise() %>% mutate(stype=i))) %>%
                        pivot_wider(id_cols=c(motif_rbps, motif_id), names_from=stype, values_from=enrichment))

enriched_motifs_prox_dist_720 <- unlist(sapply(c("Distal_Axon_720", "Proximal_Axon_720"),
                                               function(name) subset(enrichment_res[[name]], adj_p_value < 0.1 & enrichment > 1)$motif_id))

to_plot_prox_dist_720 <- to_plot %>% subset(motif_id %in% enriched_motifs_prox_dist_720)

plt <- ggplot(to_plot_prox_dist_720, 
       aes(Proximal_Axon_720, Distal_Axon_720, label=motif_id)) + 
  geom_point() +
  geom_text_repel(data=to_plot_prox_dist_720, mapping=aes(Proximal_Axon_720, Distal_Axon_720, label=motif_id), max.overlaps = 20) +
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.1) + 
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.1) + 
  xlim(0.5, 1.5) +
  ylim(0.5, 1.5) +
  theme_classic() +
  theme(text=element_text(size=5))
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/figure_axon_motif_distal_proximal_720.pdf"),
       plot=plt,
       height=80,
       width=80, units="mm", useDingbats=FALSE)
}


enriched_motifs_prox720_dist1080 <- unlist(sapply(c("Proximal_Axon_720", "Distal_Axon_1080"),
                                               function(name) subset(enrichment_res[[name]], adj_p_value < 0.1 )$motif_id))

to_plot_dist <- to_plot %>% subset(motif_id %in% enriched_motifs_prox720_dist1080)

plt <- ggplot(to_plot_dist, aes(Proximal_Axon_720, Distal_Axon_1080, label=motif_id)) + 
  geom_point() +
  geom_text_repel(data=to_plot_dist, mapping=aes(Proximal_Axon_720, Distal_Axon_1080, label=motif_id), max.overlaps = 20) +
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.1) + 
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.1) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=5))

plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/figure_axon_motif_proximal_720_vs_distal_1080.pdf"),
       plot=plt,
       height=80,
       width=80, units="mm", useDingbats=FALSE)
} 


plt <- ggplot(to_plot, aes(Distal_Axon_720, Distal_Axon_1080, label=motif_id)) + 
  geom_point() +
  geom_text_repel(data=to_plot, mapping=aes(Distal_Axon_720, Distal_Axon_1080, label=motif_id), max.overlaps = 20) +
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.1) + 
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.1) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  theme(text=element_text(size=5))
plot(plt)
if (save_pdf) {
  ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/figure_axon_motif_distal_720_vs_distal_1080.pdf"),
       height=80,
       width=80, units="mm", useDingbats=FALSE)
}



## ---------------------------------------------------------------------------------------------------------------
# load in the AU(n) information from FIMO
# same as from Soma RNA features script, probably should make a proper helper script
motif_locs <- read.table(file.path(file_path_ref$project_directory, "data/SLAM_Axon_Analysis/rna_features/20230818_FIMO_all_axon.utr_sequences_cpn_valid_labeled_genes.for_FIMO.txt"), sep="\t", header = TRUE) 

# load UTR sequences
# called utrs_filt
load(file.path(file_path_ref$project_directory, "data/SLAM_Axon_Analysis/rna_features/expressed_utr_sequences_truncated.RData"))

collapse_motif_locs <- function(peak_name) {
  
  peak_motifs <- subset(motif_locs, sequence_name == peak_name)
  peak_motif_range <- data.frame(reduce(IRanges(start=peak_motifs$start, end=peak_motifs$stop)))
  if (length(peak_motif_range) > 0) {
    return(data.frame(max_width=max(peak_motif_range$width),
                      num_motifs=nrow(peak_motif_range),
                      total_coverage=sum(peak_motif_range$width),
                      starts=paste0(peak_motif_range$start, collapse=";"),
                      ends=paste0(peak_motif_range$end, collapse=";"),
                      max_start=peak_motif_range$start[which.max(peak_motif_range$width)],
                      max_end=peak_motif_range$end[which.max(peak_motif_range$width)],
                      peak_name=peak_name))
  } else {
    return(data.frame(max_width=0,
                      num_motifs=0,
                      total_coverage=0,
                      max_start=NA,
                      max_end=NA,
                      starts="",
                      ends="",
                      peak_name=peak_name))
  }
  
}

motif_locs_info <- data.frame(do.call(rbind,
                                      lapply(unique(motif_locs$sequence_name), 
                                             collapse_motif_locs)))

#utrs_no_motif <- setdiff(quantseq_utr_annotated_filt_utr$peak_name, motif_locs_info$peak_name)
utrs_no_motif <- setdiff(utrs_filt$peak_name, motif_locs_info$peak_name)

motif_locs_info_annot <- rbind(motif_locs_info,
  data.frame(max_width=0, num_motifs=0, total_coverage=0, 
             starts=NA, ends=NA, max_start=NA, max_end=NA, 
             peak_name=utrs_no_motif)) %>%
  #inner_join(quantseq_utr_annotated_filt_utr[,c("peak_name", "utr3_trunc")]) %>%
  inner_join(utrs_filt[,c("peak_name", "utr3_trunc")]) %>%
  rowwise() %>%
  mutate(#is_dist12=peak_name %in% motif_genes_rbms$fg_motif$peak_name,
         utr3_length=nchar(utr3_trunc)) 

# also compare with all axon peaks
# remove 18h timepoint because only <=10 genes

peak_axon_sample_map <- all_glm_res_df_padj_geneid_fix %>%
  rowwise() %>%
  mutate(axon_group=ifelse(padj < 0.2 & !grepl("1080", sample_type), 
                           sample_type, "Soma")) %>%
  inner_join(utrs_filt, by=c("gene_id_fix" = "ensembl_gene_id")) %>%
  dplyr::select(c(peak_name, gene_id_fix, axon_group)) %>%
  unique()

### START FIGURE: barplot comparing proportion with length au repeat
au_repeat_bins <- c("0-6nt", "7-8nt", "9-10nt", ">10nt")
motif_locs_info_all_axon <- motif_locs_info_annot %>%
  inner_join(peak_axon_sample_map) %>%
  rowwise() %>%
  mutate(
         au_width_group=factor(au_repeat_bins[.bincode(max_width,
                                                       breaks=c(-1, 6, 8, 10, Inf), 
                                                       include.lowest = TRUE)],
                               levels=au_repeat_bins),
         au_width_norm=max_width/utr3_length*1000)

motif_locs_info_all_axon_maxwidth_bin <-
  motif_locs_info_all_axon %>%
  group_by(axon_group) %>%
  summarize(total_motifs=n()) %>%
  inner_join(motif_locs_info_all_axon, by="axon_group") %>%
  rowwise() %>%
  group_by(axon_group, au_width_group) %>%
  summarize(prop_max_width=n()/total_motifs[1]) %>%
  subset(!grepl("1080", axon_group))

plt <- motif_locs_info_all_axon_maxwidth_bin %>%
  ggplot(aes(au_width_group, prop_max_width, fill=axon_group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9)) +
  scale_fill_manual(values=c(accent_colors[2], accent_colors[4], background_colors[5])) +
  theme_classic() +
  theme(text=element_text(size=7))
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/20230818_axon_utrs_length_longest_au_motif_barplot.pdf"),
       plot=plt,
       height=50, width=80, units="mm")
}

## FIGURE: now look at relative location of motif in UTR
motif_locs_nomerge <- peak_axon_sample_map %>%
  right_join(motif_locs, by=c("peak_name" = "sequence_name" )) %>%
  inner_join(utrs_filt[,c("peak_name", "utr3_trunc")],
             by="peak_name") %>% 
  rowwise() %>%
  mutate(axon_group=ifelse(is.na(axon_group), "Soma", axon_group)) %>%
  subset(!grepl("1080", axon_group))

plt <- motif_locs_nomerge %>% ggplot(aes(start/nchar(utr3_trunc), color=axon_group)) +
  geom_density() +
  scale_color_manual(values=c(accent_colors[2], accent_colors[4], background_colors[5])) +
  theme_classic() +
  theme(text=element_text(size=7)) +
  xlab("Relative Motif Location in 3'UTR")
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/20230818_axon_utrs_location_of_AU_motifs_relative.pdf"),
       plot=plt,
       height=50, width=80, units="mm")
}

plt <- motif_locs_nomerge %>% ggplot(aes(start, color=axon_group)) +
  geom_density() +
  scale_color_manual(values=c(accent_colors[2], accent_colors[4], background_colors[5])) +
  theme_classic() +
  theme(text=element_text(size=7)) +
  xlab("Absolute Motif Location in 3'UTR") +
  scale_x_log10()
plot(plt)
if (save_pdf) {
ggsave(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/20230818_axon_utrs_location_of_AU_motifs_absolute.pdf"),
       plot=plt,
       height=50, width=80, units="mm")
}


