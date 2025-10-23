library(MASS)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
loc <- args[1]
stime <- args[2]
model_type <- args[3]

project_dir <- "/n/holyscratch01/macklis_lab/priyav/lab_notebook/SLAM_final/20230818_axon_conversion_counts"

conversion_df <- read.csv(file.path(project_dir, "20230625_axon_agg_counts.csv")) 

metadata <- read.csv(file.path(project_dir, "SLAM_final_design.csv"))

intron_median_rates <- read.csv(file.path(project_dir, 
                                          "20230625_axon_intron_median_rates.csv"))

valid_gene_df <- read.csv(file.path(project_dir, 
                           "20230818_axon_valid_genes.csv"))

meta_info_group_number <- metadata %>%
  subset(!grepl("R5_UPRT_12h_Dist|R5_Neg_12h_Dist", sample_id) & location %in% c("Distal_Axon", "Proximal_Axon")) %>% 
  group_by(location, sample_time_min) %>%
  summarize(min_num_sample_group=min(sum(IUE_plasmid== "pPV25"), sum(IUE_plasmid == "dollar_51")))

reads_per_sample <- conversion_data %>%
  group_by(sample_id) %>%
  summarize(total_reads=sum(reads),
            total_TC_reads=sum(TC_reads))

conv_info <- conversion_data %>% 
  subset(!grepl("R5_UPRT_12h_Dist|R5_Neg_12h_Dist", sample_id)) %>% 
  inner_join(meta_info_group_number, by=c("location", "sample_time_min")) %>%
  inner_join(intron_median_rates_annot[, c("sample_id", "intron_multiplier")]) %>% 
  inner_join(reads_per_sample, by="sample_id") %>%
  rowwise() %>% 
  mutate(TC_per_1k_T=nTC/nT*1000*intron_multiplier, 
         TC_reads_per_mil=TC_reads/total_reads*1e6,
         TC_reads_per_milTC=TC_reads/total_TC_reads*1e6,
         reads_per_mil=reads/total_reads*1e6,
         read_valid=reads > 50,
         sample_type=sprintf("%s_%s", location, sample_time_min))
                          
#conversion_data <- conversion_df %>% 
#  inner_join(metadata, by="sample_id") %>% 
#  inner_join(intron_median_rates, by="sample_id") 
run_fit_glm_nb <- function(gene_subs, fit_name, design, return_model=FALSE) {
  tryCatch({
    #fit3 <- glm.nb(tc_intron_scaled~IUE_plasmid, data=gene_subs)
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

run_t_read_scaled_glm <- function(geneid, loc, stime, design, input_data=conv_info, run_type="all") {
  
  gene_subs <- subset(input_data, GF == geneid & location == loc & sample_time_min == stime) %>%
    rowwise() %>% mutate(tc_scaled=round(TC_reads_per_milTC/reads_per_mil*250)) # based on the median of 250 TC reads ish across all samples
  print(gene_subs)
  fit_lst <- list()
  
  if (run_type == "all") {
    fit_lst[[1]] <- run_fit_glm_nb(gene_subs, "all", design)
  } else if (run_type == "cv") {
    
    
    # now leave one out xval
    for (i in 1:nrow(gene_subs)) {
      keep_idx <- setdiff(1:nrow(gene_subs), i)
      out_sample <- gene_subs$sample_id[i]
      glm_res <- run_fit_glm_nb(gene_subs[keep_idx, ], sprintf("%s_excluded", out_sample), design, return_model=TRUE)
      glm_res_df <- glm_res$res_df
      
      
      if (!is.null(glm_res$model)) {
        prediction <- predict(glm_res$model, newdata=data.frame(IUE_plasmid=gene_subs$IUE_plasmid[i]))
        truth <- log(gene_subs$tc_scaled[i])
        glm_res_df$prediction <- prediction
        glm_res_df$truth <- truth
      } else {
        glm_res_df$prediction <- NA
        glm_res_df$truth <- NA
      }
      
      fit_lst[[(i)]] <- glm_res_df
      
      
    }
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

run_t_intron_scaled_glm <- function(geneid, loc, stime, design, input_data=conversion_data, intron_scale_column="median_intron_convrate", run_type="all") {
  
  gene_subs <- subset(input_data, GF == geneid & location == loc & sample_time_min == stime)
  t_conv_scale <- median(gene_subs$nT)/gene_subs$nT
  intron_scale <- median(unlist(gene_subs[intron_scale_column]))/unlist(gene_subs[intron_scale_column])
  
  gene_subs$tc_t_scaled=round(gene_subs$nTC*t_conv_scale)
  gene_subs$tc_scaled=round(gene_subs$tc_t_scaled*intron_scale)
  
  fit_lst <- list()
  
  if (run_type == "all") {
    fit_lst[[1]] <- run_fit_glm_nb(gene_subs, "all", design)
  } else if (run_type == "cv") {
    
    
    # now leave one out xval
    for (i in 1:nrow(gene_subs)) {
      keep_idx <- setdiff(1:nrow(gene_subs), i)
      out_sample <- gene_subs$sample_id[i]
      glm_res <- run_fit_glm_nb(gene_subs[keep_idx, ], sprintf("%s_excluded", out_sample), design, return_model=TRUE)
      glm_res_df <- glm_res$res_df
      
      
      if (!is.null(glm_res$model)) {
        prediction <- predict(glm_res$model, newdata=data.frame(IUE_plasmid=gene_subs$IUE_plasmid[i]))
        truth <- log(gene_subs$tc_scaled[i])
        glm_res_df$prediction <- prediction
        glm_res_df$truth <- truth
      } else {
        glm_res_df$prediction <- NA
        glm_res_df$truth <- NA
      }
      
      fit_lst[[(i)]] <- glm_res_df
      
      
    }
  } else if (run_type == "bootstrap") {
    
    set.seed(67) # important to set the seed so that all genes have same random ordering
    lst_idx <- 1 #nrow(gene_subs) + 1
    for (j in 1:100) {
      random_idx <- sample(1:nrow(gene_subs), size=nrow(gene_subs), replace=TRUE) # bootstrap
      #subs_test <- gene_subs
      subs_test$tc_scaled <- gene_subs$tc_scaled[random_idx]
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

# requres a df with sample time, and loc facets
# and also pvalue column on which to do adjustment
adjust_p_values_condition <- function(df){
  res_lst <- list()
  for (l in unique(df$location)) {
    for (s in unique(df$sample_time_min)) {
      stype <- sprintf("%s_%s", l, s)
      stype_res <- subset(df, sample_time_min == s & location == l)
      stype_res$padj <- p.adjust(stype_res$pvalue, method="BH")
      stype_res$sample_type <- stype
      res_lst[[stype]] <- stype_res
    }
  }
  
  return(data.frame(do.call(rbind, res_lst)))
}

valid_gene_subs <- subset(valid_gene_df, location == loc & sample_time_min == stime)

if (model_type == "full") {
all_glm_res <- list()
for (i in 1:nrow(valid_gene_subs)) {
  if (i %% 50 == 0) {
    print(i)
  }
  g <- valid_gene_subs$GF[i]
  l <- valid_gene_subs$location[i]
  s <- valid_gene_subs$sample_time_min[i]
  all_glm_res[[i]] <- run_t_read_scaled_glm(g, l, s, "tc_scaled~IUE_plasmid")
  
}
all_glm_res_df <- data.frame(do.call(rbind, all_glm_res))
all_glm_res_df$pvalue <- 1-pnorm(all_glm_res_df$z_value_IUE_plasmidpPV25)
all_glm_res_df_padj <- adjust_p_values_condition(all_glm_res_df)

print(sprintf("Num Sig padj < 0.2: %s", nrow(subset(all_glm_res_df_padj, padj < 0.2))))
print(sprintf("Num Sig padj < 0.1: %s", nrow(subset(all_glm_res_df_padj, padj < 0.1))))
print(sprintf("Num Sig padj < 0.05: %s", nrow(subset(all_glm_res_df_padj, padj < 0.05))))

### reduced model

all_glm_res_reduced <- list()
for (i in 1:nrow(valid_gene_subs)) {
  if (i %% 50 == 0) {
    print(i)
  }
  g <- valid_gene_subs$GF[i]
  l <- valid_gene_subs$location[i]
  s <- valid_gene_subs$sample_time_min[i]
  all_glm_res_reduced[[i]] <- run_t_read_scaled_glm(g, l, s, "tc_scaled~1")
  
}
all_glm_res_reduced_df <- data.frame(do.call(rbind, all_glm_res_reduced))


all_glm_compare <- all_glm_res_df_padj %>% 
  inner_join(all_glm_res_reduced_df, by=c("sample_time_min", "location", "gene_id"),
             suffix=c(".full", ".reduced")) %>%
  rowwise() %>%
  mutate(lr_stat=-twologlik.reduced+twologlik.full) %>%
  rowwise() %>%
  mutate(lr_pvalue=pchisq(lr_stat, 1, lower.tail=FALSE))

all_glm_compare_lr_padj <- adjust_p_values_condition(
  all_glm_compare %>% dplyr::select(c(location, sample_time_min, gene_id, lr_pvalue)) %>%
    dplyr::rename(pvalue=lr_pvalue))

save(all_glm_res_df, all_glm_res_df_padj, 
     all_glm_res_reduced_df, 
     all_glm_compare_lr_padj, 
     file=file.path(project_dir, sprintf("20230818_axon_neg_binomial_GLM_results.%s_%s.%s.RData", loc, stime, model_type))) #_randomize_data.RData"))

}

if (model_type == "cv") {
### cross val model
cv_glm_res <- list()
for (i in 1:nrow(valid_gene_subs)) {
  if (i %% 50 == 0) {
    print(i)
  }
  g <- valid_gene_subs$GF[i]
  l <- valid_gene_subs$location[i]
  s <- valid_gene_subs$sample_time_min[i]
  cv_glm_res[[i]] <- run_t_read_scaled_glm(g, l, s, "tc_scaled~IUE_plasmid", run_type="cv")
}
cv_glm_res_df <- data.frame(do.call(rbind, cv_glm_res))
cv_glm_res_df <- cv_glm_res_df %>%
  rowwise() %>% 
  mutate(sample_id=gsub("_excluded", "", glm_nb_fit_name)) %>%
  inner_join(metadata, by=c("sample_id", "sample_time_min", "location"))

agg_cv_res <- function(cv_res_df) {
  agg_res <- cv_res_df %>% 
    group_by(sample_time_min, location, gene_id) %>%
    summarize(ss_resid=sum((prediction-truth)^2),
              ss_tot=sum((truth-mean(truth))^2),
              err_uprt=mean((prediction-truth)[IUE_plasmid=="pPV25"]),
              err_neg=mean((prediction-truth)[IUE_plasmid=="dollar_51"]),
              mean_uprt=mean(truth[IUE_plasmid=="pPV25"]),
              mean_neg=mean(truth[IUE_plasmid=="dollar_51"]),
              mean_express=mean(truth)) %>%
    rowwise() %>%
    mutate(valid_rsquared=1-ss_resid/ss_tot)
  
  return(agg_res)
  
}


cv_glm_res_df_rsq <- agg_cv_res(cv_glm_res_df)


save(
     cv_glm_res_df_rsq, cv_glm_res_df, 
     file=file.path(project_dir, sprintf("20230818_axon_neg_binomial_GLM_results.%s_%s.%s.RData", loc, stime, model_type))) #_randomize_data.RData"))
}

if (model_type == "boot") {
  ### bootstrap model
  boot_glm_res <- list()
  for (i in 1:nrow(valid_gene_subs)) {
    if (i %% 50 == 0) {
      print(i)
    }
    g <- valid_gene_subs$GF[i]
    l <- valid_gene_subs$location[i]
    s <- valid_gene_subs$sample_time_min[i]
    boot_glm_res[[i]] <- run_t_read_scaled_glm(g, l, s, "tc_scaled~IUE_plasmid", run_type="bootstrap")
  }
  boot_glm_res_df <- data.frame(do.call(rbind, boot_glm_res))
 
  
  save(
    boot_glm_res_df, 
    file=file.path(project_dir, sprintf("20230818_axon_neg_binomial_GLM_results.%s_%s.%s.RData", loc, stime, model_type))) #_randomize_data.RData"))
}