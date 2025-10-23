## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------------------
library(MASS)
library(dplyr)
library(data.table)
library(ggrastr)
library(purrr)
library(tidyr)
library(stats)
library(preprocessCore)

### NEED TO SET TO PROPER DIRECTORY
source("Y:/users/priyav/slam_final/scripts/helper_script.R")

plot_subdir <- file.path(file_path_ref$project_directory, "plots/SLAM_Soma_NTR_and_Kinetics")
data_subdir <- file.path(file_path_ref$project_directory, "data/SLAM_Soma_NTR_and_Kinetics")


## ---------------------------------------------------------------------------------------------------------------
# load the counts data
load(file.path(file_path_ref$project_directory, 
               file_path_ref$results$soma_counts_data))

# retain only genes with at least 50 counts in 2/3 of all samples
# and with at least 3 samples (approx one timepoint) that have greater than 10 T->C converted reads
valid_cpn_genes <- unique((counts_soma %>%
  group_by(Geneid) %>%
  summarize(valid_unlabeled=sum(reads_total > 50)/length(reads_total) > 2/3,
            valid_labeled=sum(tc_geq1_reads_total > 10) >= 3) %>%
  subset(valid_unlabeled & valid_labeled))$Geneid)

# get conversion rates per gene
cpn_conversion_rates_all <- counts_soma %>%
  rowwise() %>%
  mutate(conversion_rate=TC_count_total/T_count_total,
         sample_time_hours=sample_time_min/60)

# get conversion rates per whole sample
cpn_global_conversion_rate <- cpn_conversion_rates_all %>%
  group_by(sample_id, sample_time_hours, IUE_plasmid) %>%
  summarize(conversion_rate=sum(TC_count_total)/sum(T_count_total)) %>%
  subset(sample_time_hours %in% c(0, 2, 4, 6, 8)) %>%
  rowwise() %>%
  mutate(sample_group=paste0(sample_time_hours, IUE_plasmid))

# plot conversion rates per whole sample for Supplement 1
ggplot(cpn_global_conversion_rate, aes(sample_time_hours, conversion_rate, 
                                       group=sample_group,
                                       fill=IUE_plasmid)) + 
  geom_boxplot(linewidth=0.1) + geom_point(position=position_dodge(1.6), size=0.2) + 
  theme_classic() + 
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8)) +
  scale_fill_manual(values=c(background_colors[5], accent_colors[2])) +
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "Supplement1_conversion_rate_SLAM_soma_whole_sample.pdf"),
       height=30, width=50, units="mm")

# remove the 0 hour timepoint (no label)
cpn_conversion_rates <- subset(cpn_conversion_rates_all, 
                               IUE_plasmid == "pPV25" & sample_time_min %in% c(120, 240, 360, 480)) 



## ---------------------------------------------------------------------------------------------------------------
# define the soma metadata for New Txpt Prop --> use only UPRT+ CPN soma samples
# also only use times 0-8h
meta_soma <- subset(metadata, 
                      location == "Soma" & IUE_plasmid == "pPV25" & sample_time_min %in% c(0,120, 240, 360, 480) )
  
cB_files <- file.path(file_path_ref$project_directory, 
                            file_path_ref$results$conversion_counts_dir,
                            sprintf("%s.conversions.bakr.csv", meta_soma$sample_id))
  
cB <- data.frame(do.call(rbind, 
                         lapply(cB_files, function(x) read.csv(x) %>% mutate(sample=gsub(".conversions.bakr.csv", "", basename(x)))))) 
cB$sample_id <- NULL


# find number of reads in each sample
# with a given number TC conversions (TC) and number of T (nT) 
cB_grouped <- data.frame(cB %>%
      group_by(TC, nT, sample) %>%
      summarize(n=sum(n)))

# give a dummy gene name 
cB_grouped$XF <- "all_genes"

# estimate p_old using 0h (no label samples)
# and just use the average conversion rate in the old samples
# averaging both within a sample and between all the samples
polds <- cB_grouped %>% inner_join(meta_soma, by=c("sample" = "sample_id")) %>%
  subset(sample_time_min == 0) %>%
  group_by(sample) %>%
  summarize(convrate=sum(TC*n)/sum(nT*n))
# between the samples 
pold_mean <- mean(polds$convrate)

# Vector-compatible function to obtain binomial mixture likelihood for observing a 
# number of conversions per read given the number of t's in the read
# Arguments
# nT: integer vector number of T's in a read
# pold: float converison rate in unlabeled RNA (errors)
# TC: integer vector number of TC conversions in a read, same length as TC
# n : number of reads with the given nT and TC values
# params: float vector of length 2: c(theta, pnew) 
#         where theta is the proportion of new transcripts
#         and where pnew is the probability of observing a T->C conversion in new RNA 
# Returns
#   float log likelihood
mixed_lik_binom <- function(nT, pold, TC, n, params){
    pnew <- params[2]
    theta <- params[1]
 
    logl <- sum(n*log(theta*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(pnew^TC)*((1 -pnew)^(nT-TC)) +  (1-theta)*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(pold^TC)*((1 - pold)^(nT-TC)) ) )

          return(logl)

        }

# Vector-compatible function to obtain poisson mixture likelihood for observing a 
# number of conversions per read given the number of t's in the read
# Arguments
# nT: integer vector number of T's in a read
# pold: float converison rate in unlabeled RNA (errors)
# TC: integer vector number of TC conversions in a read, same length as TC
# n : number of reads with the given nT and TC values
# params: float vector of length 2: c(theta, pnew) 
#         where theta is the proportion of new transcripts
#         and where pnew is the probability of observing a T->C conversion in new RNA 
# Returns
#   float log likelihood
mixed_lik_pois <- function(nT, pold, TC, n, params){
  pnew <- params[2]
  theta <- params[1]
  lam_n <- nT*pnew
  lam_o <- nT*pold
  logl <- sum(n*log(theta*lam_n^TC*exp(-lam_n)/factorial(TC) + (1-theta)*(lam_o^TC)*exp(-lam_o)/factorial(TC)))
  return(logl)
}

# Function to jointly find the maximum likelihood estimates for the per-sample 
# proportion of new transcript (theta) and the probability of observing a conversion 
# in new RNA (pnew). Uses both binomial mixture model and poisson mixture model
# Arguments:
#   sample_id: sample_id, must exist in cB_grouped
#   pold: float likelihood of conversion in old RNA (errors); default is pold_mean
#   conv_df: data frame with columns nT (int), TC (int), n (int)
# Returns:
#   data frame of size 1x5, with columns containing the maximum likelihood estimates
#   colnames: sample_id, theta_binomial, pnew_binomial, theta_poisson, pnew_poisson
get_sample_stats <- function(sample_id, pold=pold_mean, conv_df=cB_grouped) {
  sample_data <- na.omit(conv_df %>% subset(sample == sample_id))
  print(sample_id)
  
  binom_res <- stats::optim(par=c(0.5, 0.01), mixed_lik_binom, 
                            nT=sample_data$nT, TC=sample_data$TC, 
                            n=sample_data$n,  pold=pold, lower=0.0001, upper=0.9999, 
                            method = "L-BFGS-B", control=list(fnscale=-1))$par
  
  pois_res <- stats::optim(par=c(0.5, 0.01), mixed_lik_pois, 
                           nT=sample_data$nT, TC=sample_data$TC, 
                           n=sample_data$n, pold=pold, lower=0.0001, upper=0.9999, 
                           method = "L-BFGS-B", control=list(fnscale=-1))$par
  
  return(data.frame(sample_id=sample_id,
                    theta_binomial=binom_res[1],
                    pnew_binomial=binom_res[2],
                    theta_poisson=pois_res[1],
                    pnew_poisson=pois_res[2]))
}


# calculate the MLEs for theta and pnew for all samples
mle_params_soma <- 
  data.frame(do.call(rbind, 
                     lapply(unique(cB_grouped$sample), get_sample_stats))) %>%
  inner_join(metadata, by="sample_id") %>%
  rowwise() %>%
  mutate(sample_time_hours = sample_time_min/60)%>%
  subset(sample_time_hours > 0)

# plot the decreasing pnew as a function of time 
mle_params_soma %>%
subset(sample_time_hours > 0) %>%
ggplot(aes(sample_time_hours, pnew_poisson)) +
geom_point() +
theme_classic() +
scale_x_continuous(breaks=c(2, 4, 6, 8)) 

ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_whole_sample_pnew_estimate.pdf"),
        height=50, width=50, units='mm')



## ---------------------------------------------------------------------------------------------------------------
# Estimate sample-wide kinetics using two approaches
# approach 1: 1-e^(-kt)
theta_global_no_scale <- nls(theta_poisson~1-exp(-l*(sample_time_hours)), data=mle_params_soma, start=c(l=0.1))

# approach 2 (scaled): 1-C*e^(-kt)
theta_global_scaled <- nls(theta_poisson~1-C*exp(-l*(sample_time_hours)), data=mle_params_soma, start=c(C=1, l=0.1))

# plot the two approaches for Supplementary Figure on estimating kinetic parameters
mle_params_soma %>%
  subset(sample_time_hours > 0) %>%
  ggplot(aes(sample_time_hours, theta_poisson)) +
  geom_point(size=0.1) +
  geom_line(data=data.frame(sample_time_hours=1:9,
                            theta_poisson=1-coef(theta_global_scaled)[1]*exp(-coef(theta_global_scaled)[2]*1:9)),
            color=accent_colors[2], linewidth=0.2) +
  geom_line(data=data.frame(sample_time_hours=1:9,
                            theta_poisson=1-exp(-coef(theta_global_no_scale)[1]*1:9)),
            color=accent_colors[8], linewidth=0.2) +
  theme_classic() +
  scale_x_continuous(breaks=c( 2, 4, 6, 8)) +
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_kinetic_fit_theta.pdf"), 
       height=50, width=50, units="mm")


## ---------------------------------------------------------------------------------------------------------------
# Vector-compatible function to obtain poisson mixture likelihood for observing a 
# number of conversions per read given the number of t's in the read
# Arguments
# nT: integer vector number of T's in a read
# pold: float probability of observing a conversion in old RNA (errors)
# pnew: float probability of observing a conversion in new RNA
# TC: integer vector number of TC conversions in a read, same length as TC
# n : number of reads with the given nT and TC values
# params: float vector of length 1: c(theta) 
#         where theta is the proportion of new transcripts
# Returns
#   float log likelihood
mixed_lik_gene_pois <- function(nT, pnew, pold, TC, n, params){
  theta <- params[1]
  lam_n <- nT*pnew
  lam_o <- nT*pold
  logl <- sum(n*log(theta*lam_n^TC*exp(-lam_n)/factorial(TC) + (1-theta)*(lam_o^TC)*exp(-lam_o)/factorial(TC)))
  return(logl)
}

# Function to find the maximum likelihood estimates for the per-sample, per gene
# proportion of new transcript (theta). Uses poisson mixture model
# Arguments:
#   sample_id: sample_id, must exist in cB
#   pold: float likelihood of conversion in old RNA (errors)
#   pnew: float probability of observing a conversion in new RNA
#   geneid : character, must exist in cB$XF
#   df_subs: data frame with columns XF (character), nT (int), TC (int), n (int)
# Returns:
#   data frame of size 1x3, with colnames: sample_id, gene_id, theta_poisson_MLE
get_gene_stats <- function(pnew, sample_id, geneid, pold, df_subs=cB) {
  sample_data <- na.omit(df_subs %>% subset(XF == geneid & sample == sample_id))
  

  
  pois_res <- stats::optim(par=0.5, mixed_lik_gene_pois, 
                           nT=sample_data$nT, TC=sample_data$TC, 
                           n=sample_data$n, pnew=pnew, pold=pold, 
                           lower=0.0001, upper=0.9999, 
                           method = "L-BFGS-B", control=list(fnscale=-1))$par
  
  return(data.frame(sample_id=sample_id,
                    geneid=geneid,
                    theta_poisson=pois_res[1]))
}

# Function to find the MLE estimates for new transcript proportion for a gene
# across all samples contained in cB
# Arguments:
#   geneid: sample_id, must exist in cB
#   pold: float likelihood of conversion in old RNA (errors); default is pold_mean
#   conv_df: data frame with columns XF (character), nT (int), TC (int), n (int); default is cB
#    soma_mle_est: data frame with columns sample_id (character), sample_time_hours (int), pnew_poisson (float) containing the pnew poisson MLE for each sample
# Returns:
#   data frame of size num_samplesx4, with colnames: sample_id (character), gene_id (character), theta_poisson_MLE (float), sample_time_hours (int)
get_gene_stats_all_sample <- function(geneid, p_old=pold_mean,
                                      conv_df=cB, soma_mle_est=mle_params_soma) {
  
  # get conversion info for gene across all samples
  subs_gene <- subset(conv_df, XF == geneid)
  
  # estimate theta per gene
  resdf <- data.frame(do.call(rbind,
                              lapply(1:nrow(soma_mle_est), 
                                     function(i) get_gene_stats(soma_mle_est$pnew_poisson[i], 
                                                                soma_mle_est$sample_id[i], geneid,
                                                                df_subs=subs_gene,
                                                                pold=p_old))
  )
  )
  # add time info
  resdf$sample_time_hours <- soma_mle_est$sample_time_hours
  return(resdf)
}

# Estimate kinetic parameters for a given gene using nonlinear least squares
# (nls) package in R 
# Arguments:
#   gene_stats_res: output of gene_gene_stats_all_sample(); data frame of size 
#      num_samplesx4, with colnames: sample_id (character), gene_id (character),
#      theta_poisson_MLE (float), sample_time_hours (int)
#   C: float value to scale the exponential model for estimating kdeg 
#      given theta and time: theta=1-C*e^(-k*t)
#   sample_time_hours_subset: integer vector of sample time hours to use to estimate
#      kinetics; defaults to c(2, 4, 6, 8)
# Returns:
#   data frame of size 1x(6+num_timepoints) with the following columns:
#   c(k_nls, C_nls, k_linear, rsquared_nls, rsquared_linear, geneid, <time>h for each <time>)
#   k_nls, C_nls correspond to k_deg and C that are optimal for the nls model
#   note that C_nls is the input C unless it is NA, in which case it is estimated
#   k_linear corresponds to the kdeg rate calculated from log(1-theta)~time+0
#   rsquared_nls and rsquared_linear are the rsquared goodness of fit coefficients
#   for each model
#   geneid is a character and self explanatory
#   <time>h refers to the average residual (model_expected-true_value) for each timepoint. This is used to find systematic trends in model performance
est_kdeg_gene_from_theta <- function(gene_stats_res, 
                                     C, sample_time_hours_subset=c(2,4,6,8)) { 
  
  # get specified time subset
  gene_stats_res <- gene_stats_res %>% 
    subset(sample_time_hours %in% 
           sample_time_hours_subset)

 tryCatch({
   if (is.na(C)) {
      poisson_res <- nls(theta_poisson~(1-C *exp(-k*sample_time_hours)), 
                       data=gene_stats_res, start=c(C=1, k=0.1))
    k <- coef(poisson_res)[2]
    C <- coef(poisson_res)[1]
   } else {
     
    poisson_res <- nls(theta_poisson~(1-C *exp(-k*sample_time_hours)), 
                       data=gene_stats_res, start=c( k=0.1))
    k <- coef(poisson_res)[1]
   }
   gene_stats_res <- gene_stats_res %>%
     inner_join()
    lm_res <- lm(log(1-theta_poisson)~sample_time_hours+0, data=gene_stats_res)
    k_linear <- -coef(lm_res)[1]
    
    stimes <- gene_stats_res$sample_time_hours
    expected_nls <- (1-C*exp(-k*stimes)) 
    expected_linear <- (1-exp(-k_linear*stimes))
    true <- gene_stats_res$theta_poisson
    
    ssr_nls <- sum((expected_nls-true)^2)
    ssr_linear <- sum((expected_linear-true)^2)
    sst <- sum((true-mean(true))^2)
    rsquared_nls <- 1-(ssr_nls/sst)
    rsquared_linear <- 1-(ssr_linear/sst)
    
        
    # get residuals by time
    resid_df <- data.frame(stime=sprintf("%sh",stimes),
                           resid=expected_nls-true) %>%
      group_by(stime) %>%
      summarize(mean_residual=mean(resid)) %>%
      pivot_wider(names_from=stime, values_from=mean_residual)

    
    return(cbind(data.frame(
                      k_nls=k,
                      C_nls=C,
                      k_linear=k_linear,
                      rsquared_nls=rsquared_nls,
                      rsquared_linear=rsquared_linear,
                      geneid=gene_stats_res$geneid[1]),
                 resid_df))
  }, error=function(e) {
    
    resid_df <- data.frame(stime=sprintf("%sh", unique(gene_stats_res$sample_time_hours)),
                           mean_residual=NA) %>%
      pivot_wider(names_from=stime, values_from=mean_residual)
    return(cbind(data.frame(
                      k_nls=NA,
                      C_nls=NA,
                      k_linear=NA,
                      rsquared_nls=NA,
                      rsquared_linear=NA,
                      geneid=gene_stats_res$geneid[1]),
                 resid_df))
  })
}



# get a list then data frame of thetas per gene per sample
# uses the function get_gene_stats_all_sample
gene_thetas <- lapply(valid_cpn_genes, get_gene_stats_all_sample)
gene_thetas_df <- data.frame(do.call(rbind, gene_thetas))

# get a list then data frame of gene kinetic parameters (k_deg) 
# first don't perform correction for decreasing labeling rate
gene_kinetics_df <- data.frame(
  do.call(rbind, lapply(gene_thetas, 
                        function(g) est_kdeg_gene_from_theta(g, 1))
  )
)

# next do perform correction for decreasing labeling rate (calculated from soma bulk timecourse)
gene_kinetics_theta_correction_df <- data.frame(
  do.call(rbind, lapply(gene_thetas, 
                        function(g) est_kdeg_gene_from_theta(g, 1.12796))))

# next do perform correction for decreasing labeling rate but estimate the correction factor per gene
gene_kinetics_theta_correction_genewise_df <- data.frame(
  do.call(rbind, lapply(gene_thetas, 
                        function(g) est_kdeg_gene_from_theta(g, NA))))



## ---------------------------------------------------------------------------------------------------------------
# plot gene thetas by sample
ggplot(gene_thetas_df, aes(factor(sample_time_hours), 
                           theta_poisson, group=sample_id, 
                           fill=factor(sample_time_hours))) + 
  geom_violin(draw_quantiles = 0.5, linewidth=0.1) +
  theme_classic() +
  scale_fill_manual(values=accent_colors[c(10, 8, 4,2)]) +
  theme(legend.position = "None")
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_theta_per_sample_violin.pdf"), 
       height=50, width=100, units="mm")



# anova residuals
theta_anova_resid <- data.frame(do.call(rbind,
                                        lapply(gene_thetas, function(df) {
  res_aov <- aov(theta_poisson ~ sample_time_hours,
                data = df)
  
  df$log_theta <- log(df$theta_poisson)
  res_aov_log <- aov(log_theta ~ sample_time_hours,
                data = df)
  
  
  df$log_theta_old <- log(1-df$theta_poisson)
  res_aov_log_theta_old <- aov(log_theta_old ~ sample_time_hours,
                data = df)
  
  return(data.frame(sample_time_hours=df$sample_time_hours,
                    residuals=res_aov$residuals,
                    residuals_log=res_aov_log$residuals,
                    residuals_log_theta_old=res_aov_log_theta_old$residuals))
})))

theta_anova_resid$sample_time_hours <- factor(theta_anova_resid$sample_time_hours)

plt <- rasterize(ggqqplot(theta_anova_resid, 
                          x="residuals", color="sample_time_hours") + 
                   facet_wrap(vars(sample_time_hours)) +
                   scale_color_manual(values=accent_colors[c(2, 4, 7, 9)]), dpi=300)

ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_theta_anova_residuals_qqplot.pdf"), 
       height=100, width=100, units="mm")

plt <- rasterize(ggqqplot(theta_anova_resid, 
                          x="residuals_log", color="sample_time_hours") + 
                   facet_wrap(vars(sample_time_hours)) +
                   scale_color_manual(values=accent_colors[c(2, 4, 7, 9)]), dpi=300)

ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_theta_anova_residuals_log_qqplot.pdf"), 
       height=100, width=100, units="mm")


plt <- rasterize(ggqqplot(theta_anova_resid, 
                          x="residuals_log_theta_old", color="sample_time_hours") + 
                   facet_wrap(vars(sample_time_hours)) +
                   scale_color_manual(values=accent_colors[c(2, 4, 7, 9)]), dpi=300)

ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_theta_anova_residuals_log_theta_old_qqplot.pdf"), 
       height=100, width=100, units="mm")


# kruskal tests
gene_theta_kw <- data.frame(do.call(rbind,
                                    lapply(gene_thetas, function(df) {
  res_kw <- kruskal.test(theta_poisson ~ sample_time_hours,
                data = df)
  return(data.frame(statistic=res_kw$statistic,
                    pvalue=res_kw$p.value))

} )))

plt <- rasterize(ggplot(gene_theta_kw, aes(pvalue)) + 
  geom_histogram(bins=100), dpi=300) + 
  theme_classic() +
  xlab("Kruskal-Wallis P-Value\nComparing NTR Over Time") + 
  ylab("Gene Count") +
  theme(text=element_text(size=7))
ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_theta_kruskal_wallis.pdf"),
       height=50, width=50, units="mm")

# now look at distribution of median thetas per gene over time
mean_gene_theta <- gene_thetas_df %>% 
  group_by(geneid, sample_time_hours) %>%
  summarize(mean_theta_poisson=mean(theta_poisson)) 

ggplot(mean_gene_theta, 
       aes(sample_time_hours, mean_theta_poisson, 
           group=sample_time_hours)) + 
  geom_violin(linewidth=0.1, draw_quantiles = 0.5, fill=accent_colors[2]) +
  theme_classic() +
  scale_x_continuous(breaks=c(2, 4, 6, 8)) 
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_mean_gene_theta_violin.pdf"),
       height=50, width=50, units="mm")

# correlation
gene_theta_cor <- data.frame(do.call(rbind,
                                     lapply(gene_thetas,
                                             
                         function(df) {
                           cor_res <- cor.test(df$sample_time_hours,
                                               df$theta_poisson, 
                                          method = "spearman")
                           return(data.frame(statistic=cor_res$statistic,
                                             rho=cor_res$estimate,
                                             pvalue=cor_res$p.value
                                             ))
                         }
                                     )
)
)

plt <- ggplot(gene_theta_cor, aes(rho)) + 
  geom_histogram(bins=100) +
  theme_classic() +
  theme(text=element_text(size=7))

ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_test_gene_correlation.pdf"),
       height=50, width=50, units="mm")



## ---------------------------------------------------------------------------------------------------------------

# plot the residuals
residual_cols <- grepl("X[0-9]+h", colnames(gene_kinetics_df))
no_scale_residuals <- gene_kinetics_df[,residual_cols] %>%
  pivot_longer(contains("X")) 
no_scale_residuals$model_type <- "no_scale"
kw_no_scale_resid <- kruskal.test(value~name, data=no_scale_residuals)

scale_residuals <- gene_kinetics_theta_correction_df[,residual_cols] %>% 
  pivot_longer(contains("X"))
scale_residuals$model_type <- "global_scale"
kw_scale_resid <- kruskal.test(value~name, data=scale_residuals)

# genewise is much worse
scale_residuals_genewise <- gene_kinetics_theta_correction_genewise_df[,residual_cols] %>% 
  pivot_longer(contains("X"))
scale_residuals_genewise$model_type <- "gene_scale"
kw_scale_genewise_resid <- kruskal.test(value~name, data=scale_residuals_genewise)

data.frame(rbind(no_scale_residuals, scale_residuals)) %>%
  ggplot(aes(model_type, value, fill=name)) + 
  geom_violin(draw_quantiles = 0.5, linewidth=0.1) + 
  scale_fill_manual(values=accent_colors[c(10, 8, 4,2)]) +
  theme_classic() + 
  ylab("Mean Residual: Expected-Observed NTR") +
  theme(legend.position = "none", text=element_text(size=7)) +
  geom_hline(yintercept=0, linetype="dashed") 
ggsave(file.path(plot_subdir, "202308_halflife_fit_residuals_scale_vs_no_scale.pdf"),
       height=30, width=50, units="mm")

## for global model plot rsquared
ggplot(gene_kinetics_df %>% inner_join(gene_kinetics_theta_correction_df, by=c("geneid"), suffix=c(".no_scale", ".global_scale")), aes(rsquared_nls.no_scale, rsquared_nls.global_scale)) + 
  geom_point(size=0.1, color='black') +
  geom_density_2d(linewidth=0.1, color='white') +
  xlim(-2, 1) + 
  ylim(-2, 1) +
  geom_abline(linewidth=0.1, linetype='dashed') + 
  theme_classic() +
  theme(text=element_text(size=7))
ggsave(file.path(plot_subdir, "202308_rsquared_expfit_scale_vs_no_scale.pdf"),
       height=30, width=30, units="mm")


## ---------------------------------------------------------------------------------------------------------------
# rename df and cols
gene_kinetics_df_final <- gene_kinetics_theta_correction_df
gene_kinetics_df_final$peak_name <- gene_kinetics_df_final$geneid

# Calculate halflife as log(2)/kdeg
gene_kinetics_df_final$halflife_Veeraraghavan_CPN <- log(2)/gene_kinetics_df_final$k_nls

# plot kdeg vs r-squared of fit to kinetic model
rasterize(ggplot(gene_kinetics_df_final, aes(k_nls, rsquared_nls)) + 
  geom_point(size=0.1), dpi=300) + 
  geom_density_2d(color=background_colors[2], linewidth=0.1) + 
  scale_x_log10() +
  theme_classic() + 
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_kdeg_vs_rsquared_kinetic_model.pdf"),
       height=50, width=50, units="mm")

# plot halflife vs r-squared of fit to kinetic model (basically same as above)
rasterize(ggplot(gene_kinetics_df_final, aes(halflife_Veeraraghavan_CPN, rsquared_nls)) + 
  geom_point(size=0.1), dpi=300) + 
  geom_density_2d(color=background_colors[2], linewidth=0.1) + 
  scale_x_log10() +
  ylim(-2, 1) +
  theme_classic() + 
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_vs_rsquared_kinetic_model.pdf"),
       height=50, width=50, units="mm")

# plot rsquared histogram
ggplot(gene_kinetics_df_final, aes(rsquared_nls)) + 
  geom_histogram(color="black", fill='black', bins=100) + 
  xlim(-2, 1) +
  theme_classic() + 
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_rsquared_kinetic_model_histogram.pdf"),
       height=50, width=50, units="mm")


# plot rsquared histogram, filtering out the genes for which the model is not informative
ggplot(gene_kinetics_df_final %>% subset(rsquared_nls > 0.1), aes(halflife_Veeraraghavan_CPN )) + 
  geom_histogram(color="black", fill='black',bins=100) + 
  scale_x_log10() +
   theme_classic() + 
  theme(text=element_text(size=7))
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_histogram_rsquared_geq10pct.pdf"),
       height=50, width=50, units="mm")

# define and plot half-life quantiles 
quantiles <- quantile(subset(gene_kinetics_df_final, rsquared_nls > 0.1)$halflife_Veeraraghavan_CPN, 
                      seq(0, 1, 0.2))[2:4]

plt <- ggplot(gene_kinetics_df_final %>% subset(rsquared_nls > 0.1), aes(halflife_Veeraraghavan_CPN )) + 
  geom_histogram(color=background_colors[4], fill=background_colors[4],bins=100) + 
  scale_x_log10() +
   theme_classic() + 
  theme(text=element_text(size=7))

for (i in quantiles) {
  plt <- plt + geom_vline(xintercept=i, color="black", linewidth=0.2)
}
ggsave(plt, file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_histogram_rsquared_geq10pct_quantiles.pdf"),
       height=50, width=50, units="mm")


## ---------------------------------------------------------------------------------------------------------------
# load the quantseq annotations
load(file.path(file_path_ref$project_directory,
               file_path_ref$results$quantseq_3end_gene_annotated))

cpn_gene_rel_label_annot <- unique(quantseq_utr_annotated_filt[c("peak_name", "utr_name")]) %>%
  inner_join(gene_kinetics_df_final %>% subset(rsquared_nls > 0.1), 
             by=c( "peak_name")) %>%
  dplyr::rename(ensembl_transcript_id=utr_name)

# get gene ids
gene_map <- getBM(c( "external_gene_name", "external_transcript_name", "ensembl_gene_id", "ensembl_transcript_id"),
                  filters="ensembl_transcript_id",
                  values=unique(cpn_gene_rel_label_annot$ensembl_transcript_id),
                  mart=grcm38_101)

cpn_halflife_filt <- cpn_gene_rel_label_annot %>% inner_join(gene_map, by="ensembl_transcript_id")


## ---------------------------------------------------------------------------------------------------------------

genes_to_plot <- list(rad21="peak_5224", bmyc="peak_9244", Tuba1a="peak_5672", Cdk5r1="peak_2404", Gnas="peak_9079", Ybx1="peak_10833", Fus="peak_14041", Rhoa="peak_15816", Rplp0="peak_11941", satb2="peak_5295", bcl11a="peak_1282", gpm6a="peak_15356", actb="peak_12290")


plot_example_gene <- function(peak_name, 
                          gene_name="",
                          theta_df=gene_thetas_df, 
                          kinetic_df=gene_kinetics_df_final) {
  
  true_data <- subset(theta_df, geneid == peak_name)
  true_data$sample_time_hours <- true_data$sample_time_hours
  stime <- 0:9
  k_nls <- kinetic_df$k_nls[kinetic_df$geneid == peak_name]
  C <- kinetic_df$C_nls[kinetic_df$geneid == peak_name]
  rsquared_linear <- kinetic_df$rsquared_linear[kinetic_df$geneid == peak_name]
  rsquared_nls <- kinetic_df$rsquared_nls[kinetic_df$geneid == peak_name]

  predict_data <- data.frame(sample_time_hours=rep(stime, 2),
                             theta_poisson=c((1-C*exp(-k_nls*(stime)))))
                                          
  plt <- ggplot(true_data, aes(sample_time_hours, theta_poisson)) + 
    geom_point(size=0.1) + 
    geom_line(data=predict_data, 
              mapping=aes(sample_time_hours, theta_poisson),
              linewidth=0.2) +
    theme_classic() +
    ylim(0, 1) +
    theme(text=element_text(size=7, color='black')) +
    scale_x_continuous(breaks=c(2, 4, 6, 8)) +
    ggtitle(sprintf("%s %s: halflife %s, rsquared %s", gene_name, peak_name, log(2)/k_nls,  rsquared_nls))
  return(plt)
}


# slow: Ybx1
# mid: ActB
# fast: RhoA, Cdk5r1

plot_example_gene(genes_to_plot$Cdk5r1, "Cdk5r1")
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_exemplar_short_cdk5r1.pdf"),
       height=45, width=45, units="mm")


plot_example_gene(genes_to_plot$Rhoa, "RhoA")
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_exemplar_short_rhoa.pdf"),
       height=45, width=45, units="mm")

plot_example_gene(genes_to_plot$actb, "Actb")
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_exemplar_mid_actb.pdf"),
       height=45, width=45, units="mm")


plot_example_gene(genes_to_plot$Ybx1, "Ybx1")
ggsave(file=file.path(plot_subdir, "202308_SLAM_Soma_halflife_exemplar_long_ybx1.pdf"),
       height=45, width=45, units="mm")





## ---------------------------------------------------------------------------------------------------------------

# plot exemplar classes
get_ensembl_go <- function(goterms) {
  
  ens_genes <- c()
  for (goterm in c(goterms)) {
  genes <- (AnnotationDbi::select(org.Mm.eg.db,  keytype="GOALL", keys=goterm, columns = "ENSEMBL") %>%
    dplyr::select(ENSEMBL) %>%
    unique())$ENSEMBL
  ens_genes <- c(ens_genes, genes)
  }
  return(ens_genes)
  
}
ribosomal_genes <- subset(cpn_halflife_filt, 
                          grepl("Rpl|Rps|Rack", external_gene_name)) %>%
  rowwise() %>%
  mutate(gene_group="Ribosomal Protein")

tfs <- subset(cpn_halflife_filt, 
              ensembl_gene_id %in% get_ensembl_go("GO:0003700") ) %>%
  rowwise() %>%
  mutate(gene_group="Transcription Factor")

tubulin_assoc <- subset(cpn_halflife_filt, 
                    ensembl_gene_id %in% get_ensembl_go("GO:0015630")) %>%
  rowwise() %>%
  mutate(gene_group="Microtubule Cytoskeleton")

sig_transduc <- subset(cpn_halflife_filt, 
                 ensembl_gene_id %in% get_ensembl_go("GO:0009966")) %>%
  rowwise() %>%
  mutate(gene_group="Regulation Signal Transduction")

mito <- subset(cpn_halflife_filt, 
               ensembl_gene_id %in% get_ensembl_go("GO:0006119")) %>% #GO:0005739")) %>%
  rowwise() %>%
  mutate(gene_group="Mitochondria (OxPhos)")

translation <- subset(cpn_halflife_filt, 
               ensembl_gene_id %in% get_ensembl_go(c("GO:0006415", "GO:0006414"))) %>%
  rowwise() %>%
  mutate(gene_group="Translation (Elongation & Termination)")

proteasome <- subset(cpn_halflife_filt, 
                     ensembl_gene_id %in% get_ensembl_go("GO:0000502")) %>%
    rowwise() %>%
  mutate(gene_group="Proteasome")

nucleosome <- subset(cpn_halflife_filt, 
                     ensembl_gene_id %in% get_ensembl_go("GO:0006334")) %>%
  rowwise() %>%
  mutate(gene_group="Nucleosome")

kinase <- subset(cpn_halflife_filt, 
                 ensembl_gene_id %in% get_ensembl_go("GO:0016301")) %>%
  rowwise() %>%
  mutate(gene_group="Kinase")

axon_guidance <- subset(cpn_halflife_filt, 
                        ensembl_gene_id %in% get_ensembl_go("GO:0007411"))%>%
  rowwise() %>%
  mutate(gene_group="Axon Guidance")

chromatin_remodel <- subset(cpn_halflife_filt, 
                        ensembl_gene_id %in% get_ensembl_go("GO:0006338"))%>%
  rowwise() %>%
  mutate(gene_group="Chromatin Remodeling")
  


all_genes <- cpn_halflife_filt %>%
                      rowwise() %>% mutate(gene_group="All")
to_plot <- data.frame(do.call(rbind,
                              list(all_genes, chromatin_remodel, kinase, ribosomal_genes, proteasome, tfs, translation, mito, sig_transduc)))
to_plot$gene_group <- factor(to_plot$gene_group, levels=c("All", "Transcription Factor", "Chromatin Remodeling", "Regulation Signal Transduction", "Kinase", "Proteasome", "Translation (Elongation & Termination)","Mitochondria (OxPhos)",  "Ribosomal Protein"))
to_plot %>%
  ggplot(aes(gene_group, halflife_Veeraraghavan_CPN, fill=gene_group)) + 
  geom_violin(draw_quantiles = 0.5, linewidth=0.1) + 
  geom_hline(yintercept=median(cpn_halflife_filt$halflife_Veeraraghavan_CPN, na.rm=TRUE),
             linetype="dashed", color=background_colors[7], linewidth=0.2) + 
  scale_y_log10() +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values=c(background_colors[5], accent_colors[3:10])) +
  theme(legend.position="None", text=element_text(size=5)) +
  xlab("Gene Set") +
  ylab("Estimated Halflife (h)")

ggsave(file.path(file_path_ref$project_directory,
                 plot_subdir, "gene_sets_halflife_distribution.pdf"),
       height=50, width=75, units="mm")



## ---------------------------------------------------------------------------------------------------------------
rna_dynamics_dir <- file.path(file_path_ref$project_directory, "data/external_data/RNA_dynamics")

ext_data_meta <- readxl::read_xlsx(file.path(rna_dynamics_dir, 
"RNA_dynamics_comparator_datasets.xlsx"), sheet="unique_from_Agarwal_Kelley")

# get mouse mappings
gene_map1 <- getBM(c("affy_mouse430_2", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id"),
                  filters="ensembl_transcript_id",
                  values=unique(cpn_gene_rel_label_annot$ensembl_transcript_id),
                  mart=grcm38_101)
# get the human mapping from a different page
gene_map2 <- getBM(c("hsapiens_homolog_associated_gene_name", "ensembl_transcript_id"),
                   filters="ensembl_transcript_id",
                  values=unique(cpn_gene_rel_label_annot$ensembl_transcript_id),
                  mart=grcm38_101)

gene_annot_info <- gene_map1 %>% inner_join(gene_map2, by="ensembl_transcript_id")



# read in data from other studies
read_halflife_df <- function(meta_row) {
  df <- readxl::read_xlsx(file.path(rna_dynamics_dir, 
                                    sprintf("%s.xlsx", meta_row$filename)), 
                          sheet=meta_row$sheet_name)
  
  # rename columns
  desired_cols <- c(meta_row$identifier_type, sprintf("halflife_est_%s", meta_row$measurement_name))
  setnames(df, old=c(meta_row$gene_id_column, meta_row$halflife_column),
           new=desired_cols) 
  
  df <- df[,desired_cols] %>% inner_join(gene_annot_info)
  return(df)
}


# load the pre-compiled comparisons from Agarwal and Kelley
ext_halflife_df <- purrr::reduce(lapply(1:nrow(ext_data_meta), 
                                        function(i)
                     read_halflife_df(ext_data_meta[i,])),
              function(x, y) x %>% full_join(y))

ext_halflife_df2 <- read.table(file.path(rna_dynamics_dir, "all_HLs_mouse.txt"), 
                              header = TRUE) %>%
  pivot_longer(cols=-c(external_gene_name, ensembl_gene_id),
               names_to = "study_name",
               values_to = "halflife") %>%
  rowwise() %>%
  mutate(study_name=sprintf("halflife_%s", gsub("_[0-9]+$", "", study_name))) %>%
  group_by(ensembl_gene_id, external_gene_name, study_name) %>%
  summarize(halflife=mean(halflife)) %>%
  pivot_wider(id_cols=c(ensembl_gene_id, external_gene_name), names_from = study_name, values_from=halflife)


# now compare to cpn
all_halflife_df <- ext_halflife_df %>%
  full_join(ext_halflife_df2) %>%
  full_join(cpn_gene_rel_label_annot)

#%>%
#  rowwise() %>%
#  mutate(rel_order_cpn=1/estimate)

# there are a ton of duplicate entries
# so go pairwise and get unique values
# and compute spearman correlation between samples
halflife_cols <- colnames(all_halflife_df)[grepl("halflife", colnames(all_halflife_df))]
n <- length(halflife_cols)
cor_mat <- matrix(nrow=n, ncol=n)
cosine_sim_mat <- matrix(nrow=n, ncol=n)

#plots <- list()
for (i in 1:n) {
  for (j in 1:n) {
    unique_vals <- unique(all_halflife_df[,c(halflife_cols[i], halflife_cols[j])])
    
    dist_i <- as.numeric(unlist(unique_vals[,1]))
    dist_j <- as.numeric(unlist(unique_vals[,2]))
    
    
    cor_mat[i,j] <- cor(dist_i, dist_j, method="spearman", use="complete.obs" ) 
    
    norm_quant <- as.matrix(data.frame(i=dist_i, j=dist_j))
    norm_quant <- log(norm_quant+0.5)
    norm_quant <- na.omit(norm_quant)
    cosine_sim_mat[i,j] <- cosine(norm_quant[,1], norm_quant[,2])
    
  }
} 

colnames(cor_mat) <- halflife_cols
rownames(cor_mat) <- halflife_cols

colnames(cosine_sim_mat) <- halflife_cols
rownames(cosine_sim_mat) <- halflife_cols

pdf(file.path(plot_subdir, "supplement_cpn_halflife_vs_4su_all.blackwhite.pdf"),
    height=10, width=10)
pheatmap(cor_mat, clustering_distance_rows = "correlation", 
          clustering_distance_cols = "correlation",
          color =colorRampPalette(background_colors[c(1, 2, 3, 6:10)])(50),
         fontsize = 5)
dev.off()




## ---------------------------------------------------------------------------------------------------------------
read_count_df <- cpn_conversion_rates %>% 
  pivot_wider(id_cols=Geneid, 
              names_from=sample_id, 
              values_from=reads_total) 

meta_soma_nonzero <- subset(meta_soma, sample_time_min > 0) %>%
  rowwise() %>%
  mutate(sample_time_hours=sample_time_min/60)

counts_data <- data.frame(read_count_df[,meta_soma_nonzero$sample_id])
rownames(counts_data) <- read_count_df$Geneid

dds <- DESeqDataSetFromMatrix(counts_data,
                       meta_soma_nonzero,
                       ~sample_time_hours)

dds <- DESeq(dds)

res <- results(dds, name="sample_time_hours")
res_df <- data.frame(res)
res_df$peak_name <- read_count_df$Geneid

res_df %>%
  subset(peak_name %in% valid_cpn_genes) %>% 
  ggplot(aes(log2FoldChange, -log10(padj))) + 
  geom_point(size=0.1) +
  theme_classic() +
  theme(text=element_text(size=7))
ggsave(file.path(file_path_ref$project_directory,
                 plot_subdir, "202308_SLAM_Soma_read_count_dge.pdf"),
       height=30, width=30, units="mm")


## ---------------------------------------------------------------------------------------------------------------
total_sample_t_content <- cpn_conversion_rates %>% group_by(sample_id) %>% summarize(total_t=sum(T_count_total), total_tc=sum(TC_count_total))
t_cpm <- cpn_conversion_rates %>% inner_join(total_sample_t_content, by="sample_id") %>% 
  rowwise() %>%
  mutate(t_cpm=T_count_total/total_t*1e6,
         tc_cpm=TC_count_total/total_t*1e-6) %>%
  group_by(Geneid) %>%
  summarize(mean_t_cpm=mean(t_cpm),
            mean_tc_cpm=mean(tc_cpm))

rasterize(cpn_halflife_filt %>% 
  inner_join(t_cpm, by=c("peak_name" ="Geneid" )) %>% 
  ggplot(aes(mean_t_cpm, halflife_Veeraraghavan_CPN)) + geom_point(size=0.1, alpha=0.5) +
  scale_x_log10() + 
  scale_y_log10(), dpi=300) + 
  theme(text=element_text(size=5)) +
  stat_cor(method='spearman') + theme_classic()

ggsave(file.path(plot_subdir, "supplement_fig1_t_content_correlation_halflife_estimate.pdf"),
       height=60, width=60, units="mm")


## ---------------------------------------------------------------------------------------------------------------
n <- 10000
expect_new <- data.frame(pois=table(rpois(n, 0.0125*20))/n) %>% 
  dplyr::rename(n_conv=pois.Var1) %>%
  mutate(read_type="new")
expect_old <- data.frame(pois=table(rpois(n, 0.0005708237*20))/n) %>%
  dplyr::rename(n_conv=pois.Var1) %>%
  mutate(read_type="old")


rbind(expect_new, expect_old, data.frame(n_conv=setdiff(expect_new$n_conv, expect_old$n_conv),
                                         pois.Freq=0,
                                         read_type="old")) %>%
  ggplot(aes(n_conv, pois.Freq, fill=read_type)) + 
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  scale_fill_manual(values=c(accent_colors[2], background_colors[4])) + 
  theme_classic() +
  theme(text=element_text(size=7))
ggsave(file.path(plot_subdir, "202308_supplement_fig1_simulation_readcount_poisson.pdf"),
       height=120, width=60, units="mm")




## ---------------------------------------------------------------------------------------------------------------
save(cpn_halflife_filt, all_halflife_df, 
     gene_thetas_df, gene_kinetics_df_final, 
     file=file.path(file_path_ref$project_directory,
                    file_path_ref$results$soma_halflife_data))

