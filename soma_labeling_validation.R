## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(dplyr)
library(ggplot2)

source("C:/Users/prv892/Documents/GitHub/slam/scripts/helper_script.R")



## ---------------------------------------------------------------------------------------------------------------
plot_supp1c <- function(plot_savefile) {
  load(file.path(file_path_ref$project_directory, 
               file_path_ref$results$soma_counts_data))

cpn_conversion_rates_all <- counts_soma %>%
  rowwise() %>%
  mutate(conversion_rate=TC_count_total/T_count_total,
         sample_time_hours=sample_time_min/60)

cpn_global_conversion_rate <- cpn_conversion_rates_all %>%
  group_by(sample_id, sample_time_hours, IUE_plasmid) %>%
  summarize(conversion_rate=sum(TC_count_total)/sum(T_count_total)) %>%
  subset(sample_time_hours %in% c(0, 2, 4, 6, 8)) %>%
  rowwise() %>%
  mutate(sample_group=paste0(sample_time_hours, IUE_plasmid))

ggplot(cpn_global_conversion_rate, aes(sample_time_hours, conversion_rate, 
                                       group=sample_group,
                                       fill=IUE_plasmid)) + 
  geom_boxplot(linewidth=0.1) + geom_point(position=position_dodge(1.6), size=0.2) + 
  theme_classic() + 
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8)) +
  scale_fill_manual(values=c(background_colors[5], accent_colors[2])) +
  theme(text=element_text(size=7))
ggsave(file=plot_savefile,
       height=30, width=60, units="mm")
}

plot_supp1c(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/Supplement1_conversion_rate_SLAM_soma_whole_sample.pdf"))


## ---------------------------------------------------------------------------------------------------------------
plot_supp1d <- function(plot_savefile) {
  conversion_rates_injections <- read.csv(file.path(file_path_ref$project_directory, 
                                                    "data/sequencing/Expt146_Test_Injections_Amplicon/conversion_rates_injections.csv"))
  conversion_rates_injections$Injection <- factor(conversion_rates_injections$Injection, levels=c("Subcutaneous", "Intraperitoneal"))
  
  conversion_rates_injections %>%
    ggplot(aes(Injection, conversion_rate, fill=Enzyme)) +
    geom_boxplot(linewidth=0.1, outlier.shape=NA) +
    geom_point(position=position_dodge(0.8), size=0.2) +
    theme_classic() +
    scale_fill_manual(values=c(background_colors[5], accent_colors[2])) +
    theme(text=element_text(size=7))
  ggsave(file=plot_savefile,
       height=30, width=60, units="mm")
}

plot_supp1d(file.path(file_path_ref$project_directory, "plots/illustrator_pdfs/Supplement1_conversion_rate_SLAM_soma_Injection_Strategy.pdf"))


