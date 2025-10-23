## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggpubr)
library(PerformanceAnalytics)
library(DESeq2)
library(ggrepel)
library(readxl)
library(pheatmap)
library(tidyverse)

background_colors <- colorRampPalette(c("white", "black"))(10)
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


## ----utrSelection-----------------------------------------------------------------------------------------------
justify <- read.csv("Expt255_utrs_rbps_halflives.csv") %>%
  mutate(utr_half = paste0(UTR, ": ", Half.Life, " Hrs"))

ggplot(justify, aes(x = Half.Life, y = Hits)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~ RBP, scales = "free")

ggplot(justify, aes(x = RBP, y = Hits)) +
  geom_col(pos = "dodge") +
  facet_wrap(~ utr_half, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_y_continuous(trans = "log10")


## ----dnaResults-------------------------------------------------------------------------------------------------
dna_counts <- read.csv("Sequencing/20240519_all_barcode_counts.csv") %>%
  select(-X) %>%
  separate("max_prob_reference_group", 
           into = c("plasmid", "exp", "number", "utr"), remove = F) %>%
  filter(sample_id == "DNA_input") %>%
  mutate(first = gsub("\\[", "", reads_reference_group)) %>%
  mutate(second = gsub("\\]", "", first)) %>%
  separate(second, into = c("a", "b", "c", "d", "e", "f", 
                            "g", "h", "i", "j", "k", "l"), 
           sep = "\\s+") %>%
  pivot_longer(cols = c("a", "b", "c", "d", "e", "f", 
                            "g", "h", "i", "j", "k", "l"),
               names_to = "position",
               values_to = "reads") %>%
  mutate(reads = as.numeric(reads))

# loading dna counts
dna_counts <- read.csv("Sequencing/20240519_all_barcode_counts.csv") %>%
  select(-X) %>%
  separate("max_prob_reference_group", 
           into = c("plasmid", "exp", "number", "utr"), remove = F) %>%
  filter(sample_id == "DNA_input") %>%
  mutate(first = gsub("\\[", "", reads_reference_group)) %>%
  mutate(second = gsub("\\]", "", first)) %>%
  separate(second, into = c("a", "b", "c", "d", "e", "f", 
                            "g", "h", "i", "j", "k", "l"), 
           sep = "\\s+") %>%
  pivot_longer(cols = c("a", "b", "c", "d", "e", "f", 
                            "g", "h", "i", "j", "k", "l"),
               names_to = "position",
               values_to = "reads") %>%
  mutate(reads = as.numeric(reads)) %>%
  group_by(barcode_trimmed, umi, bc2, min_edit_distance, 
           prop_reads_max_ref_group, max_prob_reference_group, plasmid, exp, 
           number, utr, sample_id) %>%
  summarize(reads = max(reads, na.rm = T)) %>%
  ungroup()
  

# ggplot(dna_counts, aes(x = as.factor(min_edit_distance), y = umi)) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log10") +
#   theme_classic(base_size = 10)
# ggplot(dna_counts, aes(x = prop_reads_max_ref_group)) +
#   geom_histogram() +
#   theme_classic(base_size = 10)
utr_barcodes <- ggplot(dna_counts, aes(x = umi, 
                                       y = prop_reads_max_ref_group)) +
  geom_point(size = 0.1, color  = background_colors[5]) +
  #geom_density_2d(linewidth = 0.1, color = "white") +
  # stat_regline_equation(label.y = 0.95,
  #                       label.x = 4.5,
  #                       size = 2,
  #                       aes(label = after_stat(rr.label))) +
  # geom_smooth(method='lm', formula = 'y ~ x', 
  #             color = accent_colors[2],
  #             linetype = "dashed",
  #             linewidth = 0.1,
  #             se = F) +
  scale_x_continuous(trans = "log10") +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "Number of UMIs",
       y = "Mapped Reads per Barcode (%)"); utr_barcodes
ggsave("utr_barcodes.pdf", utr_barcodes, units = "mm",
       height = 50, width = 50)

# dna_counts_grouped <- dna_counts %>%
#   unite("barcode_plus_umi", "barcode_trimmed", "umi", remove = F) %>%
#   group_by(max_prob_reference_group, plasmid, exp, number, utr, 
#            sample_id, reads) %>%
#   summarize(barcode_umi = max(umi)) %>%
#   ungroup()
# ggplot(dna_counts_grouped, aes(x = reads, y = barcode_umi, color = utr)) +
#   geom_point() +
#   scale_y_continuous(trans = "log10") +
#   scale_x_continuous(trans = "log10") +
#   #facet_wrap(~ utr) +
#   theme_classic(base_size = 10)

dna_counts_umis <- dna_counts %>%
  group_by(utr) %>%
  summarize(sum_umi = sum(umi))

utr_dna <- ggplot(dna_counts_umis %>%
                    mutate(utr = fct_relevel(utr, "Auts2", "Marcks", "Gpm6a", 
                                              "Actb", "Mapt", "Gnas", "Gad1", 
                                              "Ldlrad3")), 
                  aes(x = utr, y = sum_umi, fill = utr)) +
  geom_bar(stat = "summary", position = "dodge") +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values=c("#f46d43",
                             "#fdae61",
                             "#fee090",
                             "#e0f3f8",
                             "#abd9e9",
                             "#74add1",
                             "grey83",
                             "grey67")) +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_hline(yintercept =  0, color = "black", linewidth = 0.25) +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "",
       y = "Number of UMIs"); utr_dna
ggsave("utr_dna.pdf", utr_dna, units = "mm",
       height = 50, width = 50)


## ----dataLoad---------------------------------------------------------------------------------------------------
# loading priya's soma half-lives
print(load("Soma_Kinetics/Soma_Kinetics.RData"))
priya_soma <- load("Soma_Kinetics/Soma_Kinetics.RData")

# loading library counts
library_counts <- read.csv("Sequencing/20240519_all_barcode_counts.csv") %>%
  select(-X) %>%
  separate("max_prob_reference_group", 
           into = c("plasmid", "exp", "number", "utr")) %>%
  filter(sample_id != "DNA_input") %>%
  group_by(sample_id, utr) %>%
  summarize(umi_sum = sum(umi)) %>%
  rename(sample = sample_id) %>%
  ungroup()

# loading non-library counts
all_counts <- read.csv("Sequencing/20240519_non_library_sample_counts.csv") %>%
  select(sample_id, reference, num_umi) %>%
  rename(sample = sample_id,
         construct = reference) %>%
  mutate(construct = gsub("_cds", "", construct)) %>%
  separate(sample, into = c("region", "time", "replicate"), remove = F) %>%
  pivot_wider(names_from = "construct", values_from = "num_umi") %>%
  left_join(library_counts %>%
              group_by(sample) %>%
              summarize(library = sum(umi_sum)) %>%
              ungroup(), by = "sample") %>%
  mutate(library = case_when(is.na(library) == T ~ 0, T ~ library))


## ----allCounts--------------------------------------------------------------------------------------------------
# lengthening all counts for plots
all_counts_long <- all_counts %>%
  pivot_longer(cols = c("fluc", "rluc", "mclover", "library"),
               names_to = "construct", values_to = "num_umi") %>%
  mutate(num_umi = case_when(num_umi == 0 ~ 0.1, T ~ num_umi))

# axon samples have very few umis
ggplot(all_counts_long %>%
         unite(col = "plot", "region", "time"), 
       aes(x = plot, y = num_umi, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Number of UMIs") +
  scale_y_continuous(trans = "log10")

# removing axon samples and both luciferase constructs
all_ipsi <- all_counts %>%
  filter(region == "ipsi") %>%
  select(-contains("luc"))
all_ipsi_long <- all_counts_long %>%
  filter(region == "ipsi") %>%
  filter(construct == "mclover" | construct == "library")

# mclover > library
ggplot(all_ipsi_long %>%
         unite(col = "plot", "region", "time"), 
       aes(x = plot, y = num_umi, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Number of UMIs") +
  scale_y_continuous(trans = "log10")

utr_normBar <- ggplot(all_ipsi_long %>%
                        mutate(time = paste0(time, "h"),
                               construct = case_when(construct == "library" ~ "Library",
                                                     construct == "mclover" ~ "mClover")), 
                      aes(x = construct, y = num_umi, fill = time)) +
  scale_fill_manual(values = c(background_colors[5], accent_colors[2])) +
  geom_boxplot(outliers = F, linewidth = 0.1) +
  # geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  # geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
  #               width = 0.25) +
  geom_point(position = position_dodge(width = 0.75), size = 0.25) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  labs(x = "",
       y = "Number of UMIs",
       fill = "") +
  scale_y_continuous(trans = "log10"); utr_normBar
ggsave("utr_normBar.pdf", utr_normBar, units = "mm", height = 50, width = 50)

# mclover and library counts are strongly correlated
ggplot(all_ipsi, 
       aes(x = mclover, y = library)) +
  geom_point(aes(color = time)) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 200, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 180,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "mClover UMIs",
       y = "Library UMIs")

utr_normScatter <- ggplot(all_ipsi %>%
                            mutate(time = paste0(time, "h")), 
       aes(x = mclover, y = library)) +
  geom_point(aes(color = time), size = 0.25) +
  scale_color_manual(values = c(background_colors[5], accent_colors[2])) +
  geom_smooth(method='lm', formula = 'y ~ x', 
              color = "black",
              linetype = "dashed",
              linewidth = 0.1,
              se = F) +
  # stat_regline_equation(label.y = 200, 
  #                       aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 300,
                        aes(label = after_stat(rr.label)),
                        size = 2) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  labs(x = "mClover UMIs",
       y = "Library UMIs", color = ""); utr_normScatter
ggsave("utr_normScatter.pdf", utr_normScatter, units = "mm",
       height = 50, width = 50)


# removing ipsi samples and both luciferase constructs
all_axon <- all_counts %>%
  filter(region == "axon") %>%
  select(-contains("luc"))
all_axon_long <- all_counts_long %>%
  filter(region == "axon") %>%
  filter(construct == "mclover" | construct == "library")

# mclover > library
ggplot(all_axon_long %>%
         unite(col = "plot", "region", "time"), 
       aes(x = plot, y = num_umi, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Number of UMIs") +
  scale_y_continuous(trans = "log10")

utr_normBar_axon <- ggplot(all_axon_long %>%
                        mutate(time = paste0(time, "h"),
                               construct = case_when(construct == "library" ~ "Library",
                                                     construct == "mclover" ~ "mClover")), 
                      aes(x = construct, y = num_umi, fill = time)) +
  scale_fill_manual(values = c(background_colors[5], accent_colors[2])) +
  geom_boxplot(outliers = F, linewidth = 0.1) +
  # geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  # geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
  #               width = 0.25) +
  geom_point(position = position_dodge(width = 0.75), size = 0.25) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  labs(x = "",
       y = "Number of UMIs",
       fill = "") +
  scale_y_continuous(trans = "log10"); utr_normBar_axon
ggsave("utr_normBar_axon.pdf", utr_normBar_axon, units = "mm", 
       height = 50, width = 50)

# mclover and library counts are strongly correlated
ggplot(all_axon, 
       aes(x = mclover, y = library)) +
  geom_point(aes(color = time)) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 200, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 180,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "mClover UMIs",
       y = "Library UMIs")

utr_normScatter_axon <- ggplot(all_axon %>%
                            mutate(time = paste0(time, "h")), 
       aes(x = mclover, y = library)) +
  geom_point(aes(color = time), size = 0.25) +
  scale_color_manual(values = c(background_colors[5], accent_colors[2])) +
  geom_smooth(method='lm', formula = 'y ~ x', 
              color = "black",
              linetype = "dashed",
              linewidth = 0.1,
              se = F) +
  # stat_regline_equation(label.y = 200, 
  #                       aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 4,
                        aes(label = after_stat(rr.label)),
                        size = 2) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top") +
  labs(x = "mClover UMIs",
       y = "Library UMIs", color = ""); utr_normScatter_axon
ggsave("utr_normScatter_axon.pdf", utr_normScatter_axon, units = "mm",
       height = 50, width = 50)


## ----libraryAnalyze---------------------------------------------------------------------------------------------
# normalization factors for non-UTR counts
all_counts_norm <- all_counts %>%
  filter(region == "ipsi") %>%
  mutate(mclover_norm = mclover / mean(mclover),
         library_norm = library / mean(library),
         all_norm = (mclover + library) / 
           (mean(mclover) + mean(library)))

# normalizing library counts
library_utrs <- library_counts %>%
  left_join(all_counts_norm, by = "sample") %>%
  filter(region == "ipsi") %>%
  mutate(umi_mcloverNorm = umi_sum / mclover_norm,
         umi_libNorm = umi_sum / library_norm,
         umi_allNorm = umi_sum / all_norm)

# correlations between utr umis and non-utr umis
chart.Correlation(library_utrs %>% select(umi_sum, mclover, library))

# correlations after different normalization approaches
chart.Correlation(library_utrs %>% select(umi_sum, mclover, library,
                                          contains("Norm", ignore.case = F)))


ggplot(library_utrs, 
       aes(x = utr, y = umi_mcloverNorm, fill = time)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  stat_compare_means(aes(group = time), label = "p.format", size = 3) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "UMIs (mClover Normalized)")
ggplot(library_utrs, 
       aes(x = utr, y = umi_libNorm, fill = time)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  stat_compare_means(aes(group = time), label = "p.format", size = 3) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "UMIs (Library Normalized)")
ggplot(library_utrs, 
       aes(x = utr, y = umi_allNorm, fill = time)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  stat_compare_means(aes(group = time), label = "p.format", size = 3) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "UMIs (mClover/Library Normalized)")

# correlating rna turnover with soma half-lives
library_utrs_halfLives <- library_utrs %>%
  left_join(cpn_halflife_filt %>%
              rename(utr = external_gene_name) %>%
              distinct(utr, halflife_Veeraraghavan_CPN),
            by = "utr") %>%
  select(sample, utr, contains("Norm", ignore.case = F), 
           halflife_Veeraraghavan_CPN) %>%
  pivot_wider(names_from = "sample",
              values_from = contains("Norm")) %>%
  mutate(umi_mcloverNorm_ipsi_1_ratio = umi_mcloverNorm_ipsi_8_1 / umi_mcloverNorm_ipsi_3_1,
         umi_mcloverNorm_ipsi_2_ratio = umi_mcloverNorm_ipsi_8_2 / umi_mcloverNorm_ipsi_3_2,
         umi_mcloverNorm_ipsi_3_ratio = umi_mcloverNorm_ipsi_8_3 / umi_mcloverNorm_ipsi_3_3,
         umi_mcloverNorm_ipsi_4_ratio = umi_mcloverNorm_ipsi_8_4 / umi_mcloverNorm_ipsi_3_4,
         umi_mcloverNorm_ipsi_5_ratio = umi_mcloverNorm_ipsi_8_5 / umi_mcloverNorm_ipsi_3_5,
         
         umi_libNorm_ipsi_1_ratio = umi_libNorm_ipsi_8_1 / umi_libNorm_ipsi_3_1,
         umi_libNorm_ipsi_2_ratio = umi_libNorm_ipsi_8_2 / umi_libNorm_ipsi_3_2,
         umi_libNorm_ipsi_3_ratio = umi_libNorm_ipsi_8_3 / umi_libNorm_ipsi_3_3,
         umi_libNorm_ipsi_4_ratio = umi_libNorm_ipsi_8_4 / umi_libNorm_ipsi_3_4,
         umi_libNorm_ipsi_5_ratio = umi_libNorm_ipsi_8_5 / umi_libNorm_ipsi_3_5,
         
         umi_allNorm_ipsi_1_ratio = umi_allNorm_ipsi_8_1 / umi_allNorm_ipsi_3_1,
         umi_allNorm_ipsi_2_ratio = umi_allNorm_ipsi_8_2 / umi_allNorm_ipsi_3_2,
         umi_allNorm_ipsi_3_ratio = umi_allNorm_ipsi_8_3 / umi_allNorm_ipsi_3_3,
         umi_allNorm_ipsi_4_ratio = umi_allNorm_ipsi_8_4 / umi_allNorm_ipsi_3_4,
         umi_allNorm_ipsi_5_ratio = umi_allNorm_ipsi_8_5 / umi_allNorm_ipsi_3_5) %>%
  select(halflife_Veeraraghavan_CPN, utr, contains("ratio")) %>%
  pivot_longer(cols = contains("ratio"), names_to = "sample", values_to = "ratio") %>%
  separate(sample, into = c("umi", "metric", "region", "replicate", "drop")) %>%
  distinct(halflife_Veeraraghavan_CPN, utr, metric, replicate, ratio)
  
ggplot(library_utrs_halfLives %>%
         filter(metric == "mcloverNorm"), 
       aes(x = halflife_Veeraraghavan_CPN, y = ratio)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 5, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 4,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (mClover Normalized)") 
ggplot(library_utrs_halfLives %>%
         filter(metric == "libNorm"), 
       aes(x = halflife_Veeraraghavan_CPN, y = ratio)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 5, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 4,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (Library Normalized)") 
ggplot(library_utrs_halfLives %>%
         filter(metric == "allNorm"), 
       aes(x = halflife_Veeraraghavan_CPN, y = ratio)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 5, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 4,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (mClover/Library Normalized)") 


## ----deseqNormalize---------------------------------------------------------------------------------------------
# making counts matrix
deseq_counts <- library_counts %>%
  pivot_wider(names_from = "sample", values_from = "umi_sum") %>%
  replace(is.na(.), 0) %>%
  column_to_rownames("utr") %>%
  select(contains("ipsi"))

# making metadta matrix
deseq_metadata <- library_counts %>%
  distinct(sample) %>%
  separate(sample, into = c("region", "time", "replicate"), remove = F) %>%
  filter(region == "ipsi") %>%
  column_to_rownames("sample")

# double-checking counts and metadata matrices
all(rownames(deseq_metadata) == colnames(deseq_counts))

# running DESeq2
deseq_dds <- DESeqDataSetFromMatrix(countData = deseq_counts,
                                    colData = deseq_metadata,
                                    design = ~ time)
deseq_dds$time <- relevel(deseq_dds$time, ref = "3")
deseq <- DESeq(deseq_dds)
deseq_results <- results(deseq, name = "time_8_vs_3") %>% 
  as.data.frame()

# plotting DESeq2 results without normalization
ggplot(deseq_results %>%
         rownames_to_column("utr") %>%
         pivot_longer(cols = contains("p"), 
                      names_to = "p_type", values_to = "p_stat"), 
       aes(x = log2FoldChange, y = -log10(p_stat), label = utr)) +
  facet_wrap(~ p_type, scales = "free") +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Log2(8 Hours - 3 Hours)",
       y = "-log10(p statistic)",
       title = "No Normalization") +
  geom_text_repel()
ggplot(deseq_results %>%
         rownames_to_column("utr") %>%
         left_join(cpn_halflife_filt %>%
                     rename(utr = external_gene_name) %>%
                     distinct(utr, halflife_Veeraraghavan_CPN),
                   by = "utr"), 
       aes(x = halflife_Veeraraghavan_CPN, y = log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 1, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 0.5,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (No Normalization") 

# creating mclover normalization matrix
mclover_norm <- library_utrs %>%
  select(sample, utr, mclover) %>%
  pivot_wider(names_from = "sample", values_from = "mclover") %>%
  mutate(ipsi_3_2 = 712, ipsi_3_4 = 848) %>%
  column_to_rownames("utr") %>%
  as.matrix()
mclover_norm <- mclover_norm / exp(rowMeans(log(mclover_norm)))

# running DESeq2 with mclover normalization matrix
deseq_dds_mclover <- deseq_dds
normalizationFactors(deseq_dds_mclover) <- mclover_norm
deseq_mclover <- DESeq(deseq_dds_mclover)
deseq_results_mclover <- results(deseq_mclover, name = "time_8_vs_3") %>% 
  as.data.frame()

# plotting DESeq2 results with mclover normalization
ggplot(deseq_results_mclover %>%
         rownames_to_column("utr") %>%
         pivot_longer(cols = contains("p"), 
                      names_to = "p_type", values_to = "p_stat"), 
       aes(x = log2FoldChange, y = -log10(p_stat), label = utr)) +
  facet_wrap(~ p_type, scales = "free") +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Log2(8 Hours - 3 Hours)",
       y = "-log10(p statistic)",
       title = "mClover Normalization") +
  geom_text_repel()
ggplot(deseq_results_mclover %>%
         rownames_to_column("utr") %>%
         left_join(cpn_halflife_filt %>%
                     rename(utr = external_gene_name) %>%
                     distinct(utr, halflife_Veeraraghavan_CPN),
                   by = "utr"), 
       aes(x = halflife_Veeraraghavan_CPN, y = log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 1, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 0.5,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (mClover Normalization") 

# creating library normalization matrix
library_norm <- library_utrs %>%
  select(sample, utr, library) %>%
  pivot_wider(names_from = "sample", values_from = "library") %>%
  mutate(ipsi_3_2 = 46, ipsi_3_4 = 54) %>%
  column_to_rownames("utr") %>%
  as.matrix()
library_norm <- library_norm / exp(rowMeans(log(library_norm)))

# running DESeq2 with library normalization matrix
deseq_dds_library <- deseq_dds
normalizationFactors(deseq_dds_library) <- library_norm
deseq_library <- DESeq(deseq_dds_library)
deseq_results_library <- results(deseq_library, name = "time_8_vs_3") %>% 
  as.data.frame()

# plotting DESeq2 results with library normalization
ggplot(deseq_results_library %>%
         rownames_to_column("utr") %>%
         pivot_longer(cols = contains("p"), 
                      names_to = "p_type", values_to = "p_stat"), 
       aes(x = log2FoldChange, y = -log10(p_stat), label = utr)) +
  facet_wrap(~ p_type, scales = "free") +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Log2(8 Hours - 3 Hours)",
       y = "-log10(p statistic)",
       title = "Library Normalization") +
  geom_text_repel()
ggplot(deseq_results_library %>%
         rownames_to_column("utr") %>%
         left_join(cpn_halflife_filt %>%
                     rename(utr = external_gene_name) %>%
                     distinct(utr, halflife_Veeraraghavan_CPN),
                   by = "utr"), 
       aes(x = halflife_Veeraraghavan_CPN, y = log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 1, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 0.5,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (Library Normalization") 

# creating mclover/library normalization matrix
all_norm_calc <- mclover_norm %>%
  as.data.frame() %>%
  pivot_longer(cols = contains("ipsi"), 
               names_to = "sample", values_to = "mclover_norm") %>%
  distinct(sample, mclover_norm) %>%
  full_join(library_norm %>%
              as.data.frame() %>%
              pivot_longer(cols = contains("ipsi"),
                           names_to = "sample", values_to = "library_norm") %>%
              distinct(sample, library_norm),
            by = "sample") %>%
  mutate(all_norm = (mclover_norm + library_norm) / 2)
all_norm <- library_utrs %>%
  select(-contains("norm")) %>%
  left_join(all_norm_calc, by = "sample") %>%
  select(sample, utr, all_norm) %>%
  pivot_wider(names_from = "sample", values_from = "all_norm") %>%
  mutate(ipsi_3_2 = 0.3258201, ipsi_3_4 = 0.3854183) %>%
  column_to_rownames("utr") %>%
  as.matrix()
all_norm <- all_norm / exp(rowMeans(log(all_norm)))

# running DESeq2 with mclover/library normalization matrix
deseq_dds_all <- deseq_dds
normalizationFactors(deseq_dds_all) <- all_norm
deseq_all <- DESeq(deseq_dds_all)
deseq_results_all <- results(deseq_all, name = "time_8_vs_3") %>% 
  as.data.frame()

# normalizng with rlog
all_rld_blind <- rlog(deseq_all, blind = T)

# pca for top 500 genes
all_pca_var <- apply(assay(all_rld_blind), 1, sd)
all_pca_var_df <- assay(all_rld_blind)[order(all_pca_var, 
                       decreasing = TRUE)[seq_len(8)],]
all_pca <- prcomp(t(all_pca_var_df), scale = FALSE)
all_pca_df <- all_pca$x %>% 
  data.frame() %>% 
  rownames_to_column("sample") %>% 
  separate("sample", into = c("region", "time", "replicate"), remove = F)
all_pca_percent <- round(100 * all_pca$sdev^2/sum(all_pca$sdev^2), 1)

# plotting PC1 and PC2 for top 8 genes
utr_pca <- ggplot(all_pca_df %>%
         mutate(Time = case_when(time == 3 ~ "3h",
                                 time == 8 ~ "8h")), 
       aes(get(paste0("PC", 1)), 
                   get(paste0("PC", 2)),
                   col = Time)) + 
  labs(x = paste0("PC", 1, ": ", all_pca_percent[1], "%"), 
       y = paste0("PC", 2, ": ", all_pca_percent[2], "%"),
       color = "") + 
  coord_fixed() +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_point(size = 1) +
  scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "top"); utr_pca
ggsave("utr_pca.pdf", utr_pca, units = "mm",
       height = 50, width = 50)

# heatmap of the top 8 most variable genes
all_select <- order(rowVars(counts(deseq_all, normalized=TRUE)),
                decreasing=TRUE)[1:8]
all_select_df <- as.data.frame(colData(deseq_all)[,c("time")])

pheatmap(assay(all_rld_blind)[order(all_pca_var, 
                       decreasing = TRUE)[seq_len(8)],], 
         cluster_rows = T, show_rownames = T, 
         show_colnames = T, cluster_cols = T, 
         annotation_col = as.data.frame(colData(deseq_all)[,c("region", "time", "replicate")]))

# plotting DESeq2 results with mclover/library normalization
ggplot(deseq_results_all %>%
         rownames_to_column("utr") %>%
         pivot_longer(cols = contains("p"), 
                      names_to = "p_type", values_to = "p_stat"), 
       aes(x = log2FoldChange, y = -log10(p_stat), label = utr)) +
  facet_wrap(~ p_type, scales = "free") +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Log2(8 Hours - 3 Hours)",
       y = "-log10(p statistic)",
       title = "mClover/Library Normalization") +
  geom_text_repel()
ggplot(deseq_results_all %>%
         rownames_to_column("utr") %>%
         left_join(cpn_halflife_filt %>%
                     rename(utr = external_gene_name) %>%
                     distinct(utr, halflife_Veeraraghavan_CPN),
                   by = "utr"), 
       aes(x = halflife_Veeraraghavan_CPN, y = log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 1, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 0.5,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Priya Half Life",
       y = "UMI Ratio (mClover/Library Normalization") 

utr_diff <- ggplot(deseq_results_all %>%
                     rownames_to_column("utr") %>%
                     mutate(label_up = case_when(log2FoldChange > 0 ~ 
                                                    paste0("q = ", as.character(signif(padj, digits = 2))),
                                                  T ~ ""),
                            label_down = case_when(log2FoldChange < 0 ~ 
                                                    paste0("q = ", as.character(signif(padj, digits = 2))),
                                                  T ~ ""),
                            utr = fct_relevel(utr, "Auts2", "Marcks", "Gpm6a", 
                                              "Actb", "Mapt", "Gnas", "Gad1", 
                                              "Ldlrad3")), 
       aes(x = utr, y = log2FoldChange, fill = utr)) +
  geom_bar(stat = "summary", position = "dodge") +
  scale_fill_manual(values=c("#f46d43",
                             "#fdae61",
                             "#fee090",
                             "#e0f3f8",
                             "#abd9e9",
                             "#74add1",
                             "grey83",
                             "grey67")) +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_hline(yintercept =  0, color = "black", linewidth = 0.25) +
  geom_text(aes(label = label_up),
            position = position_dodge(width = .9), vjust = -0.5, size = 3 / .pt) +
  geom_text(aes(label = label_down),
            position = position_dodge(width = .9), vjust = 1.2, size = 3 / .pt) +
  geom_errorbar(aes(ymin = log2FoldChange-lfcSE, ymax = log2FoldChange+lfcSE),
                linewidth = 0.25, width = 0.25) +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "",
       y = "log2FoldChange (8h - 3h)"); utr_diff
ggsave("utr_diff.pdf", utr_diff, units = "mm",
       height = 50, width = 50)


## ----halfLifeComparison-----------------------------------------------------------------------------------------
# loading half-lives from Loedige et al, Mol Cell, 2023
# loedige_comp <- read_excel("Loedige_MolCell.xlsx",
#                            sheet = "Half-lives compartments")
loedige_whole <- read_excel("Loedige_MolCell.xlsx", 
                            sheet = "halflife_wt_whole_neuron")

# loading priya's soma half-lives
print(load("Soma_Kinetics/Soma_Kinetics.RData"))
priya_soma <- load("Soma_Kinetics/Soma_Kinetics.RData")

# getting y values for half-life histogram
halfLife_hist <- ggplot(cpn_halflife_filt, 
       aes(x = halflife_Veeraraghavan_CPN)) +
  geom_histogram(bins = 100); halfLife_hist
halfLife_hist_values <- ggplot_build(halfLife_hist)$data[[1]]

utr_halfLife_hist <- data.frame(utr_name = c("Marcks", "Marcks", "Marcks", "Mapt",
                                        "Mapt", "Mapt", "Gnas", "Auts2", 
                                        "Actb", "Gpm6a"),
                           half_life = c(10.287260, 10.539929, 6.444936, 
                                         16.577309, 17.804996, 17.822574,
                                         21.441875, 5.558173, 13.568570,
                                         10.958915),
                           hist_y = c(62, 68, 187,
                                         6, 7, 7,
                                         1, 96, 17,
                                         35))
utr_noHalfLife_hist <- data.frame(utr_name = c("Gad1", "Ldlrad3"),
                                  fake_half_x = c(15, 15),
                                  fake_hist_y = c(200, 175))

utr_halfLife <- ggplot(cpn_halflife_filt, 
                       aes(x = halflife_Veeraraghavan_CPN)) +
  geom_histogram(bins = 100) +
  geom_point(data = utr_halfLife_hist,
             size = 0.25,
             aes(x = half_life, y = hist_y),
             color = accent_colors[2]) +
  geom_label_repel(data = utr_halfLife_hist,
                   aes(x = half_life, y = hist_y, label = utr_name, fill = utr_name),
                   ylim = c(50, max(halfLife_hist_values$y)),
                   size = 1) +
  geom_point(data = utr_noHalfLife_hist, 
             aes(x = fake_half_x, y = fake_hist_y),
             color = "white") +
  geom_label_repel(data = utr_noHalfLife_hist,
                   aes(x = fake_half_x, y = fake_hist_y, label = utr_name, 
                       fill = utr_name),
                   #ylim = c(175, max(halfLife_hist_values$y)),
                   size = 1) +
  #scale_fill_brewer(palette = "Set1", direction = 1) +
  scale_fill_manual(values=c("#e0f3f8",
                             "#f46d43",
                             "grey83",
                             "#74add1",
                             "#fee090",
                             "grey67",
                             "#abd9e9",
                             "#fdae61")) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "Relative Halflife in Soma (h)",
       y = "Gene Count"); utr_halfLife
ggsave("utr_halfLife.pdf", utr_halfLife, units = "mm",
       height = 50, width = 50)

# dataframe to compare priya and loedige half lives
priya_loedige <- all_halflife_df %>%
  select(external_gene_name, ensembl_gene_id, 
         ensembl_transcript_id, hsapiens_homolog_associated_gene_name,
         halflife_est_Loedige_whole_PCN_SLAM, halflife_Veeraraghavan_CPN) %>%
  distinct() %>%
  mutate(loedige = case_when(halflife_est_Loedige_whole_PCN_SLAM == ">16" ~ "16",
                             T ~ halflife_est_Loedige_whole_PCN_SLAM),
         priya = halflife_Veeraraghavan_CPN,
         .keep = "unused") %>%
  mutate(loedige = as.numeric(loedige)) %>%
  filter(is.na(loedige) == F | is.na(priya) == F) %>%
  mutate(is_reporter = case_when(external_gene_name == "Actb" | # moderate half-life, CThPN gcs but enriched in somata
           external_gene_name == "Auts2" | # short half-life, enriched in CPN and CThPN gcs over somata, contains gc-rich sequence
           external_gene_name == "Gad1" | # long half-life, enriched in CPN gcs over somata
           external_gene_name == "Gnas" | # second longest half-life, loedige half-life is > 16
           external_gene_name == "Gpm6a" | # short half-life, in gcs?, no loedige half-life
           external_gene_name == "Ldlrad3" | # longest half-life, in gcs?
           external_gene_name == "Mapt" |
           external_gene_name == "Marcks" ~ T, # shortest overall half-life, in CPN GCs
           T ~ F),
         priya_over_loedige = priya / loedige,
         utr_name = case_when(is_reporter == T ~ external_gene_name)) %>%
  arrange(is_reporter)

# comparing priya and loedige half lives
ggplot(priya_loedige, aes(x = loedige,
                          y = priya)) +
  geom_point(aes(color = is_reporter)) +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 20, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 18,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Loedige Half Life",
       y = "Priya Half Life") +
  geom_label_repel(data = priya_loedige %>%
                    filter(is_reporter == T),
                  aes(label = external_gene_name, color = is_reporter))

utr_loedige <- ggplot(priya_loedige, aes(x = loedige,
                          y = priya)) +
  geom_point(aes(color = is_reporter), size = 0.1) +
  geom_density_2d(linewidth = 0.1, color = "white") +
  geom_smooth(method='lm', formula = 'y ~ x', 
              color = accent_colors[2],
              linetype = "dashed",
              linewidth = 0.1,
              se = F) +
  stat_regline_equation(label.y = 22,
                        size = 2,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  
  labs(x = "Loedige et al. Half Life (h)",
       y = "Veeraraghavan et al. Half Life (h)") +
  geom_label_repel(data = priya_loedige %>%
                    filter(is_reporter == T),
                   min.segment.length = 0.1,
                   size = 1,
                  aes(label = external_gene_name, 
                      fill = utr_name)) +
  scale_color_manual(values=c(background_colors[5], accent_colors[2])) +
  scale_fill_manual(values=c("#e0f3f8",
                             "#f46d43",
                             "#74add1",
                             "#fee090",
                             "#abd9e9",
                             "#fdae61"),
                     na.value = "grey50"); utr_loedige
ggsave("utr_loedige.pdf", utr_loedige, units = "mm",
       height = 50, width = 50)

# filtering for reporter utrs
priya_loedige_reporters <- priya_loedige %>%
  filter(is_reporter == T) %>%
  mutate(priya_expect = 4.1 + 0.93 * loedige,
         priya_resid = priya - priya_expect)

# comparing loedige and priya half-lives for reporter utrs
ggplot(priya_loedige_reporters, aes(x = loedige,
                          y = priya, 
                          label = external_gene_name)) +
  geom_point() +
  geom_smooth(method='lm', formula = 'y ~ x') +
  stat_regline_equation(label.y = 20, 
                        aes(label = after_stat(eq.label))) +
  stat_regline_equation(label.y = 18,
                        aes(label = after_stat(rr.label))) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Loedige Half Life",
       y = "Priya Half Life") +
  geom_text_repel()

# long dataframe for reporter half-life analysis
priya_loedige_reporters_long <- priya_loedige_reporters %>%
  pivot_longer(cols = c("loedige", "priya"), 
               names_to = "experiment",
               values_to = "half_life")

# reporters have longer priya half-lives than loedige half-lives
ggplot(priya_loedige_reporters_long, 
       aes(x = external_gene_name, y = half_life, fill = experiment)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Half-Life",
       fill = "Experiment")

# Gnas is the reporter with most similar half-life between priya and loedige
ggplot(priya_loedige_reporters_long %>%
         select(-experiment, -half_life) %>%
         distinct(), 
       aes(x = external_gene_name, y = priya_over_loedige)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Half-Life Ratio (Priya / Loedige)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",
             linewidth = 0.5)

# Auts2 and Marcks have least deviation from overall line of best fit
ggplot(priya_loedige_reporters_long %>%
         select(-experiment, -half_life) %>%
         distinct(), 
       aes(x = external_gene_name, y = priya_resid)) +
  geom_bar(stat = "summary", position = "dodge", fun = "mean") +
  geom_errorbar(stat='summary', position = position_dodge(width = 0.9), 
                width = 0.25) +
  geom_point(position = position_dodge(width = 0.9)) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", angle = 90, 
                                   vjust = 0.5, hjust=1),
        axis.text.y = element_text(color = "black")) +
  labs(x = "",
       y = "Half-Life Residual") +
  geom_hline(yintercept = 0, color = "black")


## ----sessionInfo------------------------------------------------------------------------------------------------
sessionInfo()

