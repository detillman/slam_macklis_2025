## Helper scripts
## SLAM setup script
library(bakR)
library(biomaRt)
library(boot)
library(rjson)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)

input <- "Y:/users/priyav/slam_final/metadata.json"

file_path_ref <- fromJSON(file = input)

metadata <- read.csv(file.path(file_path_ref$project_directory, file_path_ref$design))

grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")

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

common_genes <- list(gpm6a="peak_15356", actb="peak_12290")
