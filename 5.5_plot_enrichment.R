rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)

# Paths
# my_lib <- "/home/caf77/R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))

basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
# basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
# setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir <- file.path(basedir, "index")

summary_dir <- file.path(basedir, "summary_tables")

figpath <- file.path(basedir, "figures")
dir.create(figpath,showWarnings = F,recursive = T)

waffle_file="all_samples_internal_merged_indexed_filtered.tsv"
B_genes_non_singleton_count="Cargo_genes_pairs_multisite_bySampleID.tsv"
Unique_non_singleton_count="Unique_LGT_pairs_multisite_bySampleID.tsv"

site_col <- c(
  "humanFecal"="#56B4E9", 
  "humanSkin"="#EE2C2C", 
  "humanOral"="#EFC000", 
  "dogFecal"="#3A5FCD", 
  "dogSkin"="#7D0226", 
  "dogOral"="#9A5324"  
)

pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")
type.ord <- c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")

#######################
## Helper Functions  ##
#######################
save_plot <- function(p, name, width, height, unit="in", scale=1) {
  formats <- c("rds", "png", "pdf")
  for(fmt in formats){
    file <- file.path(figpath, paste0(name, ".", fmt))
    if(fmt == "rds") saveRDS(p, file)
    else ggsave(file, plot = p, width = width, height = height, dpi = 150, unit=unit, scale=scale) 
  }
}

################
## Load Data  ##
################
B_genes <- read.delim(file.path(summary_dir,B_genes_non_singleton_count))
LGT <- read.delim(file.path(summary_dir,Unique_non_singleton_count))

pop.df <- read.delim(file.path(input_dir, "Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(sampleid, 1, 6)) %>%
  select(SAMPLE, population, individual, type, site) %>%
  distinct(SAMPLE, .keep_all = TRUE) %>%
  merge.data.frame(data.frame(SAMPLE=substr(colnames(B_genes)[-1],1,6),
                              original_sid=colnames(B_genes)[-1]), by="SAMPLE") %>%
  mutate(
    group = case_when(
      population %in% c("BTU", "ORT", "ORS") ~ "Early-transition",
      population %in% c("APT", "TBU", "BSP") ~ "Late-transition",
      TRUE ~ "Agriculture"
    ),
    host = ifelse(grepl("x", individual), "dog", "human"),
    type = case_when(
      type == "DHF" ~ "humanFecal",
      type == "DHV" ~ "humanOral",
      type == "DHP" ~ "humanSkin",
      type == "DDF" ~ "dogFecal",
      type == "DDV" ~ "dogOral",
      type == "DDH" ~ "dogSkin"
    )
  ) %>%
  arrange(host, population) %>%
  mutate(
    individual = factor(individual, levels = unique(individual)),
    population = factor(population, levels = pop.ord),
    group = factor(group, levels = c("Early-transition", "Late-transition", "Agriculture")),
    type = factor(type, levels = type.ord),
    host = factor(host, levels = c("human", "dog")),
    site = factor(site, levels = c("Fecal", "Oral", "Skin"))
  )


################
## Load Data  ##
################

