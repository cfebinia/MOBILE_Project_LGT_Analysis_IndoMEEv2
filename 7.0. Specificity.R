rm(list=ls())
gc()
set.seed(7)

library(dplyr)
library(tidyr)
library(stringr)


# Paths
# my_lib <- "/home/caf77/R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))

basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra/"
# basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
# basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir <- file.path(basedir, "index")

summary_dir <- file.path(basedir, "summary_tables")

figpath <- file.path(basedir, "figures")
dir.create(figpath,showWarnings = F,recursive = T)

waffle_file="all_samples_internal_merged_indexed_filtered.tsv"
outlier="outliers_waafle.txt"

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
types <- c("humanFecal","humanOral","humanSkin","dogFecal","dogOral","dogSkin")
#######################
## Helper Functions  ##
#######################
save_plot <- function(p, name, width, height, unit="in", dpi=150, scale=1) {
  formats <- c("rds", "png", "pdf")
  for(fmt in formats){
    file <- file.path(figpath, paste0(name, ".", fmt))
    if(fmt == "rds") saveRDS(p, file)
    else ggsave(file, plot = p, width = width, height = height, dpi = dpi, unit=unit, scale=scale) 
  }
}


################
## Load Data  ##
################
outlier_ids <- read.delim(file.path(input_dir,outlier))$sampleid

indata <- read.delim(file.path("summary_tables",waffle_file))
indata <- indata %>% 
  mutate(SYNTENY_LEN = nchar(SYNTENY)) %>%
  filter(SYNTENY_LEN <= 35 & B_count <= 10) %>%
  filter(!original_sid %in% outlier_ids)

pop.df <- read.delim(file.path(input_dir, "Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(sampleid, 1, 6)) %>%
  select(SAMPLE, population, individual, type, site) %>%
  distinct(SAMPLE, .keep_all = TRUE) %>%
  merge.data.frame(data.frame(SAMPLE=substr(indata$original_sid,1,6), 
                              new_sid=indata$SAMPLE) %>% 
                     distinct(), 
                   by = "SAMPLE") %>%
  select(-SAMPLE) %>%
  rename(SAMPLE=new_sid) %>%
  distinct() %>%
  mutate(
    group = case_when(
      population %in% c("BTU", "ORT", "ORS") ~ "Early-transition",
      population %in% c("APT", "TBU", "BSP") ~ "Late-transition",
      TRUE ~ "Agriculture"),
    host = ifelse(grepl("x", individual), "dog", "human")) %>%
  arrange(host, population) %>%
  mutate(
    individual = factor(individual, levels = unique(individual)),
    population = factor(population, levels = pop.ord),
    group = factor(group, levels = c("Early-transition", "Late-transition", "Agriculture")),
    type = factor(type, levels = type.ord),
    host = factor(host, levels = c("human", "dog")),
    site = factor(site, levels = c("Fecal", "Oral", "Skin"))
  )



####################
# Host specificity #
####################
library(UpSetR)
library(ggplot2)

dim(indata)
length(unique(indata$SAMPLE))

df <- indata %>%
  mutate(SAMPLE=substr(original_sid,1,6)) %>%
  mutate(sandwichDist = paste("len",A_gapLen,sandwich, sep = "_")) %>%
  select(sandwichDist, SAMPLE, original_sid, CONTIG_NAME, DIRECTION, LCA) %>%
  mutate(pair_contig=paste(original_sid, CONTIG_NAME, sep="+")) %>%
  mutate(type = substr(SAMPLE,1,3)) %>%
  mutate(type = factor(type, levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH"))) %>%
  left_join(pop.df %>% 
              mutate(SAMPLE=substr(SAMPLE,1,6)) %>%
              select(SAMPLE, individual) %>% 
              distinct(), by="SAMPLE")

indata %>%
  mutate(type = substr(SAMPLE,1,3)) %>%
  mutate(type = factor(type, levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH"))) %>%
  group_by(type) %>%
  summarise(n=n_distinct(SAMPLE))
  
sandwich_nSample <- df %>%
  select(sandwichDist, type, SAMPLE, original_sid, individual) %>%
  distinct() %>%
  group_by(sandwichDist) %>%
  summarise(
    humanFecal = sum(grepl("DHF", SAMPLE)),
    humanOral = sum(grepl("DHV", SAMPLE)),
    humanSkin = sum(grepl("DHP", SAMPLE)),
    dogFecal = sum(grepl("DDF", SAMPLE)),
    dogOral = sum(grepl("DDV", SAMPLE)),
    dogSkin = sum(grepl("DDH", SAMPLE)),
    human =  sum(grepl("DHF|DHV|DHP", individual)),
    dog = sum(grepl("DDF|DDV|DDH", individual)),
    across_type = n_distinct(individual),
  .groups = "drop") %>%
  arrange(desc(across_type)) %>%
  select(-across_type)

head(sandwich_nSample)


upset_data <- sandwich_nSample %>%
  mutate(across(c(humanFecal, humanOral, humanSkin, dogFecal, dogOral, dogSkin), ~ ifelse(.x > 0, 1, 0))) 

my_upset <- upset(
  as.data.frame(upset_data),
  sets = rev(c("humanFecal", "humanOral", "humanSkin", "dogFecal", "dogOral", "dogSkin")),
  keep.order = TRUE,
  main.bar.color = "steelblue",
  sets.bar.color = rev(site_col[types]),
  point.size = 2,
  line.size = 0.5,
  shade.alpha = 0,
  scale.sets = 1,
  text.scale = 1,
  set_size.show = TRUE,
  set_size.numbers_size = 8,
  set_size.scale_max = max(colSums(upset_data[,-1]))*1.5,
  mainbar.y.label = "# Unique LGT",
  sets.x.label = "# Unique LGT by Type",
  mb.ratio = c(0.65, 0.35)
)
my_upset 

# 1. Save as PDF
figout <- file.path(figpath,"upset/upset_plot")
pdf(file = paste0(figout, ".pdf"), width = 15 / 2.54, height = 10 / 2.54)
my_upset
dev.off()

# 2. Save as PNG
png(
  filename = paste0(figout, ".png"),
  width = 15 / 2.54,
  height = 10 / 2.54,
  units = "in",
  res = 600
)
my_upset
dev.off()

# 3. Save as RDS
saveRDS(my_upset, file = paste0(figout, ".rds"))

##
# by site
###
library(gridExtra)
library(grid)
df1 <- df %>%
  left_join(
    pop.df %>%
      select(individual, population) %>%
      distinct(),
    by = "individual"
  )

for (i in unique(df1$type)) {
  df2 <- df1 %>%
    filter(type == i)
  
  upset_data <- df2 %>%
    select(sandwichDist, population, SAMPLE, original_sid, individual) %>%
    distinct() %>%
    group_by(sandwichDist, population) %>%
    summarise(
      n = as.numeric(n_distinct(individual) > 0),
      .groups = "drop"
    ) %>%
    pivot_wider(id_cols = sandwichDist, names_from = population, values_from = n, values_fill = 0) %>%
    arrange(across(all_of(unique(df2$population)), desc))
  
  myset <- pop.ord[pop.ord %in% colnames(upset_data)]
  
  my_upset <- upset(
    as.data.frame(upset_data),
    sets = rev(myset),
    keep.order = TRUE,
    main.bar.color = "steelblue",
    point.size = 2,
    line.size = 0.5,
    shade.alpha = 0,
    scale.sets = 1,
    text.scale = 1,
    set_size.show = TRUE,
    set_size.numbers_size = 8,
    set_size.scale_max = max(colSums(upset_data[,-1]))*1.5,
    mainbar.y.label = "# Distinct LGT",
    sets.x.label = paste("# Distinct LGT in", i),
    mb.ratio = c(0.6, 0.4)
  )
  
  print(my_upset)
  
  count <- df2 %>%
    filter(type == i) %>%
    group_by(population) %>%
    summarise(n = n_distinct(individual), .groups = "drop")
  
  print(i)
  print(count)
  
  figout <- file.path(figpath, paste0("upset/upset_plot_byPop_", i))
  
  i_width = ifelse(i %in% c("DHF", "DHV","DDF"),
                   25, 20) / 2.54
  i_height = ifelse(i %in% c("DHF", "DHV","DDF"),
                   10, 8) / 2.54
  pdf(file = paste0(figout, ".pdf"), width = i_width, height =i_height)
  print(my_upset)
  dev.off()
  
  png(
    filename = paste0(figout, ".png"),
    width = i_width,
    height = i_height,
    units = "in",
    res = 1200
  )
  print(my_upset)
  dev.off()
  
  saveRDS(my_upset, file = paste0(figout, ".rds"))
}

dev.off()


sink(file.path(figpath,"upset/n_sample.txt"))
print("ALL")
indata %>%
  mutate(type = substr(SAMPLE,1,3)) %>%
  mutate(type = factor(type, levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH"))) %>%
  group_by(type) %>%
  summarise(n=n_distinct(SAMPLE))

for (i in unique(df1$type)) {
  df2 <- df1 %>%
    filter(type == i)
  print(i)
  count <- df2 %>%
    filter(type == i) %>%
    group_by(population) %>%
    summarise(n = n_distinct(individual), .groups = "drop")
  print(count)
}
sink()
