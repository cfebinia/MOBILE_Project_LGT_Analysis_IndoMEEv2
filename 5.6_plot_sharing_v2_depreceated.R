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
LGT_count="nLGT_by_seqID.tsv" #"nLGT_by_sample.tsv" #
dist_df="dist/ruzicka_uniqueLGT_multisite_bySampleID.rds"
LGTmat="dist/ruzicka_uniqueLGT_multisite_bySampleID_LGTmatrix.rds"
pairs_count="Unique_LGT_pairs_multisite_bySampleID.tsv"
outlier="outliers_waafle.txt"
sPairs="all_seqID_pairs_waafle.tsv"

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
outlier_ids <- read.delim(file.path(input_dir,outlier))$sampleid

count_df <- read.delim(file.path(summary_dir,LGT_count)) %>%
  filter(!SAMPLE %in% outlier_ids) %>%
  mutate(original_sid=gsub("_","",original_sid)) %>%
  mutate(type=case_when(
    type=="DHF" ~ "humanFecal",
    type=="DHV" ~ "humanOral",
    type=="DHP" ~ "humanSkin",
    type=="DDF" ~ "dogFecal",
    type=="DDV" ~ "dogOral",
    type=="DDH" ~ "dogSkin"
  )) %>%
  mutate(population=factor(population, pop.ord),
         group=factor(group,c("Early-transition","Late-transition","Agriculture")),
         type=factor(type, type.ord),
         host=factor(host, c("human","dog")),
         site=factor(site, c("Fecal","Oral","Skin"))) %>%
  group_by(type) %>%
  mutate(is_outlier = {
    range <- quantile(count_unique, c(0.025, 0.975))
    count_unique < range[1] | count_unique > range[2]
  }) %>%
  ungroup() %>%
  distinct()

ru_dist <- readRDS(file.path(basedir, dist_df)) %>% as.matrix()
colnames(ru_dist) <- gsub("_","",colnames(ru_dist))
rownames(ru_dist) <- gsub("_","",rownames(ru_dist))
original_sid <- colnames(ru_dist)

pop.df <- read.delim(file.path(input_dir, "Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(sampleid, 1, 6)) %>%
  select(SAMPLE, population, individual, type, site) %>%
  distinct(SAMPLE, .keep_all = TRUE) %>%
  merge.data.frame(data.frame(SAMPLE = substr(original_sid, 1, 6),
                              original_sid = original_sid), by = "SAMPLE") %>%
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
  ) %>%
  distinct()

pairs_df <- read.delim(file.path(input_dir, sPairs)) %>%
  mutate(type=gsub("human|dog","",type)) %>%
  mutate(seqID_1=gsub("_","",seqID_1),
         seqID_2=gsub("_","",seqID_2),
         seqID1_2=paste(seqID_1,seqID_2,sep="_")) %>%
  distinct()

LGTmat <- readRDS(LGTmat)
rownames(LGTmat) <- gsub("_","",rownames(LGTmat))

####################
## SAMPLE FILTER  ##
####################

count_richOut <- cbind(
  table(count_df$type, count_df$is_outlier),
  round(prop.table(table(count_df$type, count_df$is_outlier), margin = 1), 2)
)
colnames(count_richOut) <- c("in_range", "out_range", "in_range_pct", "out_range_pct")
count_richOut <- data.frame(type = rownames(count_richOut), count_richOut, row.names = NULL)

count_fullSingle <- cbind(
  table(count_df$type, count_df$singleton_frac == 1),
  round(prop.table(table(count_df$type, count_df$singleton_frac == 1), margin = 1), 2)
)
colnames(count_fullSingle) <- c("not_fullSingleton", "fullSingleton", "not_fullSingleton_pct", "fullSingleton_pct")
count_fullSingle <- data.frame(type = rownames(count_fullSingle), count_fullSingle, row.names = NULL)

count_both <- cbind(
  table(count_df$type, count_df$is_outlier | count_df$singleton_frac == 1),
  round(prop.table(table(count_df$type, count_df$is_outlier | count_df$singleton_frac == 1), margin = 1), 2)
)
colnames(count_both) <- c("pass", "fail", "pass_pct", "fail_pct")
count_both <- data.frame(type = rownames(count_both), count_both, row.names = NULL)

filter_out <- data.frame(
  Type = count_both$type,
  Pass = count_both$pass,
  Pass_pct = count_both$pass_pct,
  out_range = count_richOut$out_range,
  full_singleton = count_fullSingle$fullSingleton,
  Fail_pct = count_both$fail_pct
)

filter_out

write.table(filter_out, file.path(input_dir,"passing_filter_stats.tsv"), quote = F,sep = "\t",row.names = F)

filtered_sid <- count_df$original_sid[count_df$is_outlier==TRUE | count_df$singleton_frac == 1]
filtered_sid <- unique(c(filtered_sid, outlier_ids))

keep_sid <- count_df$original_sid[!count_df$original_sid %in% filtered_sid]
writeLines(filtered_sid,file.path(input_dir,"sid_failed_filter.txt"))
writeLines(keep_sid,file.path(input_dir,"sid_passing_filter.txt"))

table(count_df$type[count_df$SAMPLE %in% keep_sid])


# filter from all
count_df <- count_df %>% filter(original_sid %in% keep_sid)
pop.df <- pop.df %>% filter(original_sid %in% keep_sid)
ru_dist <- ru_dist[rownames(ru_dist)%in%keep_sid, colnames(ru_dist)%in%keep_sid]
LGTmat <- LGTmat[rownames(LGTmat) %in% keep_sid,]
LGTmat <- LGTmat[,colSums(LGTmat)>0]


pairs_df <- pairs_df %>% filter(seqID_1 %in% keep_sid) %>% filter(seqID_2 %in% keep_sid)

length(keep_sid)
n_distinct(count_df$original_sid)
n_distinct(pop.df$original_sid)
dim(ru_dist)
nrow(LGTmat)
n_distinct(c(pairs_df$seqID_1,pairs_df$seqID_2))
keep_sid[!keep_sid %in% c(pairs_df$seqID_1,pairs_df$seqID_2)]


###############
# prepare data
###############
library(vegan)
#my_matrix <- LGTmat[,colSums(LGTmat)>2]
#my_matrix <- my_matrix[rowSums(my_matrix)>0,]
#my_matrix <- my_matrix[sort(rowSums(my_matrix)),]
# select_cols <- rowSums(convert_to_presence(my_matrix))
# my_matrix <- my_matrix[select_cols>6,]
# my_matrix <- my_matrix[select_cols>6,]

convert_to_presence <- function(species_matrix){
  x <- species_matrix
  x[x>0] <- 1
  return(x)
}

my_matrix <- LGTmat[sort(rownames(LGTmat)), ]
sample_types <- c("DHF", "DDF", "DHV", "DDV", "DHP", "DDH")

kept_cols_list <- lapply(sample_types, function(x) {
  get_samples <- substr(rownames(my_matrix), 1, 3)
  mat <- my_matrix[get_samples == x, , drop = FALSE]
  col_counts <- colSums(mat > 0)
  cols_to_keep <- names(col_counts)[col_counts > 3]
  return(cols_to_keep)
})

unique_keep_cols <- unique(unlist(kept_cols_list))
my_matrix <- my_matrix[, unique_keep_cols, drop = FALSE]
my_matrix <- my_matrix[rowSums(convert_to_presence(my_matrix))>0,]


dim(my_matrix)

# jaccard_dist <- vegdist(convert_to_presence(my_matrix), method = "jaccard")
ruzicka_dist <- vegdist(my_matrix, method = "jaccard")

sLGT_clean <- 1-as.matrix(ruzicka_dist)
# sLGT_clean <- 1-as.matrix(jaccard_dist)
sLGT_clean <- as.dist(sLGT_clean)

uniqueLGT_list <- keep_sid

# merge
sLGT_long <- sLGT_clean %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(seqID_1 = rownames(.)) %>%
  pivot_longer(
    cols = -seqID_1,
    names_to = "seqID_2",
    values_to = "shared_LGT"
  ) %>%
  filter(seqID_2 > seqID_1) %>%
  filter(seqID_1 != seqID_2) %>%
  mutate(seqID1_2 = paste(seqID_1, seqID_2, sep = "_")) %>%
  select(-seqID_1, -seqID_2) %>%
  left_join(pairs_df, by="seqID1_2", keep = FALSE)

table(sLGT_long$relation)
n_distinct(c(sLGT_long$seqID_1,sLGT_long$seqID_2))

## plot
library(dplyr)
library(ggplot2)
library(purrr)
library(patchwork)

cleaned_df_sameInd <-  sLGT_long %>%
  filter(relation %in% c("replicate", "same_person_different_type")) %>%
  mutate(
    host = factor(sub("_.*", "", host_pair), c("human", "dog")),
    sType = paste0(host, type)
  ) %>% 
  filter(seqID_1 != seqID_2) %>% 
  arrange(desc(shared_LGT))

cleaned_df <- sLGT_long %>% 
  filter(relation %in% c("same_household", "same_village")) %>% 
  filter(host_pair %in% c("human_human", "dog_dog")) %>% 
  filter(type != "MIXED") %>% 
  filter(population != "MIXED") %>% 
  mutate(
    host = factor(sub("_.*", "", host_pair), c("human", "dog")),
    sType = paste0(host, type)
  ) %>% 
  filter(seqID_1 != seqID_2) #%>% 
  #filter(sample_1 != sample_2) %>%
  #bind_rows(cleaned_df_sameInd)

# Generate individual flipped plots
raw_plot_list <- cleaned_df %>%
  group_split(type, host) %>%
  map(function(sub_df) {
    current_type <- as.character(sub_df$type[1])
    current_host <- as.character(sub_df$host[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " - ", current_host, " (n = ", n_unique, ")")
    
    ggplot(sub_df) +
      geom_point(aes(x = relation, y = shared_LGT, colour = sType), size=1.5, alpha=0.5,
                 shape = 16, position = position_jitter(width = 0.2, height = 0)) +
      geom_boxplot(aes(x = relation, y = shared_LGT), outlier.shape = NA, fill="white", alpha=0.5) +
      scale_colour_manual(values = site_col) +
      labs(title = title_text) +
      coord_flip() +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(size = 12, face = "bold"))
  })

# Strip x-axis text from all plots except the final bottom plot
n_plots <- length(raw_plot_list)
final_plot_list <- map2(raw_plot_list, seq_along(raw_plot_list), function(p, idx) {
  if (idx < n_plots) {
    p <- p + theme(
      #axis.text.x  = element_blank(),
      #axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  }
  return(p)
})

# Collate vertically using patchwork
combined_plot <- wrap_plots(final_plot_list, ncol = 1) + 
  plot_layout(guides = "collect")

combined_plot


## bar
raw_plot_list <- cleaned_df %>%
  group_split(type, host) %>%
  map(function(sub_df) {
    current_type <- as.character(sub_df$type[1])
    current_host <- as.character(sub_df$host[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " - ", current_host, " (n = ", n_unique, ")")
    
    nPairs <- table(sub_df$relation)
    
    ggplot(sub_df, aes(x = relation, y = shared_LGT, fill = sType)) +
      stat_summary(fun = "mean", geom = "bar", alpha = 0.8, width = 0.6) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      scale_fill_manual(values = site_col) +
      scale_x_discrete(labels = function(x) {
        paste0(x, "\n(nPairs = ", nPairs[x], ")")
      }) +
      labs(title = title_text) +
      coord_flip() +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(size = 12, face = "bold"))
  })

# Strip axis elements from all plots except the final bottom plot
n_plots <- length(raw_plot_list)
final_plot_list <- map2(raw_plot_list, seq_along(raw_plot_list), function(p, idx) {
  if (idx < n_plots) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  } else {
    p <- p + theme(
      axis.title.y = element_blank()
    )
  }
  return(p)
})

# Collate vertically using patchwork
combined_plot <- wrap_plots(final_plot_list, ncol = 1) + 
  plot_layout(guides = "collect") &
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(5, "mm"),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

nSandwich_retained <- ncol(my_matrix)
nSample_retained <- nrow(my_matrix)
nPairs_retained <- n_distinct(cleaned_df$seqID1_2)

filter_message <- paste0("Retained:\n(1) A set of unique LGT appearing >3 samples per sample type (",nSandwich_retained,");",
"\n(2) Samples which has at least one of the LGT set above (",nSample_retained,")",
"\n(3) Same-type samples pairs from a set that fulfill both of the above (",nPairs_retained,")")

combined_plot <- combined_plot +
  plot_annotation(
    title = "Sample Simiarity by LGT Profile",
    subtitle = "(1-D of Ruzicka's distance metric)",
    caption =  filter_message )

combined_plot

save_plot(combined_plot, "houshold_village_compare_ruzicka_dist_compare(1-D)", height = 20, width = 15, unit = "cm", scale=1.5)














###############
library(ggridges)

cleaned_df %>% 
  #filter(host_pair != "human_dog" & type != "MIXED" & population != "MIXED") %>% 
  mutate(
    population = factor(population, pop.ord),
    host  = factor(sub("_.*", "", host_pair), c("human", "dog")),
    sType = factor(paste0(host, type), type.ord)
  ) %>%
  filter(population != "LDY") %>%
  filter(seqID_1 != seqID_2) %>%
  filter(sample_1 != sample_2) %>%
  ggplot(aes(x = shared_LGT, y = relation, colour = sType, fill = sType)) +
  geom_density_ridges(alpha = 0.3) +
  facet_grid(population ~ sType, scales = "free") +
  theme_minimal(base_size = 20) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA),
    panel.grid = element_blank(),
    legend.position    = "bottom",
    axis.text.x        = element_text(size = 10)
  ) +
  scale_y_discrete(limits = rev) +
  #scale_x_continuous(breaks = seq(0, 100, by = 0.2)) +
  scale_fill_manual(values = site_col) +
  scale_colour_manual(values = site_col) +
  labs(
    y = "Community",
    x = "LGT Similarity"
  )
