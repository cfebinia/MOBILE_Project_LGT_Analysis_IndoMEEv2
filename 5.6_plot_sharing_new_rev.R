rm(list=ls())
gc()
set.seed(7)

library(dplyr)
library(tidyr)
library(stringr)
library(lmerTest)


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
LGT_count="nLGT_by_seqID.tsv" #"nLGT_by_sample.tsv" #
dist_df="dist/ruzicka_uniqueLGT_multisite_bySampleID.rds"
LGTmat="dist/ruzicka_uniqueLGT_multisite_bySampleID_LGTmatrix.rds"
pairs_count="Unique_LGT_pairs_multisite_bySampleID.tsv"
outlier="outliers_waafle.txt"
sPairs="all_seqID_pairs_waafle.tsv"

site_col <- c(
  "humanFecal"="#56B4E9", 
  "humanSkin"="#EE2C2C", 
  "humanOral"="#ef9f00", 
  "dogFecal"="#2850c8", 
  "dogSkin"="#780247", 
  "dogOral"="#9A5324"  
)

pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")
type.ord <- c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")

pLGT <- "pLGT_bySample_withRep.tsv"
nContigs <- "nContig_bySample_withRep.tsv"


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

count_df <- read.delim(file.path(summary_dir,LGT_count)) %>%
  filter(!original_sid %in% outlier_ids) %>%
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

count_df %>%
  group_by(type) %>%
  summarise(nSample=n_distinct(SAMPLE))

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

list_lgt_by_type <- sapply(sample_types, function(x){
  sub_mat <- my_matrix[grepl(x,rownames(my_matrix)),]
  sub_mat <- sub_mat[,which(colSums(sub_mat) > 0)]
  return(colnames(sub_mat))
})
list_lgt_by_type

list_sample_by_type <- sapply(sample_types, function(x){
  sub_mat <- my_matrix[grepl(x,rownames(my_matrix)),]
  return(rownames(sub_mat))
})
list_sample_by_type

# jaccard_dist <- vegdist(convert_to_presence(my_matrix), method = "jaccard")
ruzicka_dist <- vegdist(my_matrix, method = "jaccard")

sLGT_clean <- 1-as.matrix(ruzicka_dist)
# sLGT_clean <- 1-as.matrix(jaccard_dist)
sLGT_clean <- as.dist(sLGT_clean)

uniqueLGT_list <- keep_sid

# pair count
pa_mat <- (my_matrix > 0) * 1

shared_features <- pa_mat %*% t(pa_mat)
sample_pairs <- which(lower.tri(shared_features), arr.ind = TRUE)
shared_features_long <- data.frame(
  seqID_1 = pmin(rownames(shared_features)[sample_pairs[, 1]], colnames(shared_features)[sample_pairs[, 2]]),
  seqID_2 = pmax(rownames(shared_features)[sample_pairs[, 1]], colnames(shared_features)[sample_pairs[, 2]]),
  Shared_Count = shared_features[sample_pairs]
) %>%
  filter(seqID_1 != seqID_2) %>%
  mutate(seqID1_2 = paste(seqID_1, seqID_2, sep = "_"))

sample_richness <- rowSums(pa_mat)
unique_counts <- outer(sample_richness, sample_richness, FUN = "+") - (2 * shared_features)
sample_pairs <- which(lower.tri(unique_counts), arr.ind = TRUE)
unique_counts_long <- data.frame(
  seqID_1 = pmin(rownames(unique_counts)[sample_pairs[, 1]], colnames(unique_counts)[sample_pairs[, 2]]),
  seqID_2 = pmax(rownames(unique_counts)[sample_pairs[, 1]], colnames(unique_counts)[sample_pairs[, 2]]),
  Unique_Count = unique_counts[sample_pairs]
) %>%
  filter(seqID_1 != seqID_2) %>%
  mutate(seqID1_2 = paste(seqID_1, seqID_2, sep = "_"))

# merge
sLGT_long <- sLGT_clean %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(seqID_1 = rownames(.)) %>%
  pivot_longer(
    cols = -seqID_1,
    names_to = "seqID_2",
    values_to = "LGT_similarity"
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
library(coin)

cleaned_df <- sLGT_long %>% 
  filter(relation %in% c("same_household", "same_village")) %>% 
  filter(host_pair %in% c("human_human", "dog_dog")) %>% 
  filter(type != "MIXED") %>% 
  filter(population != "MIXED") %>% 
  mutate(
    host = factor(sub("_.*", "", host_pair), c("human", "dog")),
    sType = paste0(host, type)
  ) %>% 
  filter(seqID_1 != seqID_2) %>% 
  filter(sample_1 != sample_2)

# Generate individual flipped plots
raw_plot_list <- cleaned_df %>%
  group_split(type, host) %>%
  map(function(sub_df) {
    current_type <- as.character(sub_df$type[1])
    current_host <- as.character(sub_df$host[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " - ", current_host, " (n = ", n_unique, ")")
    
    ggplot(sub_df) +
      geom_point(aes(x = relation, y = LGT_similarity, colour = sType), size=1.5, alpha=0.5,
                 shape = 16, position = position_jitter(width = 0.2, height = 0)) +
      geom_boxplot(aes(x = relation, y = LGT_similarity), outlier.shape = NA, fill="white", alpha=0.5) +
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
    
    print(unique(sub_df$sType))
    sub_df$relation <- factor(sub_df$relation, levels = c("same_household", "same_village"))
    nPairs <- table(sub_df$relation)
    
    r_min <- min(table(sub_df$relation))
    r_max_valid <- max(table(sub_df$relation))
    optimal_multiplier <- round((r_min + r_max_valid) / 2)
    
    if (any(nPairs == 0) || length(unique(sub_df$LGT_similarity[!is.na(sub_df$LGT_similarity)])) < 2) {
      p_stars <- "ns"
    } else {
      # Standard Wilcoxon rank-sum test (Mann-Whitney U)
      test_res <- wilcox.test(
        LGT_similarity ~ relation,
        data = sub_df,
        alternative = "greater",
        exact=FALSE
      )
      
      print(test_res)
      
      p_val <- test_res$p.value
      p_stars <- ifelse(p_val < 0.001, "***",
                        ifelse(p_val < 0.01,  "**",
                               ifelse(p_val < 0.05,  "*", "ns")))
    }
    
    max_y <- max(mean_se(sub_df$LGT_similarity[sub_df$relation=="same_household"]), na.rm = TRUE)
    if (is.infinite(max_y)) max_y <- 0
    bracket_y <- max_y * 1.02 
    tip_y <- max_y * 1.01
    
    ggplot(sub_df, aes(x = relation, y = LGT_similarity, fill = sType)) +
      stat_summary(fun = "mean", geom = "bar", alpha = 0.8, width = 0.6) +
      annotate("path", x = c(1, 1, 2, 2), y = c(tip_y, bracket_y, bracket_y, tip_y), colour = "black", linewidth = 0.5) +
      annotate("text", x = 1.5, y = bracket_y * 1.05, label = p_stars, size = 5, vjust = 0.3, angle = 0) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      scale_fill_manual(values = site_col) +
      scale_x_discrete(expand = expansion(mult = c(0.05, 0.1)),
        labels = function(x) {
        paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
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

filter_message <- paste0("P-value from one-tailed Wilcoxon Rank-Sum test, whereby:
                         ns = not significant; * = p<0.05; ** = p<0.01; *** = p<0.001.
                         Dataset:\n(1) A set of unique LGT appearing >3 samples per sample type (",nSandwich_retained,");",
"\n(2) Samples which has at least one of the LGT set above (",nSample_retained,")",
"\n(3) Same-type samples pairs from a set that fulfill both of the above (",nPairs_retained,")")

combined_plot <- combined_plot +
  plot_annotation(
    title = "Sample Simiarity by LGT Profile",
    subtitle = "(1-D of Ruzicka's distance metric)",
    caption =  filter_message )

combined_plot

save_plot(combined_plot, "household_village_compare_ruzicka_dist_compare(1-D)", height = 20, width = 15, unit = "cm", scale=1.5, dpi=300)



### Add different-village
richness_outlier_summary <- count_df %>%
  group_by(type) %>%
  mutate(
    q_low  = quantile(count_unique, 0.025),
    q_high = quantile(count_unique, 0.975),
    is_low  = count_unique < q_low,
    is_high = count_unique > q_high
  ) %>%
  summarise(
    total_samples = n(),
    low_outliers  = sum(is_low),
    high_outliers = sum(is_high),
    total_excluded = sum(is_low | is_high),
    excluded_pct  = round(total_excluded / total_samples * 100, 2)
  )

pre_exclusion <- round(100 * (
  sum(richness_outlier_summary$total_excluded) /
    sum(richness_outlier_summary$total_samples)
),0)

pre_exclusion <- paste(", ",pre_exclusion,"%", sep="")
n_exclusion <- paste(sum(richness_outlier_summary$total_excluded), 
                     " of ",  
                     sum(richness_outlier_summary$total_samples), sep="")

n_exclusion


base_df <- cleaned_df

get_seqID <- union(base_df$seqID_1, base_df$sample_2)
get_sampleID <- union(base_df$sample_1, base_df$sample_2)
exclude_pairs_seqID <- c(base_df$seqID1_2)
exclude_pairs_sampleID <- c(base_df$sample1_2)

different_village <- sLGT_long %>%
  filter(population=="MIXED" & type != "MIXED" & host_pair != "human_dog") %>%
  mutate(host=ifelse(host_pair=="dog_dog", "dog", "human"),
         sType=paste(host,type, sep=""),
         relation="different_village")

table(different_village$host_pair)
table(different_village$population)
table(different_village$type)
table(different_village$relation)

different_village %>%
  group_by(sType) %>%
  summarise(nSample=n_distinct(union(sample_1, sample_2)))

dim(different_village)
dim(base_df)

plot_df <- bind_rows(base_df, different_village) %>%
  mutate(host=factor(host, c("human","dog")),
         type=factor(type, c("Fecal","Oral","Skin")),
         relation=factor(relation, c("different_village","same_village","same_household"))) %>%
  arrange(type, host) %>%
  mutate(sType=factor(sType, unique(sType)))

write.table(plot_df, file.path(basedir, "summary_tables/LGT_similarity_by_social_relation.tsv"), quote=F, row.names=F, sep="\t")

relation_levels <- c("different_village", "same_village", "same_household")

raw_plot_list <- plot_df %>%
  group_split(type, host) %>%
  map(function(sub_df) {
    current_type <- as.character(sub_df$type[1])
    current_host <- as.character(sub_df$host[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " - ", current_host, " (n = ", n_unique, ")")
    
    # Ensure factor ordering (1: different_village, 2: same_village, 3: same_household)
    sub_df$relation <- factor(sub_df$relation, levels = relation_levels)
    nPairs <- table(sub_df$relation)
    
    # ----------------------------------------------------
    # 1. PAIRWISE WILCOXON TESTS
    # ----------------------------------------------------
    # Filter non-empty groups with sufficient data
    valid_groups <- sub_df %>%
      group_by(relation) %>%
      summarise(n = sum(!is.na(LGT_similarity)), .groups = "drop") %>%
      filter(n > 0) %>%
      pull(relation)
    
    df_brackets <- tibble()
    
    if (length(valid_groups) >= 2 && length(unique(na.omit(sub_df$LGT_similarity))) >= 2) {
      
      # Pairwise Wilcoxon rank-sum test with p-adjustment
      pw_res <- pairwise.wilcox.test(
        x = sub_df$LGT_similarity,
        g = sub_df$relation,
        p.adjust.method = "fdr",
        exact = FALSE
      )
      
      # Convert p-value matrix into a tidy frame for plotting
      p_mat <- pw_res$p.value
      
      comp_list <- list()
      for (r in rownames(p_mat)) {
        for (c in colnames(p_mat)) {
          if (!is.na(p_mat[r, c])) {
            comp_list[[length(comp_list) + 1]] <- tibble(
              group1 = c,
              group2 = r,
              p_val  = p_mat[r, c]
            )
          }
        }
      }
      
      if (length(comp_list) > 0) {
        df_brackets <- bind_rows(comp_list) %>%
          mutate(
            x1 = as.numeric(factor(group1, levels = relation_levels)),
            x2 = as.numeric(factor(group2, levels = relation_levels)),
            p_stars = case_when(
              p_val < 0.001 ~ "***",
              p_val < 0.01  ~ "**",
              p_val < 0.05  ~ "*",
              TRUE          ~ "ns"
            )
          )
      }
    }
    
    # ----------------------------------------------------
    # 2. DYNAMIC BRACKET Y-POSITIONING
    # ----------------------------------------------------
    # Find highest bar mean + SE to position brackets safely above data
    max_bar_height <- sub_df %>%
      group_by(relation) %>%
      summarise(
        mean_val = mean(LGT_similarity, na.rm = TRUE),
        se_val   = sd(LGT_similarity, na.rm = TRUE) / sqrt(n()),
        top      = mean_val + se_val,
        .groups  = "drop"
      ) %>%
      pull(top) %>%
      max(na.rm = TRUE)
    
    if (is.infinite(max_bar_height) || is.na(max_bar_height)) max_bar_height <- 0.1
    
    # Calculate staggered bracket positions
    step_height <- max_bar_height * 0.10
    base_y      <- max_bar_height * 1.03
    
    if (nrow(df_brackets) > 0) {
      df_brackets <- df_brackets %>%
        mutate(
          bracket_y = base_y + (row_number() - 1) * step_height,
          tip_y     = bracket_y - (step_height * 0.15),
          text_y    = bracket_y + (step_height * 0.3),
          x_mid     = (x1 + x2) / 2
        )
    }
    
    # ----------------------------------------------------
    # 3. GGPLOT CONSTRUCTION
    # ----------------------------------------------------
    p <- ggplot(sub_df, aes(x = relation, y = LGT_similarity, fill = sType)) +
      stat_summary(fun = "mean", geom = "bar", width = 0.9) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      scale_fill_manual(values = site_col) +
      scale_x_discrete(
        expand = expansion(mult = c(0.1, 0.1)),
        labels = function(x) {
          paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
        }
      ) +
      labs(title = title_text) +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40")
      )
    
    # Add dynamic brackets if comparisons exist
    if (nrow(df_brackets) > 0) {
      for (i in seq_len(nrow(df_brackets))) {
        b <- df_brackets[i, ]
        
        # Draw path bracket (x1, tip_y) -> (x1, bracket_y) -> (x2, bracket_y) -> (x2, tip_y)
        p <- p + 
          annotate(
            "path",
            x = c(b$x1, b$x1, b$x2, b$x2),
            y = c(b$tip_y, b$bracket_y, b$bracket_y, b$tip_y),
            colour = "black", linewidth = 0.5
          ) +
          annotate(
            "text",
            x = b$x_mid,
            y = b$text_y,
            label = b$p_stars,
            size = 4.5,
            vjust = 0.5
          )
      }
    }
    
    return(p)
  })

# Strip x-axis titles and text/labels from all plots except the last one
n_plots <- length(raw_plot_list)

final_plot_list <- map2(raw_plot_list, seq_along(raw_plot_list), function(p, idx) {
  if (idx < n_plots) {
    # Keep tick labels and marks, but remove the axis title for top/middle plots
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 9.5, colour = "grey40")
    )
  } else {
    # Retain x-axis title on the bottom plot only
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 12, margin = margin(t = 10))
    )
  }
  return(p)
})

# Calculate summary metrics from data
nSample_retained <- length(unique(c(plot_df$seqID_1, plot_df$seqID_2)))
nPairs_retained  <- dplyr::n_distinct(plot_df$seqID1_2)

# Build multi-line caption message
filter_message <- paste0(
  "P-values from pairwise Wilcoxon Rank-Sum tests with FDR adjustment:\n",
  "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001.\n\n",
  "Dataset filtering:\n",
  "(1) Initial set of samples falling within of 95% quantile range in LGT Richness per sample type (n = ", length(keep_sid), ");\n",
  "(2) Unique LGT set appearing in >3 samples per sample type (n = ", nSandwich_retained, ");\n",
  "(3) Samples containing at least one LGT from the above set (n = ", nSample_retained, ");\n",
  "(4) Same-type pairwise combinations, e.g. Fecal-Fecal (n = ", nPairs_retained, ")."
)

# Assemble composite layout with annotations
combined_plot <- wrap_plots(final_plot_list, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Sample Similarity by LGT Profile",
    subtitle = "(1 - Ruzicka Distance Metric)",
    caption  = filter_message
  ) &
  theme(
    # Titles & Captions
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.title    = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0, margin = margin(b = 10)),
    plot.caption  = element_text(size = 9, hjust = 0, color = "grey30", margin = margin(t = 12)),
    
    # Legend formatting
    legend.title    = element_blank(),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(5, "mm"),
    legend.position = "bottom",
    legend.box      = "horizontal"
  )

# Display plot
combined_plot

save_plot(combined_plot, "household_village_compare_ruzicka_dist_compare(1-D)_v2", height = 21, width = 15, unit = "cm", scale=1.5, dpi=300)

### Final report of sample count and filtering categories
sink(file = "summary_tables/valid_social_dataset.txt", split = TRUE)
options(echo = TRUE)

all <- read.delim(file.path(summary_dir, LGT_count))
atleast_1_lgt <- length(unique(all$original_sid))
print("atleast_1_lgt:")
print(atleast_1_lgt)

print("length(outlier_ids):")
print(length(outlier_ids))

print("count_fullSingle:")
print(count_fullSingle)

print("sum(count_fullSingle$fullSingleton):")
print(sum(count_fullSingle$fullSingleton))

print("richness_outlier_summary:")
print(richness_outlier_summary)

print("sum(richness_outlier_summary$total_excluded):")
print(sum(richness_outlier_summary$total_excluded))

n_valid_lgt <- ncol(my_matrix)
print("n_valid_lgt:")
print(n_valid_lgt)

n_valid_lgt_by_Type <- sapply(list_lgt_by_type, length)
print("n_valid_lgt_by_Type:")
print(n_valid_lgt_by_Type)

n_valid_samples <- nrow(my_matrix)
print("n_valid_samples:")
print(n_valid_samples)

n_valid_samples_byType <- sapply(list_sample_by_type, length)
print("n_valid_samples_byType:")
print(n_valid_samples_byType)

print("List of valid samples:")
print(union(plot_df$seqID_1, plot_df$seqID_2))

print("List of included individuals:")
print(union(plot_df$sample_1, plot_df$sample_2))
sink()



####################
## pLGT with Dogs ##
####################
library(dplyr)
library(ggplot2)
library(purrr)
library(ggpubr)

hh_dogs <- read.delim(file.path(input_dir, "Household_counts_filledv2.tsv")) %>%
  filter(!is.na(sampleid))

sid_with_dogs <- hh_dogs$sampleid[which(hh_dogs$Dogs_HH > 0)]

nLGT_withDogs <- read.delim("summary_tables/nLGT_by_seqID.tsv") %>%
  filter(host == "human") %>%
  mutate(
    with_dogs = ifelse(individual %in% sid_with_dogs, "with_dogs", "no_dogs"),
    type      = paste0(host, site),
    type      = factor(type, levels = c("humanFecal", "humanOral", "humanSkin"))
  ) %>%
  filter(group != "Agriculture") %>%
  mutate(pLGT=pLGT*100)

wDogs_fit <- list()

for (i in unique(nLGT_withDogs$type)) {
  sub_df <- nLGT_withDogs %>%
    filter(type == i) %>%
    mutate(with_dogs = ifelse(with_dogs == "with_dogs", TRUE, FALSE)) %>%
    rename(community = population)
  
  fit <- lmer(pLGT ~ with_dogs * group + (1 | community), data = sub_df)
  
  wDogs_fit[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}

raw_plot_list <- nLGT_withDogs %>%
  group_split(type, .keep = TRUE) %>%
  keep(~ nrow(.x) > 0) %>% 
  map(function(sub_df) {
    
    title_text <- as.character(sub_df$type[1])
    
    counts <- sub_df %>% 
      count(with_dogs) %>% 
      tibble::deframe()
    
    my_comparisons <- list(c("no_dogs", "with_dogs"))
    
    # 1. Compute top of the bar + error bar for positioning
    # (Mean + Standard Error per group)
    bar_tops <- sub_df %>%
      group_by(with_dogs) %>%
      summarise(
        mean_val = mean(pLGT, na.rm = TRUE),
        se_val   = sd(pLGT, na.rm = TRUE) / sqrt(n()),
        top      = mean_val + se_val
      )
    
    max_bar_height <- max(bar_tops$top, na.rm = TRUE)
    
    # Position bracket just 5% above highest error bar
    bracket_pos <- max_bar_height * 1.05 
    
    nSample <- sub_df %>%
      group_by(with_dogs) %>%
      summarise(nSample = n_distinct(individual))
    x_labels <- setNames(paste0(nSample$with_dogs, "\n (", nSample$nSample, ")"), nSample$with_dogs)
    
    p <- ggplot(sub_df, aes(x = with_dogs, y = pLGT, fill = type)) +
      stat_summary(fun = "mean", geom = "bar", width = 0.6) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      
      # --- WILCOXON BRACKET ---
      stat_compare_means(
        method = "wilcox.test",
        comparisons = my_comparisons,
        label = "p.signif",           
        hide.ns = FALSE,             
        label.y = bracket_pos         # Placed snugly above error bar
      ) +
      
      # Expand y-limit slightly so bracket label doesn't get clipped
      scale_y_continuous(limits = c(0, max_bar_height * 1.2)) +
      
      scale_x_discrete(
        #expand = expansion(mult = c(0.2, 0.2)),
        labels = x_labels) +
      labs(title = title_text, x = "Dog Ownership", y = "pLGT (%)") +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40"),
        panel.border = element_rect(colour="grey40", linewidth = 0.5)
      ) +
      scale_fill_manual(values = site_col)
    
    return(p)
  })

names(raw_plot_list) <- c("humanFecal", "humanOral", "humanSkin")

combined_plot <- wrap_plots(raw_plot_list, ncol = 3)

combined_plot
save_plot(combined_plot, "pLGT_with_dogs", width = 15, height = 8, dpi=300,unit = "cm", scale=1.2)

wDogs_coeff_list <- lapply(names(wDogs_fit), function(model_name) {
  model <- wDogs_fit[[model_name]]
  
  df <- data.frame(type=model_name,
                   as.data.frame(coef(summary(model))))
  df$term <- rownames(df)
  rownames(df) <- NULL
  
  p_col <- grep("Pr", colnames(df), value = TRUE)[1]
  
  if (!is.null(p_col) && !is.na(p_col)) {
    df$p.stars <- symnum(
      df[[p_col]],
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  } else {
    df$p.stars <- " "
  }
  
  r2_val <- MuMIn::r.squaredGLMM(model)
  df$marginal_r2 <- round(r2_val[1, "R2m"], 3)
  df$conditional_r2 <- round(r2_val[1, "R2c"], 3)
  
  return(df)
})

wDogs_coeff_df <- do.call(rbind, wDogs_coeff_list)
wDogs_coeff_df

out_path <- file.path(basedir, "summary_tables/lmer_with_dogs_pLGT.tsv")

writeLines("# LMER model : pLGT ~ with_dogs[1/0] * group + (1 | community)", con = out_path)

write.table(
  wDogs_coeff_df,
  file = out_path,
  append = T,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)




#############################
## By Household Membership ##
#############################
library(lmerTest)

metadata <- read.delim(file.path(input_dir,"combined_ind_metadata.tsv")) %>%
  rename(sampleid=Sample) %>%
  select(sampleid, age, sex) %>%
  distinct()

hh_id <- read.delim("input_file/socials/addresses.txt")

get_mode <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA)
  ux <- unique(x_clean)
  ux[which.max(tabulate(match(x_clean, ux)))]
}

hh_member <- read.delim(file.path(input_dir, "Household_counts_filledv2.tsv")) %>%
  filter(!is.na(Total_HH)) %>%
  left_join(hh_id, by = "sampleid") %>%
  mutate(Rel_young = Rel_children + Rel_childrenstep + Rel_neicenephewinlaw + Rel_grandchildren) %>%
  select(sampleid, population, HH_ID, Total_HH, Rel_young) %>%
  left_join(metadata, by = "sampleid") %>%
  group_by(population, HH_ID) %>%
  filter(Total_HH == get_mode(Total_HH)) %>%
  ungroup() %>%
  rename(individual=sampleid) %>%
  select(-population)


nLGT_HH <- read.delim("summary_tables/nLGT_by_seqID.tsv") %>%
  filter(host == "human") %>%
  inner_join(hh_member, by = "individual") %>%
  mutate(sType = paste0(host, site)) %>%
  rename(community = population) %>%
  mutate(
    with_young = ifelse(Rel_young > 0, "with_young", "no_young"),
    type      = paste0(host, site),
    type      = factor(type, levels = c("humanFecal", "humanOral", "humanSkin"))
  ) %>%
  filter(group != "Agriculture") %>%
  mutate(pLGT=pLGT*100)

totalHH_fit <- list()

for (i in unique(nLGT_HH$sType)) {
  sub_df <- subset(nLGT_HH, sType == i)
  
  fit <- lmer(pLGT ~ Total_HH + group + (1 | community), data = sub_df)
  
  totalHH_fit[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}


relYoung_fit <- list()

for (i in unique(nLGT_HH$sType)) {
  sub_df <- subset(nLGT_HH, sType == i & age < 50)
  
  fit <- lmer(pLGT ~ Rel_young * age + group + (1 | community), data = sub_df)
  
  relYoung_fit[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}


relYoung_fit2 <- list()

for (i in unique(nLGT_HH$sType)) {
  sub_df <- subset(nLGT_HH, sType == i & age < 50)
  sub_df$Rel_young <- ifelse(sub_df$Rel_young > 0, TRUE, FALSE)
  
  fit <- lmer(pLGT ~ Rel_young * age + group + (1 | community), data = sub_df)
  
  relYoung_fit2[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}

raw_plot_list <- nLGT_HH %>%
  filter(age < 50) %>%
  group_split(type, group, .keep = TRUE) %>%
  keep(~ nrow(.x) > 0) %>% 
  map(function(sub_df) {
    
    current_group <- as.character(sub_df$group)[1]
    title_text <- as.character(sub_df$type[1])
    counts <- sub_df %>% 
      count(with_young) %>% 
      tibble::deframe()
    
    my_comparisons <- list(c("no_young", "with_young"))
    
    # 1. Compute top of the bar + error bar for positioning
    # (Mean + Standard Error per group)
    bar_tops <- sub_df %>%
      group_by(with_young) %>%
      summarise(
        mean_val = mean(pLGT, na.rm = TRUE),
        se_val   = sd(pLGT, na.rm = TRUE) / sqrt(n()),
        top      = mean_val + se_val
      )
    
    max_bar_height <- max(bar_tops$top, na.rm = TRUE)
    
    # Position bracket just 5% above highest error bar
    bracket_pos <- max_bar_height * 1.05 
    
    nSample <- sub_df %>%
      group_by(with_young) %>%
      summarise(nSample = n_distinct(individual))
    x_labels <- setNames(paste0(nSample$with_young, "\n (", nSample$nSample, ")"), nSample$with_young)
    
    p <- ggplot(sub_df, aes(x = with_young, y = pLGT, fill = type)) +
      facet_grid(type~group, scale="free")+
      stat_summary(fun = "mean", geom = "bar", width = 0.6) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      
      # --- WILCOXON BRACKET ---
      stat_compare_means(
        method = "wilcox.test",
        comparisons = my_comparisons,
        label = "p.signif",           
        hide.ns = FALSE,             
        label.y = bracket_pos         # Placed snugly above error bar
      ) +
      
      # Expand y-limit slightly so bracket label doesn't get clipped
      scale_y_continuous(limits = c(0, max_bar_height * 1.2)) +
      
      scale_x_discrete(
        #expand = expansion(mult = c(0.2, 0.2)),
        labels = x_labels
      ) +
      labs(title = title_text, y = "pLGT (%)") +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40"),
        panel.border = element_rect(colour="grey40", linewidth = 0.5)
      ) +
      scale_fill_manual(values = site_col)
    
    if(current_group == "Late-transition"){
      p <- p + theme(plot.title = element_blank())
    }
    return(p)
  })

names(raw_plot_list) <- c("humanFecal", "humanOral", "humanSkin")

combined_plot <- wrap_plots(raw_plot_list, ncol = 2)

combined_plot

save_plot(combined_plot, "pLGT_with_young", width = 15, height = 15, dpi=300,unit = "cm", scale=1.2)

relYoung_coeff_list <- lapply(names(relYoung_fit), function(model_name) {
  model <- relYoung_fit[[model_name]]
  
  df <- data.frame(type=model_name,
                   as.data.frame(coef(summary(model))))
  df$term <- rownames(df)
  df$type <- model_name
  rownames(df) <- NULL
  
  p_col <- grep("Pr", colnames(df), value = TRUE)[1]
  
  if (!is.na(p_col)) {
    df$p.stars <- symnum(
      df[[p_col]],
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  } else {
    df$p.stars <- " "
  }
  
  r2_val <- MuMIn::r.squaredGLMM(model)
  df$marginal_r2 <- round(r2_val[1, "R2m"],3)
  df$conditional_r2 <- round(r2_val[1, "R2c"],3)
  
  return(df)
})

relYoung_coeff_df <- do.call(rbind, relYoung_coeff_list)
relYoung_coeff_df


out_path <- file.path(basedir, "summary_tables/lmer_with_Relyoung_pLGT.tsv")

writeLines("# LMER model : pLGT ~ Rel_young * age + group + (1 | community)", con = out_path)

write.table(
  relYoung_coeff_df,
  file = out_path,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE, 
  append = T
)

sapply(relYoung_fit, nobs)

###########################
## sLGT By Dog Ownership ##
###########################
hh_dogs <- read.delim(file.path(input_dir, "Household_counts_filledv2.tsv")) %>%
  filter(!is.na(sampleid))

sid_with_dogs <- hh_dogs$sampleid[(!is.na(hh_dogs$Dogs_HH) & hh_dogs$Dogs_HH > 0)]

seqID_with_dogs <- pop.df %>%
  filter(host=="human") %>%
  filter(individual %in% hh_dogs$sampleid[!is.na(hh_dogs$Dogs_HH)]) %>%
  pull(SAMPLE) %>%
  unique()

sLGT_Dogs <- plot_df %>%
  filter(host == "human") %>%
  filter(population != "LDY") %>%
  filter(seqID_1 %in% seqID_with_dogs, seqID_2 %in% seqID_with_dogs) %>%
  mutate(with_dogs = case_when(
    (sample_1 %in% sid_with_dogs) & (sample_2 %in% sid_with_dogs) ~ "both_has_dogs",
    (sample_1 %in% sid_with_dogs) | (sample_2 %in% sid_with_dogs) ~ "one_has_dogs",
    TRUE ~ "no_dogs"
  )) %>%
  mutate(with_dogs = case_when(
    relation == "same_household" & with_dogs != "no_dogs" ~ "both_has_dogs",
    TRUE ~ with_dogs
  )) %>%
  mutate(
    with_dogs = factor(with_dogs, c("both_has_dogs", "one_has_dogs", "no_dogs")),
    group     = factor(relation, c("same_household", "same_village", "different_village"))
  )

raw_plot_list <- sLGT_Dogs %>%
  group_split(type, host, group) %>%
  map(function(sub_df) {
    current_type  <- as.character(sub_df$type[1])
    current_host  <- as.character(sub_df$host[1])
    current_group <- as.character(sub_df$group[1])
    
    n_unique   <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_group, " (n = ", n_unique, ")")
    
    sub_df$with_dogs <- factor(sub_df$with_dogs, levels = c("both_has_dogs", "one_has_dogs", "no_dogs"))
    nPairs           <- table(sub_df$with_dogs)
    
    if(current_group == "same_household"){
      sub_df <- sub_df %>%
        group_by(group) %>%
        complete(
          with_dogs = factor(
            c("both_has_dogs", "one_has_dogs", "no_dogs"),
            levels = c("both_has_dogs", "one_has_dogs", "no_dogs"))
        ) %>%
        mutate(LGT_similarity = replace_na(LGT_similarity, 0)) %>%
        ungroup()
    }
    
    # ----------------------------------------------------
    # 1. PAIRWISE WILCOXON TESTS
    # ----------------------------------------------------
    valid_groups <- sub_df %>%
      group_by(with_dogs) %>%
      summarise(n = sum(!is.na(LGT_similarity)), .groups = "drop") %>%
      filter(n > 0) %>%
      pull(with_dogs)
    
    if(current_group == "same_household"){
      valid_groups <- valid_groups[!grepl("one_has_dogs", valid_groups)]
    }
    
    df_brackets <- tibble()
    
    if (length(valid_groups) >= 2 && length(unique(na.omit(sub_df$LGT_similarity))) >= 2) {
      
      pw_res <- pairwise.wilcox.test(
        x = sub_df$LGT_similarity,
        g = sub_df$with_dogs,
        p.adjust.method = "none",
        #alternative = "less",
        exact = FALSE
      )
      
      p_mat <- pw_res$p.value
      comp_list <- list()
      
      for (r in rownames(p_mat)) {
        for (c in colnames(p_mat)) {
          if (!is.na(p_mat[r, c])) {
            comp_list[[length(comp_list) + 1]] <- tibble(
              group1    = c,
              group2    = r,
              p_val_raw = p_mat[r, c]
            )
          }
        }
      }
      
      if (length(comp_list) > 0) {
        df_brackets <- bind_rows(comp_list) %>%
          mutate(
            p_val_adj = p.adjust(p_val_raw, method = "fdr"),
            x1        = as.numeric(factor(group1, levels = levels(sub_df$with_dogs))),
            x2        = as.numeric(factor(group2, levels = levels(sub_df$with_dogs))),
            p_stars   = case_when(
              p_val_adj < 0.001 ~ "***",
              p_val_adj < 0.01  ~ "**",
              p_val_adj < 0.05  ~ "*",
              TRUE              ~ "ns"
            )
          )
      }
    }
    
    # ----------------------------------------------------
    # 2. BRACKET Y-POSITIONING
    # ----------------------------------------------------
    max_bar_height <- sub_df %>%
      group_by(with_dogs) %>%
      summarise(
        mean_val = mean(LGT_similarity, na.rm = TRUE),
        se_val   = sd(LGT_similarity, na.rm = TRUE) / sqrt(n()),
        top      = mean_val + se_val,
        .groups  = "drop"
      ) %>%
      pull(top) %>%
      max(na.rm = TRUE)
    
    if (is.infinite(max_bar_height) || is.na(max_bar_height)) max_bar_height <- 0.1
    
    step_height <- max_bar_height * 0.25
    base_y      <- max_bar_height * 1.08
    
    if (nrow(df_brackets) > 0) {
      df_brackets <- df_brackets %>%
        mutate(
          bracket_y = base_y + (row_number() - 1) * step_height,
          tip_y     = bracket_y - (step_height * 0.20),
          text_y    = bracket_y + (step_height * 0.35),
          x_mid     = (x1 + x2) / 2
        )
    }
    
    if(current_group == "same_household"){
      df_brackets <- df_brackets %>%
        filter(group1 != "one_has_dogs", group2 != "one_has_dogs")
    }
    
    # ----------------------------------------------------
    # 3. GGPLOT CONSTRUCTION
    # ----------------------------------------------------
    p <- ggplot(sub_df, aes(x = with_dogs, y = LGT_similarity, fill = sType)) +
      stat_summary(
        fun = "mean", 
        geom = "bar", 
        width = 0.95, 
        position = "identity"
      ) +
      stat_summary(
        fun.data = "mean_se", 
        geom = "errorbar", 
        width = 0.3, 
        linewidth=0.3,
        colour = "black", 
        position = "identity"
      ) +
      scale_fill_manual(values = site_col) +
      scale_x_discrete(
        labels = function(x) {
          paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
        }
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(subtitle = title_text, x = "", y = "LGT Similarity") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40")
      )
    
    if (nrow(df_brackets) > 0) {
      for (i in seq_len(nrow(df_brackets))) {
        b <- df_brackets[i, ]
        
        p <- p + 
          annotate(
            "path",
            x = c(b$x1, b$x1, b$x2, b$x2),
            y = c(b$tip_y, b$bracket_y, b$bracket_y, b$tip_y),
            colour = "black", linewidth = 0.3
          ) +
          annotate(
            "text",
            x = b$x_mid,
            y = b$text_y,
            label = b$p_stars,
            size = 3.5,
            vjust = 0.5
          )
      }
    }
    
    title_text <- paste0(current_host, current_type)
    
    if(current_group=="same_household"){
      p <- p+ labs(title=title_text)
    }
    
    p <- p + theme(legend.position = "none")
    return(p)
  })

# Calculate summary metrics from data
nSample_retained <- length(unique(c(sLGT_Dogs$seqID_1, sLGT_Dogs$seqID_2)))
nPairs_retained  <- dplyr::n_distinct(sLGT_Dogs$seqID1_2)
nSandwich_retained <- list_lgt_by_type[grepl("DH", names(list_lgt_by_type))] %>%
  reduce(union) %>% length()

# Build multi-line caption message
filter_message <- paste0(
  "P-values from pairwise Wilcoxon Rank-Sum tests, not adjusted FDR adjustment:\n",
  "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001.\n\n",
  "Dataset filtering:\n",
  "(1) Initial set of samples falling within of 95% quantile range in LGT Richness per sample type (n = ", length(keep_sid), ");\n",
  "(2) Unique LGT set appearing in >3 samples per sample type (n = ", nSandwich_retained, ");\n",
  "(3) Samples containing at least one LGT from the above set, and ownership data (n = ", nSample_retained, ");\n",
  "(4) Same-type pairwise combinations, e.g. Fecal-Fecal (n = ", nPairs_retained, ")."
)

# Assemble composite layout with annotations
combined_plot <- wrap_plots(raw_plot_list, nrow = 3) +
  plot_layout(
    guides = "collect",
    axis_titles   = "collect"
  ) +
  plot_annotation(
    title    = "Sample Similarity by LGT Profile by Dog Ownership",
    subtitle = "(1 - Ruzicka Distance Metric)",
    caption  = filter_message
  ) &
  theme(
    # Titles & Captions
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.title    = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0, margin = margin(b = 10)),
    plot.caption  = element_text(size = 9, hjust = 0, color = "grey30", margin = margin(t = 12)),
    axis.text = element_text(size = 8),
    
    # Legend formatting
    #legend.title    = element_blank(),
    #legend.text     = element_text(size = 10),
    #legend.key.size = unit(5, "mm"),
    #legend.position = "bottom",
    #legend.box      = "horizontal"
    legend.position = "none"
  )

# Display plot
combined_plot

save_plot(combined_plot, "household_village_compare_ruzicka_dist_compare(1-D)_Ownership", height = 13, width = 15, unit = "cm", scale=1.5, dpi=600)

## LM 
sLGT_withDogs_fit <- list()
nSamples <- list()
for (i in unique(sLGT_Dogs$sType)) {
  sub_df <- subset(sLGT_Dogs, sType == i) %>%
    mutate(with_dogs = factor(with_dogs, levels = c("no_dogs", "one_has_dogs", "both_has_dogs"))) %>%
    rename(community = population) %>%
    left_join(unique_counts_long %>%
                select(seqID1_2, Unique_Count),
              by="seqID1_2") %>%
    mutate(with_dogs=relevel(with_dogs, ref = "no_dogs"))
  
  nSamples[[as.character(i)]] <- length(union(sub_df$seqID_1, sub_df$seqID_2))
  
  fit <- lm(LGT_similarity ~ relation * with_dogs + log10(Unique_Count+1), data = sub_df)
  # fit <- lm(LGT_similarity ~ relation * with_dogs, data = sub_df)
  
  sLGT_withDogs_fit[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}

sLGT_dogsPairs_coef_list <- lapply(names(sLGT_withDogs_fit), function(model_name) {
  model <- sLGT_withDogs_fit[[model_name]]
  mod_sum <- summary(model)
  
  df <- as.data.frame(coef(mod_sum))
  df$term <- rownames(df)
  df$type <- model_name
  rownames(df) <- NULL
  
  p_col <- grep("Pr", colnames(df), value = TRUE)[1]
  
  if (!is.null(p_col) && !is.na(p_col)) {
    df$p.stars <- symnum(
      df[[p_col]],
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  } else {
    df$p.stars <- " "
  }
  
  df$r_squared <- round(mod_sum$r.squared, 3)
  df$adj_r_squared <- round(mod_sum$adj.r.squared, 3)
  df$nPairs <- nobs(model)
  df$nSample <- nSamples[[model_name]]
  
  return(df)
})

sLGT_dogsPairs_coef_df <- do.call(rbind, sLGT_dogsPairs_coef_list)
sLGT_dogsPairs_coef_df

out_path <- file.path(basedir, "summary_tables/lm_sLGT_dogsPairs.tsv")

writeLines("# LM model : LGT_similarity ~ relation * with_dogs + log10(Unique_Count+1)", con = out_path)

write.table(
  sLGT_dogsPairs_coef_df,
  file = out_path,
  append = TRUE,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

sapply(sLGT_withDogs_fit, nobs)

sLGT_Dogs %>%
  group_by(type) %>%
  summarise(nSample = n_distinct(c(sample_1, sample_2)), .groups = "drop")

#######################
## sLGT by Rel young ##
#######################
sid_with_young <- hh_member %>%
  filter(age < 50) %>%
  filter(Rel_young > 0) %>%
  #right_join(pop.df %>% filter(host=="human") %>% select(individual, SAMPLE) %>% distinct(), 
  #           by="individual") %>%
  pull(individual) %>%
  unique()

median(hh_member$Total_HH)
median(hh_member$Rel_young[hh_member$age<50])

seqID_under50 <- hh_member %>%
  filter(age < 50) %>%
  right_join(pop.df %>% filter(host=="human") %>% select(individual, SAMPLE) %>% distinct(), 
             by="individual") %>%
  pull(SAMPLE) %>%
  unique()

sLGT_HH <- plot_df %>%
  filter(host == "human") %>%
  filter(seqID_1 %in%  seqID_under50, seqID_2 %in% seqID_under50) %>%
  mutate(with_young = case_when(
    (sample_1 %in% sid_with_young) & (sample_2 %in% sid_with_young) ~ "both_has_young",
    (sample_1 %in% sid_with_young) | (sample_2 %in% sid_with_young) ~ "one_has_young",
    TRUE ~ "no_young"
  )) %>%
  mutate(with_young = case_when(
    relation == "same_household" & with_young != "no_young" ~ "both_has_young",
    TRUE ~ with_young
  )) %>%
  mutate(
    with_young = factor(with_young, c("both_has_young", "one_has_young", "no_young")),
    group     = factor(relation, c("same_household", "same_village", "different_village"))
  )

raw_plot_list <- sLGT_HH %>%
  group_split(type, host, group) %>%
  map(function(sub_df) {
    current_type  <- as.character(sub_df$type[1])
    current_host  <- as.character(sub_df$host[1])
    current_group <- as.character(sub_df$group[1])
    
    n_unique   <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_group, " (n = ", n_unique, ")")
    
    sub_df$with_young <- factor(sub_df$with_young, levels = c("both_has_young", "one_has_young", "no_young"))
    nPairs           <- table(sub_df$with_young)
    
    if(current_group == "same_household"){
      sub_df <- sub_df %>%
        group_by(group) %>%
        complete(
          with_young = factor(
            c("both_has_young", "one_has_young", "no_young"),
            levels = c("both_has_young", "one_has_young", "no_young"))
        ) %>%
        mutate(LGT_similarity = replace_na(LGT_similarity, 0)) %>%
        ungroup()
    }
    
    # ----------------------------------------------------
    # 1. PAIRWISE WILCOXON TESTS
    # ----------------------------------------------------
    valid_groups <- sub_df %>%
      group_by(with_young) %>%
      summarise(n = sum(!is.na(LGT_similarity)), .groups = "drop") %>%
      filter(n > 0) %>%
      pull(with_young)
    
    if(current_group == "same_household"){
      valid_groups <- valid_groups[!grepl("one_has_young", valid_groups)]
    }
    
    df_brackets <- tibble()
    
    if (length(valid_groups) >= 2 && length(unique(na.omit(sub_df$LGT_similarity))) >= 2) {
      
      pw_res <- pairwise.wilcox.test(
        x = sub_df$LGT_similarity,
        g = sub_df$with_young,
        p.adjust.method = "none",
        #alternative = "less",
        exact = FALSE
      )
      
      p_mat <- pw_res$p.value
      comp_list <- list()
      
      for (r in rownames(p_mat)) {
        for (c in colnames(p_mat)) {
          if (!is.na(p_mat[r, c])) {
            comp_list[[length(comp_list) + 1]] <- tibble(
              group1    = c,
              group2    = r,
              p_val_raw = p_mat[r, c]
            )
          }
        }
      }
      
      if (length(comp_list) > 0) {
        df_brackets <- bind_rows(comp_list) %>%
          mutate(
            p_val_adj = p.adjust(p_val_raw, method = "fdr"),
            x1        = as.numeric(factor(group1, levels = levels(sub_df$with_young))),
            x2        = as.numeric(factor(group2, levels = levels(sub_df$with_young))),
            p_stars   = case_when(
              p_val_adj < 0.001 ~ "***",
              p_val_adj < 0.01  ~ "**",
              p_val_adj < 0.05  ~ "*",
              TRUE              ~ "ns"
            )
          )
      }
    }
    
    # ----------------------------------------------------
    # 2. BRACKET Y-POSITIONING
    # ----------------------------------------------------
    max_bar_height <- sub_df %>%
      group_by(with_young) %>%
      summarise(
        mean_val = mean(LGT_similarity, na.rm = TRUE),
        se_val   = sd(LGT_similarity, na.rm = TRUE) / sqrt(n()),
        top      = mean_val + se_val,
        .groups  = "drop"
      ) %>%
      pull(top) %>%
      max(na.rm = TRUE)
    
    if (is.infinite(max_bar_height) || is.na(max_bar_height)) max_bar_height <- 0.1
    
    step_height <- max_bar_height * 0.25
    base_y      <- max_bar_height * 1.08
    
    if (nrow(df_brackets) > 0) {
      df_brackets <- df_brackets %>%
        mutate(
          bracket_y = base_y + (row_number() - 1) * step_height,
          tip_y     = bracket_y - (step_height * 0.20),
          text_y    = bracket_y + (step_height * 0.35),
          x_mid     = (x1 + x2) / 2
        )
    }
    
    if(current_group == "same_household"){
      df_brackets <- df_brackets %>%
        filter(group1 != "one_has_young", group2 != "one_has_young")
    }
    
    # ----------------------------------------------------
    # 3. GGPLOT CONSTRUCTION
    # ----------------------------------------------------
    p <- ggplot(sub_df, aes(x = with_young, y = LGT_similarity, fill = sType)) +
      stat_summary(
        fun = "mean", 
        geom = "bar", 
        width = 0.95, 
        position = "identity"
      ) +
      stat_summary(
        fun.data = "mean_se", 
        geom = "errorbar", 
        width = 0.3, 
        linewidth=0.3,
        colour = "black", 
        position = "identity"
      ) +
      scale_fill_manual(values = site_col) +
      scale_x_discrete(
        labels = function(x) {
          paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
        }
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(subtitle = title_text, x = "", y = "LGT Similarity") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40")
      )
    
    if (nrow(df_brackets) > 0) {
      for (i in seq_len(nrow(df_brackets))) {
        b <- df_brackets[i, ]
        
        p <- p + 
          annotate(
            "path",
            x = c(b$x1, b$x1, b$x2, b$x2),
            y = c(b$tip_y, b$bracket_y, b$bracket_y, b$tip_y),
            colour = "black", linewidth = 0.3
          ) +
          annotate(
            "text",
            x = b$x_mid,
            y = b$text_y,
            label = b$p_stars,
            size = 3.5,
            vjust = 0.5
          )
      }
    }
    
    title_text <- paste0(current_host, current_type)
    
    if(current_group=="same_household"){
      p <- p+ labs(title=title_text)
    }
    
    p <- p + theme(legend.position = "none")
    return(p)
  })

#raw_plot_list

# Calculate summary metrics from data
nSample_retained <- length(unique(c(sLGT_HH$seqID_1, sLGT_HH$seqID_2)))
nPairs_retained  <- dplyr::n_distinct(sLGT_HH$seqID1_2)
nSandwich_retained <- list_lgt_by_type[grepl("DH", names(list_lgt_by_type))] %>%
  reduce(union) %>% length()

# Build multi-line caption message
filter_message <- paste0(
  "P-values from pairwise Wilcoxon Rank-Sum tests, not adjusted FDR adjustment:\n",
  "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001.\n\n",
  "Dataset filtering:\n",
  "(1) Initial set of samples falling within of 95% quantile range in LGT Richness per sample type (n = ", length(keep_sid), ");\n",
  "(2) Unique LGT set appearing in >3 samples per sample type (n = ", nSandwich_retained, ");\n",
  "(3) Samples containing at least one LGT from the above set, age <50 with relative data (n = ", nSample_retained, ");\n",
  "(4) Same-type pairwise combinations, e.g. Fecal-Fecal (n = ", nPairs_retained, ")."
)

# Assemble composite layout with annotations
combined_plot <- wrap_plots(raw_plot_list, nrow = 3) +
  plot_layout(
    guides = "collect",
    axis_titles   = "collect"
  ) +
  plot_annotation(
    title    = "Sample Similarity by LGT Profile by Presence of Young Relatives",
    subtitle = "(1 - Ruzicka Distance Metric)",
    caption  = filter_message
  ) &
  theme(
    # Titles & Captions
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.title    = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0, margin = margin(b = 10)),
    plot.caption  = element_text(size = 9, hjust = 0, color = "grey30", margin = margin(t = 12)),
    axis.text = element_text(size = 8),
    
    # Legend formatting
    #legend.title    = element_blank(),
    #legend.text     = element_text(size = 10),
    #legend.key.size = unit(5, "mm"),
    #legend.position = "bottom",
    #legend.box      = "horizontal"
    legend.position = "none"
  )

# Display plot
combined_plot

save_plot(combined_plot, "household_village_compare_ruzicka_dist_compare(1-D)_Young", height = 13, width = 15, unit = "cm", scale=1.5, dpi=600)

## LM 
sLGT_withyoung_fit <- list()
nSamples <- list()
for (i in unique(sLGT_HH$sType)) {
  sub_df <- subset(sLGT_HH, sType == i) %>%
    mutate(with_young = factor(with_young, levels = c("no_young", "one_has_young", "both_has_young"))) %>%
    rename(community = population) %>%
    left_join(unique_counts_long %>%
                select(seqID1_2, Unique_Count),
              by="seqID1_2") %>%
    left_join(metadata %>% select(sampleid, age) %>% rename(sample_1=sampleid, age_1=age),
              by="sample_1") %>%
    left_join(metadata %>% select(sampleid, age) %>% rename(sample_2=sampleid, age_2=age),
              by="sample_2") %>%
    rowwise() %>%
    mutate(mean_age=mean(c(age_1, age_2), na.rm=T)) %>%
    ungroup()
  
  nSamples[[as.character(i)]] <- length(union(sub_df$seqID_1, sub_df$seqID_2))
  
  fit <- lm(LGT_similarity ~ relation * with_young + log10(Unique_Count+1) + mean_age, data = sub_df)
  # fit <- lm(LGT_similarity ~ relation * with_young + log10(Unique_Count+1), data = sub_df)
  # fit <- lm(LGT_similarity ~ relation * with_young, data = sub_df)
  
  sLGT_withyoung_fit[[as.character(i)]] <- fit
  
  print(i)
  print(summary(fit))
}

sLGT_youngPairs_coef_list <- lapply(names(sLGT_withyoung_fit), function(model_name) {
  model <- sLGT_withyoung_fit[[model_name]]
  mod_sum <- summary(model)
  
  df <- data.frame(type=model_name, as.data.frame(coef(mod_sum)))
  df$term <- rownames(df)
  rownames(df) <- NULL
  
  p_col <- grep("Pr", colnames(df), value = TRUE)[1]
  
  if (!is.null(p_col) && !is.na(p_col)) {
    df$p.stars <- symnum(
      df[[p_col]],
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  } else {
    df$p.stars <- " "
  }
  
  df$r_squared <- round(mod_sum$r.squared, 3)
  df$adj_r_squared <- round(mod_sum$adj.r.squared, 3)
  df$nPairs <- nobs(model)
  df$nSample <- nSamples[[model_name]]
  
  return(df)
})

sLGT_youngPairs_coef_df <- do.call(rbind, sLGT_youngPairs_coef_list) %>%
  mutate(term = gsub("relation|with_young", "", term))
sLGT_youngPairs_coef_df

out_path <- file.path(basedir, "summary_tables/lm_sLGT_youngPairs.tsv")

writeLines("# LM model : LGT_similarity ~ relation * with_young + log10(Unique_Count+1) + pair_mean_age; filtered for Age < 50", con = out_path)

write.table(
  sLGT_youngPairs_coef_df,
  file = out_path,
  append = TRUE,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

sapply(sLGT_withyoung_fit, nobs)

sLGT_HH %>%
  group_by(type) %>%
  summarise(nSample = n_distinct(c(sample_1, sample_2)), .groups = "drop")

