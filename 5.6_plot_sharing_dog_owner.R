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

type_col <- c(
  "Fecal"="#56B4E9",
  "Oral"="#EFC000", 
  "Skin"="#EE2C2C"
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
# initial filter by select pairs_df
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


# # initial filter: select pairs
selected_pairs <- pairs_df %>%
  filter(host_pair %in% c("human_dog","human_human", "dog_dog")) %>%
  filter(relation %in% c("same_household","same_village")) %>%
  filter(type != "MIXED")

table(selected_pairs$host_pair, selected_pairs$relation)
table(selected_pairs$type, selected_pairs$relation)
table(selected_pairs$type, selected_pairs$host_pair)

select_pairs <- unique(selected_pairs$seqID1_2) %>% as.character() %>% sort()
select_ids_pairs <- unique(c(selected_pairs$seqID_1, selected_pairs$seqID_2)) %>% as.character() %>% sort()

keep_sid <- keep_sid[keep_sid %in% select_ids_pairs]

writeLines(filtered_sid,file.path(input_dir,"sid_failed_filter.txt"))
writeLines(keep_sid,file.path(input_dir,"sid_passing_filter.txt"))

table(count_df$type[count_df$SAMPLE %in% keep_sid])


# filter from all
count_df <- count_df %>% filter(original_sid %in% keep_sid) %>%
  mutate(type=gsub("human|dog","",type))
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

convert_to_presence <- function(species_matrix) {
  x <- species_matrix
  x[x > 0] <- 1
  return(x)
}

my_matrix <- LGTmat[sort(rownames(LGTmat)), ]
sample_types <- c("DHF", "DDF", "DHV", "DDV", "DHP", "DDH")

kept_cols_list <- lapply(sample_types, function(x) {
  get_samples <- substr(rownames(my_matrix), 1, 3)
  selected <- grepl(pattern = x, get_samples)
  mat <- my_matrix[selected, , drop = FALSE]
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
  filter(seqID1_2 %in% select_pairs) %>%
  filter(seqID_1 != seqID_2) %>% 
  filter(sample_1 != sample_2) %>% 
  mutate(
    host = ifelse(host_pair == "human_dog", "MIXED", host_pair),
    host = ifelse(host != "MIXED", sub("_.*", "", host), host),
    host = factor(host, levels = c("human", "dog", "MIXED")),
    relation2 = paste(host, relation, sep = "_")
  ) %>%
  mutate(host=factor(host,c("human","dog","MIXED")),
         relation=factor(relation, c("same_household","same_village"))) %>%
  arrange(host, relation) %>%
  mutate(relation2 = factor(relation2, unique(relation2))) %>%
  arrange(relation2)

table(cleaned_df$relation2)

table(cleaned_df$relation2, cleaned_df$type)

## bar
raw_plot_list <- cleaned_df %>%
  group_split(type) %>%
  map(function(sub_df) {
    sub_df <- sub_df %>% filter(sub_df$host != "MIXED")
    if (nrow(sub_df) == 0) return(NULL)
    
    current_type <- as.character(sub_df$type[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " (n = ", n_unique, ")")
    
    sub_df$relation <- factor(sub_df$relation, levels = c("same_household", "same_village"))
    sub_df$host <- droplevels(factor(sub_df$host, levels = c("human", "dog", "MIXED")))
    sub_df$relation2 <- droplevels(sub_df$relation2, levels = levels(cleaned_df$relation2))
    
    nPairs <- table(sub_df$relation2)
    
    if (any(nPairs == 0) || length(unique(sub_df$LGT_similarity[!is.na(sub_df$LGT_similarity)])) < 2) {
      return(NULL)
    } else {
      p <- ggplot(sub_df, aes(x = relation2, y = LGT_similarity, fill = type)) +
        facet_grid(host ~ type, space = "free", scales = "free") +
        stat_summary(fun = "mean", geom = "bar", width = 0.9) +
        stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
        scale_fill_manual(values = type_col) +
        scale_x_discrete(expand = expansion(mult = c(0.05, 0.1)),
                         labels = function(x) {
                           paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
                         }) +
        labs(title = title_text) +
        coord_flip() +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(size = 12, face = "bold"),
              panel.spacing = unit(1, units = "cm"))
      
      unique_hosts <- unique(sub_df$host)
      
      all_pval <- do.call(bind_rows, lapply(unique_hosts, function(x) {
        test_df <- subset(sub_df, host == x)
        
        n_pair <- table(test_df$relation)
        if (length(n_pair) < 2 || any(n_pair == 0)) return(NULL)
        
        optimal_multiplier <- round((max(n_pair) / min(n_pair)) * 100)
        
        test_df <- test_df %>%
          group_by(relation) %>%
          mutate(
            weight_int = ifelse(relation == "same_household", optimal_multiplier, 1)
          ) %>%
          ungroup()
        
        test_res <- coin::independence_test(
          LGT_similarity ~ relation,
          data = test_df, alternative = "greater",
          weights = ~ weight_int, 
          distribution = coin::approximate(nresample = 999)
        )
        
        print(paste(current_type, x, sep = " - "))
        print(test_res)
        p_val <- round(as.numeric(coin::pvalue(test_res)), 10)
        p_val <- ifelse(p_val < 0.001, 0.0009, p_val)
        data.frame(host = x, type = current_type, p_val = p_val)
      }))
      
      if (is.null(all_pval) || nrow(all_pval) == 0) return(p)
      
      all_pval$host <- droplevels(factor(all_pval$host, levels = c("human", "dog", "MIXED")))
      all_pval$type <- factor(all_pval$type, levels = levels(factor(sub_df$type)))
      
      annotation_df <- sub_df %>%
        group_by(host, type) %>%
        summarise(
          se_max = if(any(relation == "same_household")) {
            ggplot2::mean_se(LGT_similarity[relation == "same_household"])$ymax
          } else 0,
          .groups = "drop"
        ) %>%
        left_join(
          sub_df %>%
            group_by(type) %>%
            summarise(global_max = max(LGT_similarity, na.rm = TRUE), .groups = "drop"),
          by = "type"
        ) %>%
        left_join(all_pval, by = c("host", "type")) %>%
        filter(!is.na(p_val)) %>%
        mutate(
          se_max = ifelse(is.infinite(se_max) | is.na(se_max), 0, se_max),
          bracket_y = se_max + (global_max * 0.01),
          tip_y = bracket_y - (global_max * 0.005),
          text_y = bracket_y + (global_max * 0.01),
          relation2 = paste(host, "same_household", sep = "_"),
          p_adj = p.adjust(p_val, method = "BH"),
          p_stars = ifelse(p_adj < 0.001, "***",
                           ifelse(p_adj < 0.01,  "**",
                                  ifelse(p_adj < 0.05,  "*", "ns")))
        )
      
      if (nrow(annotation_df) > 0) {
        p <- p + 
          geom_segment(data = annotation_df, aes(x = 1, xend = 1, y = tip_y, yend = bracket_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_segment(data = annotation_df, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_segment(data = annotation_df, aes(x = 2, xend = 2, y = bracket_y, yend = tip_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_text(data = annotation_df, aes(x = 1.5, y = text_y, label = p_stars, fill = NULL), size = 5, vjust = 0.3, angle = 0) 
      }
      
      return(p)
    }
  })
    
raw_plot_list


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
    strip.text = element_blank(),
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

save_plot(combined_plot, "household_village_compare_ruzicka_dist_compare(1-D)", height = 20, width = 15, unit = "cm", scale=1.5)
saveRDS(raw_plot_list,  file.path(figpath,"household_village_compare_ruzicka_dist_compare(1-D).rds"))




#######################
# by Dog Ownership
######################
socials_dogs <- read.delim("input_file/socials/dog-owner_pairs.tsv") %>% 
  mutate(owner_id = gsub("-", "", owner_id), dog_id = gsub("-", "", dog_id)) %>%
  filter(grepl(";", dog_id)) %>%
  select(-household, -address_original, -population) %>%
  distinct() %>%
  mutate(
    owner_id = strsplit(as.character(owner_id), ";"),
    dog_id   = strsplit(as.character(dog_id), ";")
  ) %>%
  unnest(owner_id) %>%
  mutate(owner_id = trimws(owner_id)) %>%
  unnest(dog_id) %>%
  mutate(dog_id = trimws(dog_id)) %>%
  mutate(
    sample1_2 = paste(owner_id, dog_id, sep = "_"),
    relation  = "same_household",
    host_pair = "human_dog",
    population = substr(address, 1, 3)
  ) %>%
  distinct()

owners_id <- unique(socials_dogs$owner_id)

cleaned_df_byDogs <- cleaned_df %>%
  filter(host=="human") %>%
  select(-relation2) %>%
  mutate(with_dogs=ifelse(sample_1 %in% owners_id | sample_2 %in% owners_id, "withDogs", "noDogs")) %>%
  mutate(relation2=paste(relation, with_dogs, sep="_")) %>%
  mutate(type=factor(type, c("Fecal","Oral","Skin")),
         with_dogs=factor(with_dogs, c("withDogs","noDogs")),
         relation=factor(relation, c("same_household","same_village")))%>%
  arrange(type,relation) %>%
  mutate(relation2=factor(relation2, unique(relation2)))
  

table(cleaned_df_byDogs$relation2, cleaned_df_byDogs$type)

cleaned_df_byDogs <- cleaned_df %>%
  filter(host == "human") %>%
  select(-any_of("relation2")) %>%
  mutate(with_dogs = ifelse(sample_1 %in% owners_id | sample_2 %in% owners_id, "withDogs", "noDogs")) %>%
  mutate(relation2 = paste(relation, with_dogs, sep = "_")) %>%
  mutate(type = factor(type, c("Fecal", "Oral", "Skin")),
         with_dogs = factor(with_dogs, c("withDogs", "noDogs")),
         relation = factor(relation, c("same_household", "same_village"))) %>%
  arrange(type, relation) %>%
  mutate(relation2 = factor(relation2, unique(relation2)))

table(cleaned_df_byDogs$relation2, cleaned_df_byDogs$type)

raw_plot_list <- cleaned_df_byDogs %>%
  group_split(type) %>%
  map(function(sub_df) {
    if (nrow(sub_df) == 0) return(NULL)
    
    current_type <- as.character(sub_df$type[1])
    
    n_unique <- length(unique(c(sub_df$seqID_1, sub_df$seqID_2)))
    title_text <- paste0(current_type, " (n = ", n_unique, ")")
    
    sub_df$relation <- factor(sub_df$relation, levels = c("same_household", "same_village"))
    sub_df$host <- factor(sub_df$host, levels = c("human", "dog", "MIXED"))
    sub_df$host <- droplevels(sub_df$host)
    sub_df$relation2 <- factor(sub_df$relation2, levels = levels(cleaned_df_byDogs$relation2))
    sub_df$relation2 <- droplevels(sub_df$relation2)
    
    nPairs <- table(sub_df$relation2)
    
    if (any(nPairs == 0) || length(unique(sub_df$LGT_similarity[!is.na(sub_df$LGT_similarity)])) < 2) {
      return(NULL)
    } else {
      p <- ggplot(sub_df, aes(x = relation2, y = LGT_similarity, fill = type)) +
        facet_grid(with_dogs ~ type, space = "free", scales = "free") +
        stat_summary(fun = "mean", geom = "bar", width = 0.9) +
        stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
        scale_fill_manual(values = type_col) +
        scale_x_discrete(expand = expansion(mult = c(0.05, 0.1)),
                         labels = function(x) {
                           paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
                         }) +
        labs(title = title_text) +
        coord_flip() +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(size = 12, face = "bold"),
              panel.spacing = unit(1, units = "cm"))
      
      unique_dogs <- unique(sub_df$with_dogs)
      
      all_pval <- do.call(bind_rows, lapply(unique_dogs, function(x) {
        test_df <- subset(sub_df, with_dogs == x)
        
        n_pair <- table(test_df$relation)
        if (length(n_pair) < 2 || any(n_pair == 0)) return(NULL)
        
        optimal_multiplier <- round((max(n_pair) / min(n_pair)) * 100)
        
        test_df <- test_df %>%
          mutate(
            weight_int = ifelse(relation == "same_household", optimal_multiplier, 1)
          )
        
        test_res <- coin::independence_test(
          LGT_similarity ~ relation,
          data = test_df, alternative = "greater",
          weights = ~ weight_int, 
          distribution = coin::approximate(nresample = 999)
        )
        
        print(paste(current_type, x, sep = " - "))
        print(test_res)
        p_val <- round(as.numeric(coin::pvalue(test_res)), 10)
        p_val <- ifelse(p_val < 0.001, 0.0009, p_val)
        data.frame(host = as.character(test_df$host[1]), type = current_type, with_dogs = x, p_val = p_val)
      }))
      
      if (is.null(all_pval) || nrow(all_pval) == 0) return(p)
      
      all_pval$host <- factor(all_pval$host, levels = c("human", "dog", "MIXED"))
      all_pval$host <- droplevels(all_pval$host)
      all_pval$type <- factor(all_pval$type, levels = levels(factor(sub_df$type)))
      
      annotation_df <- sub_df %>%
        group_by(host, type, with_dogs) %>%
        summarise(
          se_max = if(any(relation == "same_household")) {
            ggplot2::mean_se(LGT_similarity[relation == "same_household"])$ymax
          } else 0,
          .groups = "drop"
        ) %>%
        left_join(
          sub_df %>%
            group_by(type) %>%
            summarise(global_max = max(LGT_similarity, na.rm = TRUE), .groups = "drop"),
          by = "type"
        ) %>%
        left_join(all_pval, by = c("host", "type", "with_dogs")) %>%
        filter(!is.na(p_val)) %>%
        mutate(
          se_max = ifelse(is.infinite(se_max) | is.na(se_max), 0, se_max),
          bracket_y = se_max + (global_max * 0.01),
          tip_y = bracket_y - (global_max * 0.005),
          text_y = bracket_y + (global_max * 0.01),
          relation2 = paste(host, "same_household", sep = "_"),
          p_adj = p.adjust(p_val, method = "BH"),
          p_stars = ifelse(p_adj < 0.001, "***",
                           ifelse(p_adj < 0.01,  "**",
                                  ifelse(p_adj < 0.05,  "*", 
                                         ifelse(p_adj < 0.1,  "m", "ns"))))
        )
      
      if (nrow(annotation_df) > 0) {
        p <- p + 
          geom_segment(data = annotation_df, aes(x = 1, xend = 1, y = tip_y, yend = bracket_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_segment(data = annotation_df, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_segment(data = annotation_df, aes(x = 2, xend = 2, y = bracket_y, yend = tip_y, fill = NULL), colour = "black", linewidth = 0.5) +
          geom_text(data = annotation_df, aes(x = 1.5, y = text_y, label = p_stars, fill = NULL), size = 5, vjust = 0.3, angle = 0) 
      }
      
      return(p)
    }
  })

raw_plot_list

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
    strip.text = element_blank(),
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

save_plot(combined_plot, "household_village_compare_wDogs_ruzicka_dist_compare(1-D)", height = 20, width = 15, unit = "cm", scale=1.5)
saveRDS(raw_plot_list, file.path(figpath,"household_village_compare_wDogs_ruzicka_dist_compare(1-D).rds"))

