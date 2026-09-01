rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)

# Paths
# my_lib <- "R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))

basedir="Imaging_Lab_PC1/WAAFLE_Extra"
# basedir="Imaging_Lab_PC1/WAAFLE_Extra"
# setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir <- file.path(basedir, "index")

summary_dir <- file.path(basedir, "summary_tables")

figpath <- file.path(basedir, "figures")
dir.create(figpath,showWarnings = F,recursive = T)

waffle_file="all_samples_internal_merged_indexed_filtered.tsv"
dist_df="dist/ruzicka_uniqueLGT_multisite_bySampleID.rds"
pairs_count="Unique_LGT_pairs_multisite_bySampleID.tsv"
outlier="outliers_waafle.txt"
sPairs="all_household_village_seqID_pairs.tsv"

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

ru_dist <- readRDS(dist_df) %>% as.matrix()
sLGT <- 1-ru_dist
sLGT_clean <- sLGT[!rownames(sLGT) %in% outlier_ids, !colnames(sLGT) %in% outlier_ids]
sLGT_clean[1:5,1:5]
sLGT_clean <- as.dist(sLGT_clean)

uniqueLGT_list <-  read.delim(file.path(summary_dir, pairs_count)) %>%
  select(-any_of(outlier_ids))

original_sid <- colnames(as.matrix(sLGT_clean))

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
  )

pairs_df <- read.delim(file.path(input_dir, sPairs)) %>%
  mutate(across(everything(), ~ gsub("_POOLED", "", .x))) %>%
  mutate(
    relation  = factor(relation, c("replicate", "same_household", "same_village")),
    host_pair = factor(host_pair, c("human_human", "dog_dog", "human_dog")),
    type      = factor(type, c("Fecal", "Oral", "Skin"))
  ) %>%
  arrange(type, host_pair, relation, seqID1_2) %>%
  distinct(seqID_1, seqID_2, .keep_all = TRUE) %>%
  mutate(
    key_1     = ifelse(seqID_2 > seqID_1, seqID_1, seqID_2),
    key_2     = ifelse(seqID_2 > seqID_1, seqID_2, seqID_1),
    seqID1_2  = paste(key_1, key_2, sep = "_")
  ) %>%
  select(-c(seqID_1, seqID_2)) %>%
  rename(seqID_1 = key_1, seqID_2 = key_2) %>%
  mutate(
    key_3     = ifelse(sample_2 > sample_1, sample_1, sample_2),
    key_4     = ifelse(sample_2 > sample_1, sample_2, sample_1),
    sample1_2 = paste(key_3, key_4, sep = "_")
  ) %>%
  select(-c(sample_1, sample_2)) %>%
  rename(sample_1 = key_3, sample_2 = key_4) %>%
  filter(!is.na(relation)) %>%
  filter(seqID_1 != seqID_2)

unique(pairs_df$relation)
unique(pairs_df$host_pair)
unique(pairs_df$type)

dim(pairs_df)

#######################
## Compare by Pairs  ##
#######################
missing_seqID <- original_sid[!original_sid %in% unique(c(pairs_df$seqID_1, pairs_df$seqID_2))]

missing_pop_info <- pop.df %>%
  filter(original_sid %in% missing_seqID)

missing_pairs_constructed <- expand.grid(
  seqID_1 = missing_pop_info$original_sid,
  seqID_2 = missing_pop_info$original_sid,
  stringsAsFactors = FALSE
) %>%
  filter(seqID_2 > seqID_1) %>%
  filter(seqID_1 != seqID_2) %>%
  left_join(
    missing_pop_info %>% select(seqID_1 = original_sid, sample_1 = individual, pop_1 = population, host_1 = host, type_1 = type),
    by = "seqID_1"
  ) %>%
  left_join(
    missing_pop_info %>% select(seqID_2 = original_sid, sample_2 = individual, pop_2 = population, host_2 = host, type_2 = type),
    by = "seqID_2"
  ) %>%
  filter(pop_1 == pop_2) %>%
  filter(type_1 == type_2) %>%
  mutate(
    sample1_2 = paste(sample_1, sample_2, sep = "_"),
    seqID1_2  = paste(seqID_1, seqID_2, sep = "_"),
    relation  = "same_village",
    study     = "with_Pilot",
    population = pop_1,
    type       = type_1,
    host_pair  = paste(host_1, host_2, sep = "_")
  ) %>%
  select(colnames(pairs_df))

pairs_df_complete <- rbind(pairs_df, missing_pairs_constructed) %>%
  filter(seqID_1 %in% original_sid & seqID_2 %in% original_sid)

unrelated_pairs_constructed <- expand.grid(
  seqID_1 = pop.df$original_sid,
  seqID_2 = pop.df$original_sid,
  stringsAsFactors = FALSE
) %>%
  filter(seqID_2 > seqID_1) %>%
  filter(seqID_1 != seqID_2) %>%
  mutate(seqID1_2  = paste(seqID_1, seqID_2, sep = "_")) %>%
  left_join(
    pop.df %>% select(seqID_1 = original_sid, sample_1 = individual, pop_1 = population, host_1 = host, type_1 = type),
    by = "seqID_1"
  ) %>%
  left_join(
    pop.df %>% select(seqID_2 = original_sid, sample_2 = individual, pop_2 = population, host_2 = host, type_2 = type),
    by = "seqID_2"
  ) %>%
  filter(!seqID1_2 %in% pairs_df_complete$seqID1_2) %>%
  mutate(
    same_population = ifelse(pop_1 == pop_2, "same_pop", "FALSE"),
    same_individual = ifelse(sample_1 == sample_2, "same_indv", "FALSE"),
    same_type       = ifelse(type_1 == type_2, "same_type", "FALSE"),
  ) %>%
  mutate(host_pair=case_when(
    (host_1 == "human" | host_2 == "human") & host_1 == host_2 ~ "human_human",
    (host_1 == "dog" | host_2 == "dog") & host_1 == host_2 ~ "dog_dog",
    TRUE ~ "human_dog"
  )) %>%
  mutate(
    relation = case_when(
      same_population == "FALSE" & same_individual == "FALSE" ~ "unrelated",
      same_individual == "same_indv" & same_type == "FALSE"~ "same_person_different_type",
      same_individual == "same_indv" & same_type == "same_type"~ "replicate",
      same_population == "same_pop" & same_type == "same_type" ~ "same_village",
      same_population == "same_pop" & same_type == "FALSE" ~ "same_village_different_type",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(sample1_2  = paste(sample_1, sample_2, sep = "_"),
         population = case_when(
           same_population == "same_pop" ~ substr(sample_1,1,3),
           TRUE ~ "MIXED"
         ),
         type = case_when(
           same_type == "same_type" ~ gsub("human|dog","",type_1),
           TRUE ~ "MIXED"
         ),
         host_pair=paste(host_1, host_2, sep="_"),
         study="MOBILE_only") %>%
  arrange(same_population, same_individual, same_type) %>%
  select(all_of(colnames(pairs_df_complete)))


table(unrelated_pairs_constructed$relation)
table(unrelated_pairs_constructed$host_pair[unrelated_pairs_constructed$relation=="same_village"])

pairs_df_all_relations <- bind_rows(pairs_df_complete, unrelated_pairs_constructed) %>%
  distinct(seqID1_2, .keep_all = T) %>%
  filter(seqID_1 != seqID_2)


# get same household list
socials <- read.delim("input_file/socials/dog-owner_pairs.tsv") %>% 
  mutate(owner_id = gsub("-", "", owner_id), dog_id = gsub("-", "", dog_id)) %>%
  filter(grepl(";", dog_id))

household_dog_pairs <- socials %>%
  select(-owner_id, -household, -address_original, -population) %>%
  distinct() %>%
  mutate(dog_id = strsplit(as.character(dog_id), ";")) %>%
  unnest(dog_id) %>%
  inner_join(., ., by = "address", suffix = c("_1", "_2"), relationship = "many-to-many") %>%
  filter(dog_id_2 > dog_id_1) %>%
  filter(dog_id_1 != dog_id_2) %>%
  mutate(
    dog1_2    = paste(dog_id_1, dog_id_2, sep = "_"),
    relation  = "same_household",
    host_pair = "dog_dog",
    population = substr(address, 1,3)
  ) %>%
  distinct()

socials2 <- readRDS("input_file/socials/model_data_long_allPops.rds") %>%
  mutate(
    sample_1 = substr(sample1_2, 1, 6),
    sample_2 = substr(sample1_2, 8, 13)
  ) %>%
  mutate(
    key_1 = ifelse(sample_1 > sample_2, sample_1, sample_2),
    key_2 = ifelse(sample_1 > sample_2, sample_2, sample_1),
    sample1_2 = paste(key_1, key_2, sep = "_")
  ) %>%
  select(-c(sample_1, sample_2)) %>%
  rename(sample_1 = key_1, sample_2 = key_2) %>%
  filter(sample_1 != sample_2) %>%
  filter(House_share == 2)

household_human_dog_pairs <- socials %>%
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

same_household_pairs <- unique(c(household_dog_pairs$dog1_2, socials2$sample1_2, household_human_dog_pairs$sample1_2))

# fix pair errors
check_indv_pairs <- pairs_df_all_relations %>%
  filter(relation == "same_village") %>%
  select(sample1_2) %>%
  unlist() %>%
  unique() %>%
  as.character() %>%
  sort()

table(check_indv_pairs %in% same_household_pairs)  

to_fix1 <- check_indv_pairs[check_indv_pairs %in% same_household_pairs] # same_household


check_indv_pairs2 <- pairs_df_all_relations %>%
  filter(relation  == "same_village_different_type") %>%
  select(sample1_2) %>%
  unlist() %>%
  unique() %>%
  as.character() %>%
  sort()

table(check_indv_pairs2 %in% same_household_pairs)  

to_fix2 <- check_indv_pairs2[check_indv_pairs2 %in% same_household_pairs] # same_household_different_type

##

pairs_df_all_relations$relation   [pairs_df_all_relations$sample1_2 %in% to_fix1 &
                                     pairs_df_all_relations$relation == "same_village"] <- "same_household"

pairs_df_all_relations$relation   [pairs_df_all_relations$sample1_2 %in% to_fix2 &
                                     pairs_df_all_relations$relation == "same_village_different_type"] <- "same_household_different_type"


table(pairs_df_all_relations$relation)


write.table(pairs_df_all_relations, file.path(input_dir,"all_seqID_pairs_waafle.tsv"),quote = F,row.names = F,sep = "\t")
