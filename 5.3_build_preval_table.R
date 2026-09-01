rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# Paths
# my_lib <- "/home/caf77/R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))
my_lib <- "C:/Users/caf77_Local/Documents/Imaging_PC1_Rlib/R-4.3.0/library"
.libPaths(c(my_lib, .libPaths()))

# basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir=file.path(basedir, "index")
summary_dir=file.path(basedir, "summary_tables")

out_dir=file.path(basedir, "summary_tables")
# dir.create(out_dir,showWarnings = F,recursive = T)

figpath <- file.path(basedir, "figures")
#dir.create(figpath,showWarnings = F,recursive = T)

waafle_file="all_samples_internal_merged_indexed_filtered.tsv"


################
## Load Data  ##
################
type <- c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")
indata <- read.delim(file.path(summary_dir,waafle_file))

input_contigs <- read.delim(file.path(input_dir, "input_contigs_compile.tsv")) %>%
  mutate(SAMPLE=substr(file_name,1,6)) %>%
  group_by(SAMPLE) %>%
  summarise(n_contig=sum(n_contig))

input_contigs_sid <- read.delim(file.path(input_dir, "input_contigs_compile.tsv")) %>%
  mutate(SAMPLE=substr(file_name,1,6)) %>%
  group_by(file_name) %>%
  summarise(n_contig=sum(n_contig)) %>%
  rename(original_sid = file_name)

pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")

pop.df <- read.delim(file.path(input_dir, "Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(sampleid, 1, 6)) %>%
  select(SAMPLE, population, individual) %>%
  distinct() %>%
  mutate(group = case_when(
    population %in% c("BTU", "ORT", "ORS") ~ "Early-transition",
    population %in% c("APT", "TBU", "BSP") ~ "Late-transition",
    TRUE ~ "Agriculture"
  )) %>%
  filter(SAMPLE %in% indata$SAMPLE) %>%
  mutate(host = ifelse(grepl("x", individual), "dog", "human")) %>%
  mutate(
    host = factor(host, levels = c("human", "dog")),
    population = factor(population, levels = pop.ord)
  ) %>%
  arrange(host, population) %>%
  mutate(individual = factor(individual, levels = unique(individual)))

table(pop.df$group)

######################
## Cargo genes List ##
######################

# formula if replicates are averaged
df0 <- indata %>%
  filter(DIRECTION=="B>A") %>%
  group_by(B_genes) %>%
  summarise(
    nContig = {
      sum(vapply(unique(SAMPLE), function(x) {
        n_distinct(no[SAMPLE == x]) / n_distinct(original_sid[SAMPLE == x])
        }, FUN.VALUE = numeric(1)))
      },
    .groups = "drop") %>%
  arrange(desc(nContig))


# the formula I ended up using
df <- indata %>%
  filter(DIRECTION=="B>A") %>%
  merge.data.frame(pop.df, by="SAMPLE") %>%
  group_by(B_genes) %>%
  summarise(
    nContig = n_distinct(no),
    nSample = n_distinct(SAMPLE),
    nRecipient = n_distinct(CLADE_A),
    nDonor = n_distinct(CLADE_B),
    nSandwich = n_distinct(sandwich),
    nSandDist = n_distinct(paste(sandwich, A_gapLen)),
    nIndv = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual)
      length(ind)
    },
    human = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[prefixes == "DH"])
      length(ind)
    },
    dog = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[prefixes == "DD"])
      length(ind)
    },
    DHF = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHF"]
      length(ind)
    },
    DHV = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHV"]
      length(ind)
    },
    DHP = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHP"]
      length(ind)
    },
    DDF = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDF"]
      length(ind)
    },
    DDV = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDV"]
      length(ind)
    },
    DDH = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDH"]
      length(ind)
    },
    hEarlyTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Early-transition" & prefixes == "DH"])
      length(ind)
    },
    hLateTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Late-transition"& prefixes == "DH"])
      length(ind)
    },
    hAgriculture = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Agriculture" & prefixes == "DH"])
      length(ind)
    },
    dEarlyTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Early-transition" & prefixes == "DD"])
      length(ind)
    },
    dLateTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Late-transition"& prefixes == "DD"])
      length(ind)
    },
    dAgriculture = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Agriculture" & prefixes == "DD"])
      length(ind)
    },
    .groups = "drop") %>%
  arrange(desc(nSample))

# View(df[grep("unknown",df$B_genes),])

# write df
write.table(df, file.path(out_dir, "Cargo_genes_preval.tsv"), 
            sep="\t", quote=F, row.names=F)


#########################
## Cargo genes by indv ##
#########################
ddf <- df$B_genes[df$nIndv > 1]
col_order <- levels(pop.df$individual)

pair_df <- indata %>%
  filter(DIRECTION=="B>A") %>%
  distinct(CONTIG_NAME,CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, SAMPLE, .keep_all = TRUE) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  filter(B_genes %in% ddf) %>%
  group_by(B_genes, individual) %>%
  summarise(
    count = length(unique(individual)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = individual, values_from = count, values_fill = 0) %>%
  arrange(desc(rowSums(.[, -1], na.rm = TRUE))) %>%
  select(B_genes, any_of(col_order))

pair_df_multisite <- indata %>%
  filter(DIRECTION=="B>A") %>%
  distinct(CONTIG_NAME,CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, SAMPLE, .keep_all = TRUE) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  filter(B_genes %in% ddf) %>%
  group_by(B_genes, individual) %>%
  summarise(
    count = length(unique(SAMPLE)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = individual, values_from = count, values_fill = 0) %>%
  arrange(desc(rowSums(.[, -1], na.rm = TRUE))) %>%
  select(B_genes, any_of(col_order))

Bgenes_bySample <- indata %>%
  filter(DIRECTION=="B>A") %>%
  merge(pop.df %>% 
          mutate(type = factor(substr(SAMPLE, 1, 3), levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH"))), by = "SAMPLE") %>%
  group_by(B_genes, individual) %>%
  summarise(
    count = length(unique(SAMPLE)),
    hostype = {
      vec <- unique(substr(SAMPLE,1,3))
      paste(vec, collapse = ";")
    },
    .groups = "drop"
  ) %>%
  arrange(desc(count))

write.table(pair_df, file.path(out_dir, "Cargo_genes_pairs_anysite.tsv"), 
            sep="\t", quote=F, row.names=F)
write.table(pair_df_multisite, file.path(out_dir, "Cargo_genes_pairs_multisite.tsv"), 
            sep="\t", quote=F, row.names=F)
write.table(Bgenes_bySample, file.path(out_dir, "Cargo_genes_indv_count.tsv"), 
            sep="\t", quote=F, row.names=F)


#####################
## Unique LGT List ##
#####################
df2 <- indata %>%
  mutate(sandwichDist=paste(sandwich,A_gapLen, sep="_")) %>%
  merge.data.frame(pop.df, by="SAMPLE") %>%
  group_by(sandwichDist) %>%
  summarise(
    nContig = n_distinct(no),
    nSample = n_distinct(SAMPLE),
    nRecipient = n_distinct(CLADE_A),
    nDonor = n_distinct(CLADE_B),
    nIndv = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual)
      length(ind)
    },
    human = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[prefixes == "DH"])
      length(ind)
    },
    dog = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[prefixes == "DD"])
      length(ind)
    },
    DHF = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHF"]
      length(ind)
    },
    DHV = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHV"]
      length(ind)
    },
    DHP = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DHP"]
      length(ind)
    },
    DDF = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDF"]
      length(ind)
    },
    DDV = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDV"]
      length(ind)
    },
    DDH = {
      prefixes <- substr(SAMPLE, 1, 3)
      ind <- prefixes[prefixes == "DDH"]
      length(ind)
    },
    hEarlyTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Early-transition" & prefixes == "DH"])
      length(ind)
    },
    hLateTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Late-transition"& prefixes == "DH"])
      length(ind)
    },
    hAgriculture = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Agriculture" & prefixes == "DH"])
      length(ind)
    },
    dEarlyTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Early-transition" & prefixes == "DD"])
      length(ind)
    },
    dLateTransition = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Late-transition"& prefixes == "DD"])
      length(ind)
    },
    dAgriculture = {
      prefixes <- substr(SAMPLE, 1, 2)
      ind <- unique(individual[group == "Agriculture" & prefixes == "DD"])
      length(ind)
    }
  ) %>%
  arrange(desc(nSample))


# View(df[grep("unknown",df$B_genes),])

# write df
write.table(df2, file.path(out_dir, "Unique_LGT_preval.tsv"), 
            sep="\t", quote=F, row.names=F)

##########################
## Unique LGT by indv ##
##########################
ddf <- df2$sandwichDist[df2$nIndv > 1]
col_order <- levels(pop.df$individual)

pair_df <- indata %>%
  distinct(CONTIG_NAME,CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, SAMPLE, .keep_all = TRUE) %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep="_")) %>%
  merge(pop.df, by = "SAMPLE") %>%
  filter(sandwichDist %in% ddf) %>%
  group_by(sandwichDist, sandwich, A_gapLen, individual) %>%
  summarise(
    count = length(unique(individual)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = individual, values_from = count, values_fill = 0) %>%
  mutate(row_sum = rowSums(across(-c(sandwichDist, sandwich, A_gapLen)), na.rm = TRUE)) %>%
  arrange(desc(row_sum)) %>%
  select(sandwichDist, sandwich, A_gapLen, any_of(col_order))

pair_df_multisite <-  indata %>%
  distinct(CONTIG_NAME,CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, SAMPLE, .keep_all = TRUE) %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep="_")) %>%
  merge(pop.df, by = "SAMPLE") %>%
  filter(sandwichDist %in% ddf) %>%
  group_by(sandwichDist, sandwich, A_gapLen, individual) %>%
  summarise(
    count = length(unique(SAMPLE)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = individual, values_from = count, values_fill = 0) %>%
  mutate(row_sum = rowSums(across(-c(sandwichDist, sandwich, A_gapLen)), na.rm = TRUE)) %>%
  arrange(desc(row_sum)) %>%
  select(sandwichDist, sandwich, A_gapLen, any_of(col_order))

LGT_bySample <- indata %>%
  mutate(sandwichDist=paste(sandwich,A_gapLen, sep="_")) %>%
  merge(pop.df %>% 
          mutate(type = factor(substr(SAMPLE, 1, 3), levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH"))), by = "SAMPLE") %>%
  group_by(sandwichDist, sandwich, A_gapLen, individual) %>%
  summarise(
    count = length(unique(SAMPLE)),
    hostype = {
      vec <- unique(substr(SAMPLE,1,3))
      paste(vec, collapse = ";")
    },
    .groups = "drop"
  ) %>%
  arrange(desc(count))

write.table(pair_df, file.path(out_dir, "Unique_LGT_pairs_anysite.tsv"), 
            sep="\t", quote=F, row.names=F)
write.table(pair_df_multisite, file.path(out_dir, "Unique_LGT_pairs_multisite.tsv"), 
            sep="\t", quote=F, row.names=F)
write.table(LGT_bySample, file.path(out_dir, "Unique_LGT_indv_count.tsv"), 
            sep="\t", quote=F, row.names=F)


###########################
## Unique genes by sampleID #
###########################
ddfB <- df$B_genes[df$nSample > 1]
col_order <- levels(pop.df$individual)

pair_df_multisite <- indata %>%
  filter(DIRECTION=="B>A") %>%
  select(-type) %>%
  left_join(pop.df, by = "SAMPLE", keep=FALSE) %>%
  mutate(population = factor(population, pop.ord)) %>%
  mutate(original_sid = factor(original_sid, unique(original_sid))) %>%
  distinct(CONTIG_NAME, CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, original_sid, .keep_all = TRUE) %>%
  filter(B_genes %in% ddfB) %>%
  group_by(B_genes, original_sid) %>%
  summarise(
    count = n_distinct(original_sid),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = original_sid, values_from = count, values_fill = 0) %>%
  mutate(row_sum = rowSums(across(-c(B_genes)), na.rm = TRUE)) %>%
  arrange(desc(row_sum)) %>%
  select(-row_sum)

write.table(pair_df_multisite, file.path(out_dir, "Cargo_genes_pairs_multisite_bySampleID.tsv"), 
            sep="\t", quote=F, row.names=F) 


ddf <- df2$sandwichDist[df2$nSample > 1]
col_order <- levels(pop.df$individual)

pair_df_multisite <- indata %>%
  select(-type) %>%
  left_join(pop.df, by = "SAMPLE", keep=FALSE) %>%
  mutate(population = factor(population, pop.ord)) %>%
  mutate(original_sid = factor(original_sid, unique(original_sid))) %>%
  distinct(CONTIG_NAME, CONTIG_LENGTH, SYNTENY, Gene_SYNTENY, original_sid, .keep_all = TRUE) %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  filter(sandwichDist %in% ddf) %>%
  group_by(sandwichDist, sandwich, A_gapLen, original_sid) %>%
  summarise(
    count = n_distinct(original_sid),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = original_sid, values_from = count, values_fill = 0) %>%
  mutate(row_sum = rowSums(across(-c(sandwichDist, sandwich, A_gapLen)), na.rm = TRUE)) %>%
  arrange(desc(row_sum)) %>%
  select(-row_sum)

write.table(pair_df_multisite, file.path(out_dir, "Unique_LGT_pairs_multisite_bySampleID.tsv"), 
            sep="\t", quote=F, row.names=F)



###########################
## By total input Contig ##
###########################
pop.df <- pop.df %>% 
  mutate(
    type = factor(substr(SAMPLE, 1, 3), levels = c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH")),
    site = case_when(
      type %in% c("DHF","DDF") ~ "Fecal",
      type %in% c("DHV","DDV") ~ "Oral",
      type %in% c("DHP","DDH") ~ "Skin",
      TRUE ~ NA_character_
    )) %>%
  mutate(site=factor(site, c("Fecal","Oral","Skin")))

#global singletons across all samples
global_singletons <- indata %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(sandwichDist) %>%
  filter(n_distinct(SAMPLE) == 1) %>%
  pull(sandwichDist) %>%
  unique()

df3 <- indata %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(SAMPLE) %>%
  summarise(
    count           = n_distinct(CONTIG_NAME, CONTIG_LENGTH, original_sid),
    count_unique    = n_distinct(sandwichDist),
    count_singleton = n_distinct(sandwichDist[sandwichDist %in% global_singletons]),
    .groups = "drop"
  ) %>%
  merge.data.frame(input_contigs, by = "SAMPLE") %>%
  rename(total_contig = n_contig) %>%
  mutate(pLGT = count / total_contig,
         singleton_frac = count_singleton/count) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  arrange(desc(pLGT))
  
write.table(df3, file.path(out_dir, "nLGT_by_sample.tsv"), 
            sep="\t", quote=F, row.names=F)


df3 <- indata %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(original_sid) %>%
  summarise(
    count           = n_distinct(CONTIG_NAME, CONTIG_LENGTH, original_sid),
    count_unique    = n_distinct(sandwichDist),
    count_singleton = n_distinct(sandwichDist[sandwichDist %in% global_singletons]),
    .groups = "drop"
  ) %>%
  merge.data.frame(input_contigs_sid, by = "original_sid") %>%
  rename(total_contig = n_contig) %>%
  mutate(pLGT = count / total_contig,
         singleton_frac = count_singleton/count) %>%
  mutate(SAMPLE=substr(original_sid,1,6)) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  arrange(desc(pLGT))

write.table(df3, file.path(out_dir, "nLGT_by_seqID.tsv"), 
            sep="\t", quote=F, row.names=F)


df4 <- df3 %>%
  group_by(host, site, type) %>%
  summarise(
    n=n_distinct(SAMPLE),
    nContig_mio=sum(total_contig)/1e6,
    mean_count=mean(count),
    se_count=sd(count)/sqrt(n_distinct(SAMPLE)),
    min_count=min(count),
    max_count=max(count),
    mean_pLGT=mean(pLGT),
    se_pLGT=sd(pLGT)/sqrt(n_distinct(SAMPLE)),
    min_pLGT=min(pLGT),
    max_pLGT=min(pLGT),
  )

write.table(df4, file.path(out_dir, "nLGT_by_type.tsv"), 
            sep="\t", quote=F, row.names=F)


df5 <- indata %>%
  select(-type) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(sandwichDist) %>%
  summarise(
    nSample = n_distinct(SAMPLE),
    .groups = "drop"
  ) %>%
  mutate(singletons = case_when(
    nSample == 1 ~ "singleton",
    nSample > 1  ~ "multi"
  )) %>%
  arrange(desc(nSample)) %>%
  summarise(
    singleton = sum(singletons == "singleton", na.rm = TRUE),
    multi     = sum(singletons == "multi", na.rm = TRUE)
  )

df6 <- indata %>%
  select(-type) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(sandwichDist, type) %>%
  summarise(
    nSample = n_distinct(SAMPLE),
    .groups = "drop"
  ) %>%
  mutate(singletons = case_when(
    nSample == 1 ~ "singleton",
    nSample > 1  ~ "multi"
  )) %>%
  arrange(desc(nSample)) %>%
  group_by(type) %>%
  summarise(
    singleton = sum(singletons == "singleton", na.rm = TRUE),
    multi     = sum(singletons == "multi", na.rm = TRUE)
  )


df7 <- indata %>%
  select(-type) %>%
  merge.data.frame(pop.df, by = "SAMPLE") %>%
  mutate(sandwichDist = paste(sandwich, A_gapLen, sep = "_")) %>%
  group_by(sandwichDist, host) %>%
  summarise(
    nSample = n_distinct(SAMPLE),
    .groups = "drop"
  ) %>%
  mutate(singletons = case_when(
    nSample == 1 ~ "singleton",
    nSample > 1  ~ "multi"
  )) %>%
  arrange(desc(nSample)) %>%
  group_by(host) %>%
  summarise(
    singleton = sum(singletons == "singleton", na.rm = TRUE),
    multi     = sum(singletons == "multi", na.rm = TRUE)
  )

df8 <- bind_rows(df5 %>% mutate(group="ALL"), df6 %>% rename(group=type), df7 %>% rename(group=host)) %>%
  mutate(singleton_frac=singleton/(singleton+multi))

write.table(df8, file.path(out_dir, "singleton_fraction_byGroup.tsv"), 
            sep="\t", quote=F, row.names=F)


