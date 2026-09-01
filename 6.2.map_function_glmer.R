gc()
rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(lmerTest)

# basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

summary_dir=file.path(basedir, "summary_tables")
figpath=file.path(basedir, "figures")
out_dir=file.path(basedir, "summary_tables")

waafle_file <- "all_samples_internal_merged_indexed_filtered.tsv"

nContigs_file <- "input_file/input_contigs_compile.tsv"

cargo_preval <- "summary_tables/Cargo_genes_preval.tsv"
cargo_matrix <- "summary_tables/Cargo_genes_pairs_multisite_bySampleID.tsv"
LGT_preval <- "summary_tables/Unique_LGT_preval.tsv"
LGT_matrix <- "summary_tables/Unique_LGT_pairs_multisite_bySampleID.tsv"

top10_cargo <- "summary_tables/top_10_cargo_byContig.tsv"

annotation_file <- "input_file/annotations/consolidated_gene_annotation_filled_mapping.rds"
master_all_annot <- "input_file/annotations/master_all_annotations_filled_withGO.rds"
annot_details <- "input_file/annotations/gene_annotation_mapping_filled_draft.rds"
annot_GO <- "input_file/annotations/InterproID_to_GO_category.tsv"

# other annotation files, but use only to debug:
#external_annot <- "input_file/annotations/UPI_UniProt_info_unannotated_only.tsv"
#master_amr <- "input_file/annotations/master_amr.tsv"
#master_cazy <- "input_file/annotations/master_cazy.tsv"
#master_kofam <- "input_file/annotations/master_kofam_hits.tsv"

types <- c("humanFecal","humanOral","humanSkin","dogFecal","dogOral","dogSkin")
pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")

outliers <- "input_file/outliers_waafle.txt"

#######################
## Helper Functions  ##
#######################
save_plot <- function(p, name, width, height, unit="in", scale=1, dpi=150) {
  formats <- c("rds", "png", "pdf")
  for(fmt in formats){
    file <- file.path(figpath, paste0(name, ".", fmt))
    if(fmt == "rds") saveRDS(p, file)
    else ggsave(file, plot = p, width = width, height = height, dpi = dpi, unit=unit, scale=scale) 
  }
}

site_col <- c(
  "humanFecal"="#56B4E9", 
  "humanSkin"="#EE2C2C", 
  "humanOral"="#EFC000", 
  "dogFecal"="#3A5FCD", 
  "dogSkin"="#7D0226", 
  "dogOral"="#9A5324",
  "human" = "grey70",
  "dog" = "grey40",
  "All" = "grey10"
)

################
## Load Data  ##
################
outlier_ids <- read.delim(file.path(basedir,outliers))$sampleid

indata <- read.delim(file.path(summary_dir,waafle_file)) %>%
  rowwise() %>%
  mutate(n_Bgenes = length(unlist(strsplit(B_genes, split = "\\|")))) %>%
  ungroup() %>%
  filter(n_Bgenes <= 10 & nchar(SYNTENY) <= 35)

#top10_byContig <- read.delim(top10_cargo)

cargo.counts <- read.delim(file.path(basedir, cargo_preval))
cargo.mat <- read.delim(file.path(basedir, cargo_matrix), row.names = 1) %>%
  mutate(across(everything(), as.numeric)) %>%
  select(all_of(sort(names(.)))) %>%
  as.matrix()
cargo.mat <- cargo.mat[,!colnames(cargo.mat) %in% outlier_ids]

LGT.counts <- read.delim(file.path(basedir, LGT_preval))
LGT.mat <- read.delim(file.path(basedir, LGT_matrix), row.names = 1) %>%
  select(-sandwich, -A_gapLen) %>%
  mutate(across(everything(), as.numeric)) %>%
  select(all_of(sort(names(.)))) %>%
  as.matrix()
LGT.mat <- LGT.mat[,!colnames(LGT.mat) %in% outlier_ids]

match_ids <- data.frame(original_sid=colnames(cargo.mat), SAMPLE=substr(colnames(cargo.mat),1,6))

pop.df <- read.delim(file.path(basedir, "input_file/Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(seqID,1,6)) %>%
  select(SAMPLE, population, individual) %>%
  distinct() %>%
  right_join(match_ids, by="SAMPLE", keep=FALSE) %>%
  select(-SAMPLE) %>%
  rename(SAMPLE=original_sid) %>%
  distinct() %>%
  mutate(group = case_when(
    population %in% c("BTU", "ORT", "ORS") ~ "Early-transition",
    population %in% c("APT", "TBU", "BSP") ~ "Late-transition",
    TRUE ~ "Agriculture"
  )) %>%
  mutate(host = ifelse(grepl("x", individual), "dog", "human")) %>%
  mutate(
    host = factor(host, levels = c("human", "dog")),
    population = factor(population, levels = pop.ord)
  ) %>%
  arrange(host, population, SAMPLE) %>%
  mutate(individual = factor(individual, levels = unique(individual)))

nContigs <- read.delim(file.path(basedir, nContigs_file)) %>%
  rename(original_sid=file_name) %>%
  right_join(match_ids, by="original_sid", keep=FALSE) %>%
  select(-SAMPLE) %>%
  rename(SAMPLE=original_sid) %>%
  group_by(SAMPLE) %>%
  summarise(n_contig=round(mean(n_contig),0)) %>%
  distinct()

metadata <- read.delim("input_file/2023_MOBILE_AnthropometryData - CORRECTED.tsv") %>%
  rename(individual=sampleid) %>%
  left_join(pop.df %>% select(SAMPLE, individual) %>% distinct(),
             by="individual") %>%
  group_by(population) %>%
  mutate(age = ifelse(is.na(age), median(age, na.rm = TRUE), age),
         sex = ifelse(is.na(sex), names(which.max(table(sex))), sex)) %>%
  ungroup() %>%
  select(-population, -SAMPLE) %>%
  distinct()


annot_byGenes <- readRDS(file.path(basedir, annotation_file))
annot_bySAMPLE <- readRDS(master_all_annot)
annot_details <- readRDS(file.path(basedir, annot_details))


dim(cargo.mat)
dim(LGT.mat) # n sample is lower than cargo.mat, because of no singletons
dim(pop.df)
dim(nContigs)

types <- setNames(c("humanFecal","humanOral","humanSkin","dogFecal","dogOral","dogSkin"),
                     c("DHF","DHV","DHP","DDF","DDV","DDH"))

##############
## Analyse  ##
##############
pLGT <- indata %>%
  filter(!original_sid %in% outlier_ids) %>%
  mutate(
    B_genes = paste("len", A_gapLen, gsub("UniRef90_", "", B_genes), sep = "_"),
    joint_contig = paste(original_sid, CONTIG_NAME, sep = "+"),
    SAMPLE = original_sid,
    type = factor(as.character(types[substr(SAMPLE, 1, 3)]), levels = types)
  ) %>%
  group_by(SAMPLE, type) %>%
  summarise(
    count_LGT = n_distinct(joint_contig),
    .groups = "drop"
  ) %>%
  left_join(nContigs, by = "SAMPLE") %>%
  mutate(pLGT = 100 * count_LGT / n_contig)

pLGT_outliers <- pLGT %>%
  group_by(type) %>%
  mutate(
    q5 = quantile(pLGT, 0.025, na.rm = TRUE),
    q95 = quantile(pLGT, 0.975, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(pLGT < q5 | pLGT > q95) %>%
  pull(SAMPLE) %>%
  unique()

sample_uniqueLGT <- indata %>%
  filter(!original_sid %in% outlier_ids) %>%
  mutate(
    B_genes = paste("len", A_gapLen, gsub("UniRef90_", "", B_genes), sep = "_"),
    joint_contig = paste(original_sid, CONTIG_NAME, sep = "+"),
    SAMPLE = original_sid,
    type = factor(as.character(types[substr(SAMPLE, 1, 3)]), levels = types)
  ) %>%
  group_by(SAMPLE, type) %>%
  summarise(
    sample_totalLGT = n_distinct(B_genes),
    .groups = "drop"
  ) %>% 
  select (-type)

#### SINGLETONS
singletons <- indata %>%
  mutate(B_genes = paste("len", A_gapLen, gsub("UniRef90_", "", B_genes), sep = "_")) %>%
  group_by(B_genes) %>%
  summarise(nSample = n_distinct(original_sid), .groups = "drop") %>%
  mutate(group = case_when(
    nSample == 1 ~ "singletons",
    nSample == 2 ~ "doubletons",
    nSample > 2 ~ ">2 samples",
    TRUE ~ NA_character_
  ))

table(singletons$nSample)

p <- ggplot(singletons, aes(x = group)) +
  geom_bar(fill="coral") +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5, size=5) +
  scale_x_discrete(limits = c("singletons", "doubletons", ">2 samples")) +
  theme_bw(base_size = 14) +
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,75000))+
  labs(title="Unique LGT prevalence", y="# Unique LGT")

save_plot(p, "Unique_LGT_Prevalence", width = 8, height = 8, dpi = 300, scale=0.5)

singleton_list <- singletons %>%
  filter(group == "singletons") %>%
  pull(B_genes) %>%
  unique()

singletons <- indata %>%
  mutate(B_genes = paste("len", A_gapLen, gsub("UniRef90_", "", B_genes), sep = "_")) %>%
  mutate(type = case_when(
    grepl("DHF", SAMPLE) ~ "humanFecal",
    grepl("DHV", SAMPLE) ~ "humanOral",
    grepl("DHP", SAMPLE) ~ "humanSkin",
    grepl("DDF", SAMPLE) ~ "dogFecal",
    grepl("DDV", SAMPLE) ~ "dogOral",
    grepl("DDH", SAMPLE) ~ "dogSkin"
  )) %>% 
  mutate(
    host = ifelse(grepl("human", type), "human", "dog"),
    site = gsub("human|dog", "", type)
  ) %>%
  mutate(
    type = droplevels(factor(type, names(site_col))),
    host = factor(host, c("human", "dog")),
    site = factor(site, c("Fecal", "Oral", "Skin"))
  ) %>%
  group_by(type, host, site, B_genes) %>%
  summarise(nSample = n_distinct(original_sid), .groups = "drop") %>%
  mutate(group = case_when(
    nSample == 1 ~ "singletons",
    nSample == 2 ~ "doubletons",
    nSample > 2 ~ ">2 samples",
    TRUE ~ NA_character_
  ))

p <- ggplot(singletons, aes(x = group)) +
  facet_grid(site ~ host, scales = "free_y") +
  geom_bar(aes(fill = type), width = 0.9) +
  stat_count(geom = "text", aes(label = after_stat(count)), vjust = -0.5, size = 3.5) +
  scale_x_discrete(limits = c("singletons", "doubletons", ">2 samples")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face="bold", size=12),
    axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=12))+
  scale_fill_manual(values = site_col) +
  labs(title = "Unique LGT Prevalence", y = "# Unique LGT")

p

save_plot(p, "Unique_LGT_Prevalence_type", width = 15, height = 10, dpi = 300, scale=0.5)




############
## Cargo  ##
############
B_sigletons <- indata %>%  
  filter(!original_sid %in% outlier_ids) %>%
  filter(!original_sid %in% pLGT_outliers) %>%
  filter(DIRECTION == "B>A") %>%
  select(original_sid, B_genes) %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(B_genes) %>%
  summarise(nSample=n_distinct(original_sid), .groups = "drop") %>%
  arrange(desc(nSample)) %>%
  filter(nSample <= 1) %>%
  pull(B_genes) %>%
  unique()

B_total <- indata %>%
  filter(!original_sid %in% outlier_ids) %>%
  filter(!original_sid %in% pLGT_outliers) %>%
  select(original_sid, B_genes) %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(original_sid) %>%
  summarise(B_total = n_distinct(B_genes)) %>%
  arrange(desc(B_total)) %>%
  rename(SAMPLE=original_sid)

top_Bgenes <-  indata %>%  
  filter(!original_sid %in% outlier_ids) %>%
  filter(!original_sid %in% pLGT_outliers) %>%
  filter(DIRECTION == "B>A") %>%
  mutate(contig_join=paste(original_sid,CONTIG_NAME, sep="+")) %>%
  select(contig_join, type, B_genes) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(type, B_genes) %>%
  summarise(nContig = n_distinct(contig_join), .groups = "drop") %>%
  arrange(type, desc(nContig)) %>% 
  group_by(type) %>%
  slice_max(order_by = nContig, n = 10)

## Gene Function Proportions
View(annot_bySAMPLE)
View(indata)

df <- indata %>%
  filter(DIRECTION == "B>A", !original_sid %in% outlier_ids) %>%
  select(CONTIG_NAME, B_genes, B_genes, SAMPLE = original_sid) %>% 
  mutate(
    type = factor(types[substr(SAMPLE, 1, 3)], levels = types),
    index = row_number()
  ) %>%
  separate_rows(B_genes, sep = "\\|") %>%
  group_by(index) %>%
  mutate(
    g_index = row_number(),
    pair = paste(SAMPLE, CONTIG_NAME, B_genes, sep = "+")
  ) %>%
  ungroup()

to_annot <- annot_bySAMPLE %>%
  rename(B_genes=ANNOTATIONS.UNIREF90_fixed) %>%
  mutate(pair = paste(SAMPLE, CONTIG_NAME,B_genes, sep = "+")) %>%
  filter(pair %in% unique(df$pair)) %>%
  select(
    pair, is_amr, what_amr, 
    is_cazy, what_cazy, 
    ACCESSION, 
    Consolidated_Function,
    GO_BIOLOGICAL_PROCESS,
    GO_MOLECULAR_FUNCTION,
    GO_CELLULAR_COMPONENT
  ) %>%
  distinct()
  

annotated_df <- df %>%
  left_join(to_annot,
    by = "pair"
  ) %>%
  select(-pair)

trimm_df <- annotated_df %>%
  group_by(index, B_genes, CONTIG_NAME, SAMPLE, type) %>%
  summarise(
    is_amr = any(is_amr == TRUE, na.rm = TRUE),
    is_cazy = any(is_cazy == TRUE, na.rm = TRUE),
    what_gene = paste(na.omit(unique(ACCESSION)), collapse = "; "),
    Consolidated_Function = paste(na.omit(unique(Consolidated_Function)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    is_amr & is_cazy ~ "Any is both",
    is_amr ~ "Any is ARG",
    is_cazy ~ "Any is CAZymes",
    TRUE ~ "Other"
  )) %>%
  left_join(
    pop.df %>%
      bind_rows(
        data.frame(
          individual = c("APTDx011", "APTDx012", "TBUDx010"),
          SAMPLE = c("DDV031", "DDH032", "DDH062"),
          group = "Late-transition",
          host = "dog",
          population = c("APT", "APT", "TBU")
        )
      ),
    by = "SAMPLE"
  )

table(ifelse(trimm_df$is_amr,"AMR","Other"), ifelse(trimm_df$is_cazy,"CAZY","Other"))

LGT_by_types <- trimm_df %>%
  group_by(type) %>%
  summarise(total_LGT=n_distinct(index),
            any_ARG=100*sum(is_amr)/total_LGT,
            any_CAZyme=100*sum(is_cazy)/total_LGT,
            other=100 - any_ARG - any_CAZyme
            ) 

LGT_by_popTypes <- trimm_df %>%
  group_by(population, type) %>%
  summarise(total_LGT=n_distinct(index),
            any_ARG=100*sum(is_amr)/total_LGT,
            any_CAZyme=100*sum(is_cazy)/total_LGT,
            other=100 - any_ARG - any_CAZyme,
            .groups = "drop"
  )



prop_byType <- LGT_by_types %>%
  pivot_longer(
    cols = c(any_ARG, any_CAZyme, other),
    names_to = "category",
    values_to = "percentage"
  ) %>% 
  filter(category != "other") %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog"))) %>%
  ggplot(aes(x = category, y = percentage, fill = type)) +
  facet_grid(host~body_site, scale="free_y") +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = {
      lab=sprintf("%.2f%%", percentage);
      ifelse(lab=="0.00%","NA",lab)}),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  labs(
    y = "Percentage of LGT (%)",
    fill = "Category",
    caption = "an LGT can have both ARG and CAZymes"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values=site_col) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour="black"),
        strip.text = element_text(face="bold", size=10),
        panel.spacing.y = unit(1, "cm"))

prop_byType

save_plot(prop_byType, "LGT_function_prop_byType", width = 15, height = 15, unit = "cm",dpi=600, scale=1.2)



arg_byPop <- LGT_by_popTypes %>%
  pivot_longer(
    cols = c(any_ARG, any_CAZyme, other),
    names_to = "category",
    values_to = "percentage"
  ) %>% 
  filter(category == "any_ARG") %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog")),
         population=factor(population, pop.ord)) %>%
  ggplot(aes(x=population, y=percentage, fill=type)) +
  facet_grid(host~body_site, scale="free_y") +
  geom_col(position = "dodge") +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values=site_col) +
  labs(
    title = "Antimicrobial Resistance Genes",
    y = "Percentage of LGT with ARG (%)",
    fill = "Category",
    caption = "an LGT can have both ARG and CAZymes"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour="black"),
        strip.text = element_text(face="bold", size=10),
        panel.spacing.y = unit(1, "cm"))

arg_byPop

save_plot(arg_byPop, "LGT_function_prop_byPop_ARG", width = 15, height = 15,unit = "cm",dpi=600, scale=1.2)


cazy_byPop <-  LGT_by_popTypes %>%
  pivot_longer(
    cols = c(any_ARG, any_CAZyme, other),
    names_to = "category",
    values_to = "percentage"
  ) %>% 
  filter(category == "any_CAZyme") %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog")),
         population=factor(population, pop.ord)) %>%
  ggplot(aes(x=population, y=percentage, fill=type)) +
  facet_grid(host~body_site, scale="free_y") +
  geom_col(position = "dodge") +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values=site_col) +
  labs(
    title = "CAZymes",
    y = "Percentage of LGT with CAZymes (%)",
    fill = "Category",
    caption = "an LGT can have both ARG and CAZymes"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour="black"),
        strip.text = element_text(face="bold", size=10),
        panel.spacing.y = unit(1, "cm"))

cazy_byPop
save_plot(cazy_byPop, "LGT_function_prop_byPop_CAZY", width = 15, height = 15, unit = "cm",dpi=600, scale=1.2)



#### TOP CARGO
library(patchwork)

types <- c("humanFecal", "dogFecal", "humanOral","dogOral", "humanSkin", "dogSkin")

summary_cargo <- annotated_df %>%
  group_by(type,B_genes, 
           is_amr, is_cazy, 
           Consolidated_Function, what_amr, what_cazy, 
           GO_MOLECULAR_FUNCTION,
           GO_BIOLOGICAL_PROCESS, GO_CELLULAR_COMPONENT) %>%
  summarise(
    ACCESSION = {
      v <- na.omit(ACCESSION)
      v <- unlist(strsplit(v, split = "; "))
      v <- sort(unique(v))
      paste(v, collapse = ";")
    },
    Freq=n(),
    unique_B_genes=n_distinct(index),
    .groups = "drop") %>%
  group_by(type) %>%
  #mutate(percentage = Freq / sum(unique_B_genes), .keep = "all") %>%
  mutate(percentage = Freq / sum(Freq), .keep = "all") %>%
  ungroup() %>%
  arrange(type, desc(percentage))

plot_data <- summary_cargo %>%
  rename(Cargo=B_genes) %>%
  group_by(type) %>%
  slice_max(order_by = percentage, n = 5, with_ties = FALSE) %>%
  mutate(Cargo=gsub("UniRef90_","",Cargo))
  
plot_list <- plot_data %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(type = factor(type, types),
         body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog"))) %>%
  group_by(type) %>%
  group_map(\(sub_df, key) {
    rownames(sub_df) <- NULL
    
    itype <-  unique(key$type)
    
    sub_df <- sub_df %>%
      select(type, Cargo, percentage) %>%
      arrange(desc(percentage)) %>%
      mutate(Cargo = factor(Cargo, unique(Cargo)))
    
    p <- ggplot(sub_df, aes(x = percentage, y = Cargo, fill = type)) +
      geom_col(show.legend = FALSE) +
      theme_minimal() +
      labs(
        title = key$type,
        x = "Percentage in total # Cargo Genes (%)",
        y = "Cargo Genes"
      ) +
      theme(
        plot.title = element_text(face="bold", size=10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5)
      ) +
      scale_fill_manual(values = site_col) +
      scale_y_discrete(limit=rev)
    
    return(p)
  }, .keep = TRUE)
    
combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(
    title = "Most Transferred Genes by Cohort",
    theme = theme(
      title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
  )

combined_plot 
    
save_plot(combined_plot, "cargo_relab", height = 15, width = 15, unit = "cm", scale = 1.2)

## TABLE OUT
out <- plot_data %>%
  mutate(DB = case_when(
    grepl("(^|;)K", ACCESSION) ~ "KO",
    grepl("(^|;)IPR", ACCESSION) ~ "IP",
    grepl("(^|;)PF", ACCESSION) ~ "PF",
    grepl("(^|;)PTHR?", ACCESSION) ~ "PT",
    grepl("(^|;)SSF", ACCESSION) ~ "SF",
    grepl("(^|;)G3D", ACCESSION) ~ "GD",
    TRUE ~ NA_character_
  )) %>%
  mutate(ACCESSION = ifelse(ACCESSION == "", NA_character_, ACCESSION)) %>%
  mutate(percentage = round(percentage, 3)) %>%
  mutate(type = case_when (
    type=="humanFecal" ~ "HF",
    type=="humanOral" ~ "HV",
    type=="humanSkin" ~ "HS",
    type=="dogFecal" ~ "DF",
    type=="dogOral" ~ "DV",
    TRUE ~ "DS"
  )) %>%
  select(type, Cargo, percentage, ACCESSION, DB, Consolidated_Function)

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/top_5_cargo_byType.tsv"

writeLines(header_comment, con = file_path)
write.table(out, file = file_path, quote = FALSE, sep = "\t", row.names = FALSE, append = TRUE)





#### TOP AMR
plot_data <- summary_cargo %>%
  rename(Cargo=B_genes) %>%
  filter(is_amr) %>%
  group_by(type) %>%
  slice_max(order_by = percentage, n = 10, with_ties = FALSE) %>%
  mutate(Cargo=gsub("UniRef90_","",Cargo))

plot_list <- plot_data %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(type = factor(type, types),
         body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog"))) %>%
  group_by(type) %>%
  group_map(\(sub_df, key) {
    rownames(sub_df) <- NULL
    
    itype <-  unique(key$type)
    
    sub_df <- sub_df %>%
      select(type, Cargo, what_amr, percentage) %>%
      mutate(Cargo = paste(Cargo, paste("(",what_amr,")",sep=""), sep="\n")) %>%
      mutate(ppm = round(percentage * 1e4,1)) %>%
      arrange(desc(percentage)) %>%
      mutate(Cargo = factor(Cargo, unique(Cargo)))
    
    p <- ggplot(sub_df, aes(x = ppm, y = Cargo, fill = type)) +
      geom_col(show.legend = FALSE) +
      theme_minimal() +
      labs(
        title = key$type,
        x = "PPM in total # Cargo Genes (x1e-6)",
        y = "Cargo Genes"
      ) +
      theme(
        plot.title = element_text(face="bold", size=10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5)
      ) +
      scale_fill_manual(values = site_col) +
      scale_y_discrete(limit=rev)
    
    return(p)
  }, .keep = TRUE)

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(
    title = "Most Transferred Antimicrobial Resistance Genes by Cohort",
    theme = theme(
      title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
  )

combined_plot 

save_plot(combined_plot, "amr_relab", height = 15, width = 15, unit = "cm", scale = 1.2)

## TABLE OUT
out2 <- plot_data %>%
  mutate(DB = case_when(
    grepl("(^|;)K", ACCESSION) ~ "KO",
    grepl("(^|;)IPR", ACCESSION) ~ "IP",
    grepl("(^|;)PF", ACCESSION) ~ "PF",
    grepl("(^|;)PTHR?", ACCESSION) ~ "PT",
    grepl("(^|;)SSF", ACCESSION) ~ "SF",
    grepl("(^|;)G3D", ACCESSION) ~ "GD",
    TRUE ~ NA_character_
  )) %>%
  mutate(ACCESSION = ifelse(ACCESSION == "", NA_character_, ACCESSION)) %>%
  mutate(ppm = round(percentage*1e4,1)) %>%
  mutate(type = case_when (
    type=="humanFecal" ~ "HF",
    type=="humanOral" ~ "HV",
    type=="humanSkin" ~ "HS",
    type=="dogFecal" ~ "DF",
    type=="dogOral" ~ "DV",
    TRUE ~ "DS"
  )) %>%
  select(type, Cargo, ppm, ACCESSION, DB, Consolidated_Function)

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/top_10_cargo_amr_byType.tsv"

writeLines(header_comment, con = file_path)
write.table(out2, file = file_path, quote = FALSE, sep = "\t", row.names = FALSE, append = TRUE)




#### TOP CAZY
plot_data <- summary_cargo %>%
  rename(Cargo=B_genes) %>%
  filter(is_cazy) %>%
  group_by(type) %>%
  slice_max(order_by = percentage, n = 10, with_ties = FALSE) %>%
  mutate(Cargo=gsub("UniRef90_","",Cargo))

plot_list <- plot_data %>%
  mutate(body_site=gsub("dog|human","",type),
         host=gsub("Fecal|Oral|Skin","",type)) %>%
  mutate(type = factor(type, types),
         body_site=factor(body_site, c("Fecal","Oral","Skin")),
         host=factor(host, c("human","dog"))) %>%
  group_by(type) %>%
  group_map(\(sub_df, key) {
    rownames(sub_df) <- NULL
    
    itype <-  unique(key$type)
    
    sub_df <- sub_df %>%
      select(type, Cargo, what_cazy, percentage) %>%
      mutate(Cargo = paste(Cargo, paste("(",what_cazy,")",sep=""), sep=" ")) %>%
      mutate(ppm = round(percentage * 1e4,1)) %>%
      arrange(desc(percentage)) %>%
      mutate(Cargo = factor(Cargo, unique(Cargo)))
    
    p <- ggplot(sub_df, aes(x = ppm, y = Cargo, fill = type)) +
      geom_col(show.legend = FALSE) +
      theme_minimal() +
      labs(
        title = key$type,
        x = "PPM in total # Cargo Genes (x 1e-6)",
        y = "Cargo Genes"
      ) +
      theme(
        plot.title = element_text(face="bold", size=10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5)
      ) +
      scale_fill_manual(values = site_col) +
      scale_y_discrete(limit=rev)
    
    return(p)
  }, .keep = TRUE)

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(
    title = "Most Transferred CAZymes by Cohort",
    caption = "GT = glycosyltransferases; GH = glycoside hydrolases;\nPL = polysaccharide lyases; CE = carbohydrate esterases;\nAA = auxiliary activities; CBM = carbohydrate-binding modules",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

combined_plot

save_plot(combined_plot, "cazy_relab", height = 15, width = 15, unit = "cm", scale = 1.2)

## TABLE OUT
out3 <- plot_data %>%
  mutate(DB = case_when(
    grepl("(^|;)K", ACCESSION) ~ "KO",
    grepl("(^|;)IPR", ACCESSION) ~ "IP",
    grepl("(^|;)PF", ACCESSION) ~ "PF",
    grepl("(^|;)PTHR?", ACCESSION) ~ "PT",
    grepl("(^|;)SSF", ACCESSION) ~ "SF",
    grepl("(^|;)G3D", ACCESSION) ~ "GD",
    TRUE ~ NA_character_
  )) %>%
  mutate(ACCESSION = ifelse(ACCESSION == "", NA_character_, ACCESSION)) %>%
  mutate(ppm = round(percentage*1e4,1)) %>%
  mutate(type = case_when (
    type=="humanFecal" ~ "HF",
    type=="humanOral" ~ "HV",
    type=="humanSkin" ~ "HS",
    type=="dogFecal" ~ "DF",
    type=="dogOral" ~ "DV",
    TRUE ~ "DS"
  )) %>%
  select(type, Cargo, percentage, ACCESSION, DB, Consolidated_Function)

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/top_10_cargo_cazy_byType.tsv"

writeLines(header_comment, con = file_path)
write.table(out, file = file_path, quote = FALSE, sep = "\t", row.names = FALSE, append = TRUE)



