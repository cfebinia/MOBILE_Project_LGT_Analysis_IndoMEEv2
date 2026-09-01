gc()
rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)

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

annotation_file <- "input_file/annotations/consolidated_gene_annotation_filled_mapping.rds"
annot_details <- "input_file/annotations/gene_annotation_mapping_filled_draft.rds"
master_all_annot <- "input_file/annotations/master_all_annotations_filled_withGO.rds"

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
save_plot <- function(p, name, width, height, unit="in", scale=1) {
  formats <- c("rds", "png", "pdf")
  for(fmt in formats){
    file <- file.path(figpath, paste0(name, ".", fmt))
    if(fmt == "rds") saveRDS(p, file)
    else ggsave(file, plot = p, width = width, height = height, dpi = 150, unit=unit, scale=scale) 
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

indata <- read.delim(file.path(summary_dir,waafle_file))

B_list <- unique(indata$B_genes)
LGT_list <- unique(paste(indata$sandwich, indata$A_gapLen, sep="_"))

LGT.counts <- read.delim(file.path(basedir, LGT_preval))
LGT.mat <- read.delim(file.path(basedir, LGT_matrix), row.names = 1) %>%
  select(-sandwich, -A_gapLen) %>%
  mutate(across(everything(), as.numeric)) %>%
  select(all_of(sort(names(.)))) %>%
  as.matrix()
LGT.mat <- LGT.mat[,!colnames(LGT.mat) %in% outlier_ids]
LGT.mat <- LGT.mat[,colSums(LGT.mat)>0]

match_ids <- data.frame(original_sid=colnames(LGT.mat), SAMPLE=substr(colnames(LGT.mat),1,6))

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


annot_byGenes <- readRDS(file.path(basedir, annotation_file))
annot_bySAMPLE <- readRDS(master_all_annot)
annot_details <- readRDS(file.path(basedir, annot_details))

dim(LGT.mat)
dim(pop.df)
dim(nContigs)

##############
## Analyse  ##
##############
total_counts <- sum(nContigs$n_contig) / 1e6

type_sums <- with(nContigs, tapply(n_contig, factor(substr(SAMPLE, 1, 3), c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH")), sum)) / 1e6
host_sums <- with(nContigs, tapply(n_contig, factor(ifelse(substr(SAMPLE, 1, 2) == "DH", "human", "dog"), c("human", "dog")), sum)) / 1e6

total_nContig <- round(c(type_sums, host_sums, "All" = total_counts), 1)

names(total_nContig) <- c(types, "human","dog", "All")


top_byContig <- LGT.counts %>%
  mutate(human = DHF + DHV + DHP,
         dog = DDF + DDV + DDH) %>%
  rename(All = nContig,
         humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) %>%
  select(all_of(c("sandwichDist", names(total_nContig))))

cols_to_normalize <- names(total_nContig)
sum_per_cols <- colSums(top_byContig[,cols_to_normalize])
top_byContig[, cols_to_normalize] <- sweep(top_byContig[, cols_to_normalize], 2, sum_per_cols, `/`) * 100
top_byContig[, cols_to_normalize] <- sweep(top_byContig[, cols_to_normalize], 2, total_nContig, `/`)

head(top_byContig)

apply(top_byContig[,-1], 2, function(x) sum(x))

# TABLE OUT
write.table(top_byContig, "summary_tables/sandwich_by_Contig.tsv", quote = F,sep = "\t",row.names = F)

#############
# Top Genes #
#############
plot_data <- top_byContig %>%
  pivot_longer(
    cols = all_of(names(total_nContig)),
    names_to = "Group",
    values_to = "pLGT"
  ) %>%
  group_by(Group) %>%
  slice_max(order_by = pLGT, n = 10, with_ties = FALSE) %>%
  arrange(Group, desc(pLGT)) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  mutate(joint_id = paste(Group, sandwichDist , sep = "__")) %>%
  mutate(joint_id = factor(joint_id, levels = unique(joint_id))) %>%
  mutate(Group = factor(Group, names(total_nContig)))

plot_pLGT <- plot_data %>%
  filter(Group !="All") %>%
  mutate(Group=factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  ggplot(aes(x = pLGT, y = joint_id, fill = Group)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Group, scales = "free", ncol=2) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
  theme_minimal() +
  labs(
    title = "Top 10 Most Abundant Sandwich Gene Synteny",
    x = "Relative Abundance (%) / Total Contigs",
    y = "Gene Synteny"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 7.5)
  ) +
  scale_fill_manual(values = site_col)

plot_pLGT

save_plot(plot_pLGT, "sandwich_relab", height = 8,width = 15, unit = "cm", scale = 1.8)


# TABLE OUT
top10_sandwich <- plot_data
write.table(top10_sandwich, "summary_tables/top_10_sandwich_byContig.tsv", quote = F,sep = "\t",row.names = F)

annotated_top10 <- top10_sandwich %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  left_join(annot_byGenes, by = "ANNOTATIONS.UNIREF90", keep = FALSE) #

out <- top_byContig %>%
  right_join(annotated_top10 %>% select(-pLGT, -joint_id),
            by="sandwichDist", keep=FALSE) %>%
  filter(!Group %in% c("human","dog")) %>%
  mutate(Group=factor(Group,c("All",types))) %>%
  arrange(Group, Rank, index) %>%
  select(-index) %>%
  rename(sandwich = sandwichDist , Group_Rank = Rank) %>%
  select(Group, Group_Rank, ANNOTATIONS.UNIREF90, All,
         humanFecal, humanOral, humanSkin,
         dogFecal, dogOral, dogSkin, 
         ACCESSION, Consolidated_Function, 
         is_amr, is_cazy, sandwich) %>%
  mutate(across(all_of(c(types, "All")), \(x) ifelse(x==0,NA,round(x*100, 2)))) %>%
  mutate(across(all_of(c(types, "All")), \(x) as.character(x))) %>%
  mutate(across(all_of(c(types, "All")), \(x) ifelse(is.na(x), "--", x)))

write.table(out, "summary_tables/top_10_sandwich_genes_by_Contig_full.tsv",
            quote = F, sep="\t", row.names = F)
  
short_out <-  out %>%
  select(Group, Group_Rank, ANNOTATIONS.UNIREF90,ACCESSION,Consolidated_Function, All) %>%
  rowwise() %>%
  mutate(ACCESSION = {
    v <- unlist(strsplit(ACCESSION, split = "; ")) ;
    n <- length(v) ;
    ifelse(n > 3, "Various", ACCESSION)
  },
  Consolidated_Function = {
    v <- unlist(strsplit(ACCESSION, split = "; ")) ;
    n <- length(v) ;
    ifelse(n > 3, "See Supplementary Table", Consolidated_Function)
  }) %>%
  ungroup() %>%
  rename(total_pLGT=All)


cat("## If Group_Rank is tied, the gene came with the same sandwich\n",
    file = "summary_tables/top_10_sandwich_genes_by_Contig_short.tsv", append = F)
write.table(short_out, "summary_tables/top_10_sandwich_genes_by_Contig_short.tsv",
            quote = F, sep="\t", row.names = F, append = T)

###########
# Top AMR #
###########
amr_data <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  inner_join(
    annot_byGenes %>% 
      filter(is_amr) %>% 
      select(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_amr) %>%
      distinct(),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  pivot_longer(
    cols = all_of(names(total_nContig)),
    names_to = "Group",
    values_to = "pLGT"
  ) %>%
  mutate(Group = factor(Group, c(types, "human", "dog", "All"))) %>%
  group_by(Group) %>%
  group_map(\(sub_df, key) {
    sub_df <- sub_df %>%
      distinct(index, ANNOTATIONS.UNIREF90, .keep_all = TRUE) %>%
      group_by(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_amr) %>%
      summarise(pLGT = sum(pLGT, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(pLGT)) %>%
      mutate(pLGT = as.numeric(pLGT)) %>%
      filter(pLGT != 0)
    
    if (nrow(sub_df) == 0) {
      return(NULL)
    }
    
    sub_df <- sub_df %>%
      mutate(ANNOTATIONS.UNIREF90 = factor(ANNOTATIONS.UNIREF90, unique(ANNOTATIONS.UNIREF90))) %>%
      mutate(Group_Rank = row_number()) %>%
      slice_max(order_by = pLGT, n = 10, with_ties = FALSE)
    
    return(data.frame(Group = key$Group, sub_df))
  }, .keep = TRUE) %>% 
  bind_rows()


library(patchwork)

plot_list <- amr_data %>%
  filter(Group != "All") %>%
  mutate(Group=factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  # mutate(Group = factor(Group, c(types, "human", "dog", "All"))) %>%
  group_by(Group) %>%
  group_map(\(sub_df, key) {
    rownames(sub_df) <- NULL
    sub_df <- sub_df %>%
      group_by(Group, ANNOTATIONS.UNIREF90) %>%
      summarise(pLGT=sum(pLGT), .groups = "drop") %>%
      arrange(desc(pLGT)) %>%
      mutate(pLGT = as.numeric(pLGT)) %>%
      filter(pLGT != 0) %>%
      mutate(ANNOTATIONS.UNIREF90 = factor(ANNOTATIONS.UNIREF90, unique(ANNOTATIONS.UNIREF90)))
    
    p <- ggplot(sub_df, aes(x = pLGT, y = ANNOTATIONS.UNIREF90, fill = Group)) +
      geom_col(show.legend = FALSE) +
      scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
      theme_minimal() +
      labs(
        title = key$Group,
        x = "Relative Abundance (%) / Total Contigs",
        y = "sandwich Genes"
      ) +
      theme(
        plot.title = element_text(face="bold", size=10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5)
      ) +
      scale_fill_manual(values = site_col)
    
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

save_plot(combined_plot, "AMR_relab", height = 8,width = 15, unit = "cm", scale = 1.8)


## by AMR class
amrclass_data <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  inner_join(
    annot_byGenes %>% filter(is_amr),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  select(-is_cazy, -what_cazy) %>%
  group_by(sandwichDist ) %>%
  summarise(
    across(all_of(c(types, "All", "human", "dog")), \(x) as.numeric(unique(x))),
    what_amr = paste(unique(what_amr[what_amr != ""]), collapse = "+"),
    .groups = "drop"
  ) %>%
  group_by(what_amr) %>%
  summarise(
    across(all_of(c(types, "All", "human", "dog")), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(names(total_nContig)),
    names_to = "Group",
    values_to = "pLGT"
  ) %>%
  arrange(Group, desc(pLGT)) %>%
  group_by(Group) %>%
  mutate(Rank = row_number()) %>%
  mutate(Rank = ifelse(pLGT == 0, NA_integer_, Rank)) %>%
  filter(!is.na(Rank)) %>%
  ungroup() %>%
  distinct(what_amr, Group, Rank, .keep_all = TRUE) %>%
  mutate(joint_id = paste(Group, what_amr, sep = "__")) %>%
  mutate(joint_id = factor(joint_id, levels = unique(joint_id))) %>%
  mutate(Group = factor(Group, names(total_nContig)))

plot_pLGT <- amrclass_data %>%
  filter(Group != "All") %>%
  mutate(Group=factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  ggplot(aes(x = pLGT, y = joint_id, fill = Group)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Group, scales = "free", ncol=2) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
  theme_minimal() +
  labs(
    title = "Top 10 Most Transferred ARG Class",
    x = "Relative Abundance (%) / Total Contigs",
    y = "sandwich Genes"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 7.5)
  ) +
  scale_fill_manual(values = site_col)

plot_pLGT

save_plot(plot_pLGT, "AMRclass_relab", height = 8,width = 15, unit = "cm", scale = 1.8)


# TABLE OUT
out2 <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  filter(ANNOTATIONS.UNIREF90 %in% annot_byGenes$ANNOTATIONS.UNIREF90[annot_byGenes$is_amr]) %>%
  right_join(
    amr_data %>% select(ANNOTATIONS.UNIREF90, Group, Group_Rank),
    by = "ANNOTATIONS.UNIREF90",
    relationship = "many-to-many"
  ) %>%
  arrange(Group, desc(Group_Rank)) %>%
  group_by(Group, ANNOTATIONS.UNIREF90) %>%
  mutate(n_sandwich_variant=n_distinct(sandwichDist )) %>%
  ungroup() %>%
  inner_join(
    annot_byGenes %>% 
      filter(is_amr) %>% 
      select(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_amr) %>%
      distinct(),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  filter(!Group %in% c("human","dog")) %>%
  mutate(Group=factor(Group,c("All",types))) %>%
  arrange(Group, Group_Rank, index) %>%
  select(-index) %>%
  rename(sandwich = sandwichDist ) %>%
  select(Group, Group_Rank, ANNOTATIONS.UNIREF90, n_sandwich_variant,
         All, humanFecal, humanOral, humanSkin,
         dogFecal, dogOral, dogSkin, 
         ACCESSION, Consolidated_Function, 
         what_amr, sandwich) 


long_out <- out2 %>%
  distinct() %>%
  select(-sandwich) %>%
  group_by(Group, Group_Rank, n_sandwich_variant, 
           ANNOTATIONS.UNIREF90, ACCESSION, 
           Consolidated_Function, what_amr) %>%
  summarise(
    across(all_of(c(types, "All")), \(x) {
      val <- sum(as.numeric(x), na.rm = TRUE) * 100
      if (val == 0) NA else round(val, 2)
    }),
    .groups = "drop"
  ) %>%
  mutate(across(all_of(c(types, "All")), \(x) as.character(x))) %>%
  mutate(across(all_of(c(types, "All")), \(x) ifelse(is.na(x), "--", x)))

write.table(long_out, "summary_tables/top_10_ARG_by_Contig_short.tsv",
            quote = F, sep="\t", row.names = F)

short_out <- long_out %>%
  select(-all_of(c(types))) %>%
  rename(total_pLGT=All)

cat("## total_pLGT is sum of pLGT across all samples and relevant sandwich variant\n",
    file = "summary_tables/top_10_sandwich_genes_by_Contig_short.tsv", append = F)
write.table(short_out, "summary_tables/top_10_ARG_by_Contig_short.tsv",
            quote = F, sep="\t", row.names = F)


############
# Top CAZY #
############
cazy_data <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  inner_join(
    annot_byGenes %>% 
      filter(is_cazy) %>% 
      select(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_cazy) %>%
      distinct(),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  pivot_longer(
    cols = all_of(names(total_nContig)),
    names_to = "Group",
    values_to = "pLGT"
  ) %>%
  mutate(Group = factor(Group, c(types, "human", "dog", "All"))) %>%
  group_by(Group) %>%
  group_map(\(sub_df, key) {
    sub_df <- sub_df %>%
      distinct(index, ANNOTATIONS.UNIREF90, .keep_all = TRUE) %>%
      group_by(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_cazy) %>%
      summarise(pLGT = sum(pLGT, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(pLGT)) %>%
      mutate(pLGT = as.numeric(pLGT)) %>%
      filter(pLGT != 0)
    
    if (nrow(sub_df) == 0) {
      return(NULL)
    }
    
    sub_df <- sub_df %>%
      mutate(ANNOTATIONS.UNIREF90 = factor(ANNOTATIONS.UNIREF90, unique(ANNOTATIONS.UNIREF90))) %>%
      mutate(Group_Rank = row_number()) %>%
      slice_max(order_by = pLGT, n = 10, with_ties = FALSE)
    
    return(data.frame(Group = key$Group, sub_df))
  }, .keep = TRUE) %>% 
  bind_rows()


library(patchwork)

plot_list <- cazy_data %>%
  filter(Group != "All") %>%
  mutate(Group=factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  group_by(Group) %>%
  group_map(\(sub_df, key) {
    rownames(sub_df) <- NULL
    sub_df <- sub_df %>%
      group_by(Group, ANNOTATIONS.UNIREF90) %>%
      summarise(pLGT=sum(pLGT), .groups = "drop") %>%
      arrange(desc(pLGT)) %>%
      mutate(pLGT = as.numeric(pLGT)) %>%
      filter(pLGT != 0) %>%
      mutate(ANNOTATIONS.UNIREF90 = factor(ANNOTATIONS.UNIREF90, unique(ANNOTATIONS.UNIREF90)))
    
    p <- ggplot(sub_df, aes(x = pLGT, y = ANNOTATIONS.UNIREF90, fill = Group)) +
      geom_col(show.legend = FALSE) +
      scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
      theme_minimal() +
      labs(
        title = key$Group,
        x = "Relative Abundance (%) / Total Contigs",
        y = "sandwich Genes"
      ) +
      theme(
        plot.title = element_text(face="bold", size=10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5)
      ) +
      scale_fill_manual(values = site_col)
    
    return(p)
  }, .keep = TRUE)

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(
    title = "Most Transferred CaZy Genes by Cohort",
    theme = theme(
      title = element_text(face = "bold", size = 11, hjust = 0.5)
    )
  )

combined_plot

save_plot(combined_plot, "CAZY_relab", height = 8,width = 15, unit = "cm", scale = 1.8)



## by CAZY class
cazyclass_data <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  inner_join(
    annot_byGenes %>% filter(is_cazy),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  select(-is_amr, -what_amr) %>%
  mutate(what_cazy = gsub("; ", "|", what_cazy)) %>%
  rowwise() %>%
  mutate(cazy_class = {
    v <- unlist(strsplit(what_cazy, split = "\\|"))
    v <- unique(v)
    v <- gsub("[[:digit:]]", "", v)
    paste(v, collapse = "|")
  }) %>%
  ungroup() %>%
  group_by(sandwichDist ) %>%
  summarise(
    across(all_of(c(types, "All", "human", "dog")), \(x) as.numeric(unique(x))),
    what_cazy = paste(unique(what_cazy[what_cazy != ""]), collapse = "+"),
    .groups = "drop"
  ) %>%
  mutate(what_cazy=gsub("; ","+", what_cazy)) %>%
  group_by(what_cazy) %>%
  summarise(
    across(all_of(c(types, "All", "human", "dog")), \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(names(total_nContig)),
    names_to = "Group",
    values_to = "pLGT"
  ) %>%
  arrange(Group, desc(pLGT)) %>%
  group_by(Group) %>%
  mutate(Rank = row_number()) %>%
  mutate(Rank = ifelse(pLGT == 0, NA_integer_, Rank)) %>%
  filter(!is.na(Rank)) %>%
  filter(Rank <10) %>%
  ungroup() %>%
  distinct(what_cazy, Group, Rank, .keep_all = TRUE) %>%
  mutate(joint_id = paste(Group, what_cazy, sep = "__")) %>%
  mutate(joint_id = factor(joint_id, levels = unique(joint_id))) %>%
  mutate(Group = factor(Group, names(total_nContig)))

plot_pLGT <- cazyclass_data %>%
  filter(Group != "All") %>%
  mutate(Group=factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  ggplot(aes(x = pLGT, y = joint_id, fill = Group)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Group, scales = "free", ncol=2 ) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
  theme_minimal() +
  labs(
    x = "Relative Abundance (%) / Total Contigs",
    y = "Family",
    caption = "GT = glycosyltransferases; GH = glucoside hydrolases; 
    PL = polysaccharide lyases; = carbohydrate esterases; 
    AA = auzillary activities; CBM = carbohydrate-binding modiles"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 7.5)
  ) +
  scale_fill_manual(values = site_col)

plot_pLGT

save_plot(plot_pLGT, "cazyclass_relab", height = 18,width = 15, unit = "cm", scale = 1.8)


# TABLE OUT
out3 <- top_byContig %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(sandwichDist , pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  filter(ANNOTATIONS.UNIREF90 %in% annot_byGenes$ANNOTATIONS.UNIREF90[annot_byGenes$is_cazy]) %>%
  right_join(
    cazy_data %>% select(ANNOTATIONS.UNIREF90, Group, Group_Rank),
    by = "ANNOTATIONS.UNIREF90",
    relationship = "many-to-many"
  ) %>%
  arrange(Group, desc(Group_Rank)) %>%
  group_by(Group, ANNOTATIONS.UNIREF90) %>%
  mutate(n_sandwich_variant=n_distinct(sandwichDist )) %>%
  ungroup() %>%
  inner_join(
    annot_byGenes %>% 
      filter(is_cazy) %>% 
      select(ANNOTATIONS.UNIREF90, ACCESSION, Consolidated_Function, what_cazy) %>%
      distinct(),
    by = "ANNOTATIONS.UNIREF90"
  ) %>%
  filter(!Group %in% c("human","dog")) %>%
  mutate(Group=factor(Group,c("All",types))) %>%
  arrange(Group, Group_Rank, index) %>%
  select(-index) %>%
  rename(sandwich = sandwichDist ) %>%
  mutate(what_cazy = gsub("; ", "|", what_cazy)) %>%
  select(Group, Group_Rank, ANNOTATIONS.UNIREF90, n_sandwich_variant,
         All, humanFecal, humanOral, humanSkin,
         dogFecal, dogOral, dogSkin, 
         ACCESSION, Consolidated_Function, 
         what_cazy, sandwich) 

long_out <- out3 %>%
  distinct() %>%
  select(-sandwich) %>%
  group_by(Group, Group_Rank, n_sandwich_variant, 
           ANNOTATIONS.UNIREF90, ACCESSION, 
           Consolidated_Function, what_cazy) %>%
  summarise(
    across(all_of(c(types, "All")), \(x) {
      val <- sum(as.numeric(x), na.rm = TRUE) * 100
      if (val == 0) NA else round(val, 2)
    }),
    .groups = "drop"
  ) %>%
  mutate(across(all_of(c(types, "All")), \(x) as.character(x))) %>%
  mutate(across(all_of(c(types, "All")), \(x) ifelse(is.na(x), "--", x)))

write.table(long_out, "summary_tables/top_10_CAZY_by_Contig_short.tsv",
            quote = F, sep="\t", row.names = F)

short_out <- long_out %>%
  select(-all_of(c(types))) %>%
  rename(total_pLGT=All)

cat("## total_pLGT is sum of pLGT across all samples and relevant sandwich variant\n",
    file = "summary_tables/top_10_sandwich_genes_by_Contig_short.tsv", append = F)
write.table(short_out, "summary_tables/top_10_CAZY_by_Contig_short.tsv",
            quote = F, sep="\t", row.names = F)


