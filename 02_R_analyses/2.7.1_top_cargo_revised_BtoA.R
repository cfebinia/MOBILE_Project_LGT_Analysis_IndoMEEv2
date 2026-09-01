gc()
rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)

# basedir="Imaging_Lab_PC1/WAAFLE_Extra"
basedir="Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

summary_dir=file.path(basedir, "summary_tables")
figpath=file.path(basedir, "figures")
out_dir=file.path(basedir, "summary_tables")

waafle_file <- "all_samples_internal_merged_indexed_filtered.tsv"
contigs_stats <- "input_file/contigs/summary_contig_counts_bySample.tsv"
contigs_N50 <- "input_file/contigs/summary_n50_bySample.tsv"
contigs_stats_byType <- "input_file/contigs/summary_contig_stats_byType.tsv"

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

contigs_stats <- "input_file/contigs/summary_contig_counts_bySample.tsv"
contigs_N50 <- "input_file/contigs/summary_n50_bySample.tsv"
contigs_stats_byType <- "input_file/contigs/summary_contig_stats_byType.tsv"

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
  "humanOral"="#ef9f00", 
  "dogFecal"="#2850c8", 
  "dogSkin"="#780247", 
  "dogOral"="#9A5324", 
  "human" = "grey70",
  "dog" = "grey40",
  "All" = "grey10"
)

################
## Load Data  ##
################
outlier_ids <- read.delim(file.path(basedir,outliers))$sampleid
final_samples <- readRDS("input_file/final_sample_list.rds")

indata <- read.delim(file.path(summary_dir,waafle_file)) %>%
  rowwise() %>%
  mutate(n_Bgenes = length(unlist(strsplit(B_genes, split = "\\|")))) %>%
  ungroup() %>%
  filter(n_Bgenes <= 10 & nchar(SYNTENY) <= 35) %>%
  filter(original_sid %in% final_samples) 

B_list <- indata %>%
  filter(DIRECTION=="B>A") %>%
  pull(B_genes) %>%
  unique()

cargo.counts <- read.delim(file.path(basedir, cargo_preval)) %>%
  filter(B_genes %in% B_list) %>%
  rowwise() %>%
  mutate(Blen = length(unlist(strsplit(B_genes, split = "\\|")))) %>%
  ungroup()

cargo.mat <- read.delim(file.path(basedir, cargo_matrix), row.names = 1) %>%
  mutate(across(everything(), as.numeric)) %>%
  select(all_of(sort(names(.)))) %>%
  as.matrix()
cargo.mat <- cargo.mat[,colnames(cargo.mat) %in% final_samples]
cargo.mat <- cargo.mat[rownames(cargo.mat) %in% B_list,]
cargo.mat <- cargo.mat[,colSums(cargo.mat)>0]

LGT.counts <- read.delim(file.path(basedir, LGT_preval))
LGT.mat <- read.delim(file.path(basedir, LGT_matrix), row.names = 1) %>%
  select(-sandwich, -A_gapLen) %>%
  mutate(across(everything(), as.numeric)) %>%
  select(all_of(sort(names(.)))) %>%
  as.matrix()
LGT.mat <- LGT.mat[,colnames(LGT.mat) %in% final_samples]
LGT.mat <- LGT.mat[,colSums(LGT.mat)>0]

match_ids <- data.frame(original_sid=colnames(cargo.mat), SAMPLE=substr(colnames(cargo.mat),1,6))

pop.df <- read.delim(file.path(basedir, "input_file/Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(original_sid = gsub("_POOLED|_TP","", seqID)) %>%
  select(original_sid, population, individual) %>%
  distinct() %>%
  rename(SAMPLE=original_sid) %>%
  mutate(short_sample = substr(SAMPLE, 1,6)) %>%
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

nContigs <- read.delim(contigs_stats, skip=1) %>%
  left_join(read.delim(contigs_N50, skip=1), by="sampleid") %>%
  filter(sampleid %in% final_samples) %>%
  mutate(short_sample = substr(sampleid, 1,6)) %>%
  rename(SAMPLE=sampleid) %>%
  select(SAMPLE, short_sample, n_contig=count)

annot_byGenes <- readRDS(file.path(basedir, annotation_file))
annot_bySAMPLE <- readRDS(master_all_annot)
annot_details <- readRDS(file.path(basedir, annot_details))

dim(cargo.mat)
dim(LGT.mat) # n sample is lower than cargo.mat, because of no singletons
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

total_genes_bySample <- indata %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         jointContig=paste(original_sid, CONTIG_NAME)) %>%
  select(original_sid, type, B_genes, sandwichDist, jointContig) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(original_sid, type) %>%
  summarise(n_Bgenes=n_distinct(B_genes),
            n_Contig=n_distinct(jointContig),
            n_Sandwich=n_distinct(sandwichDist), 
            .groups = "drop")

total_genes_byType <- indata %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         jointContig=paste(original_sid, CONTIG_NAME),
         type = case_when(
           type == "DHF" ~ "humanFecal",
           type == "DHV" ~ "humanOral",
           type == "DHP" ~ "humanSkin",
           type == "DDF" ~ "dogFecal",
           type == "DDV" ~ "dogOral",
           TRUE ~ "dogSkin"
         )) %>%
  select(original_sid, type, B_genes, sandwichDist, jointContig) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(type) %>%
  summarise(n_Bgenes=n_distinct(B_genes),
            n_Contig=n_distinct(jointContig),
            n_Sandwich=n_distinct(sandwichDist), 
            .groups = "drop") %>%
  rename(Group=type)

total_genes_byHost <- indata %>%
  filter(original_sid %in% final_samples) %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         jointContig=paste(original_sid, CONTIG_NAME),
         host=case_when(
           type %in% c("DHF","DHV","DHP") ~ "human",
           TRUE ~ "dog")) %>%
  select(original_sid, host, B_genes, sandwichDist, jointContig) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(host) %>%
  summarise(n_Bgenes=n_distinct(B_genes),
            n_Contig=n_distinct(jointContig),
            n_Sandwich=n_distinct(sandwichDist), 
            .groups = "drop") %>%
  rename(Group=host)

total_genes_all <- indata %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         jointContig=paste(original_sid, CONTIG_NAME),
         host=case_when(
           type %in% c("DHF","DHV","DHP") ~ "human",
           TRUE ~ "dog")) %>%
  select(original_sid, host, B_genes, sandwichDist, jointContig) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  summarise(Group = "All",
            n_Bgenes=n_distinct(B_genes),
            n_Contig=n_distinct(jointContig),
            n_Sandwich=n_distinct(sandwichDist), 
            .groups = "drop")

top_byContig <- cargo.counts %>%
  mutate(human = DHF + DHV + DHP,
         dog = DDF + DDV + DDH) %>%
  rename(All = nContig,
         humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) %>%
  select(all_of(c("B_genes", names(total_nContig))))

cols_to_normalize <- names(total_nContig)
sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Contig, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[cols_to_normalize]

top_byContig[, cols_to_normalize] <- sweep(top_byContig[, cols_to_normalize], 2, sum_per_cols, `/`) *100

head(top_byContig)

# TABLE OUT
write.table(top_byContig, "summary_tables/Cargo_synteny_abund_rev.tsv", quote = F,sep = "\t",row.names = F)

#############
# Top Genes #
#############
sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Sandwich, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[cols_to_normalize]

top_bySandwich <- indata %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         jointContig=paste(original_sid, CONTIG_NAME),
         host=case_when(
           type %in% c("DHF","DHV","DHP") ~ "human",
           TRUE ~ "dog")) %>%
  select(original_sid, type, B_genes, sandwichDist) %>%
  distinct() %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  group_by(B_genes, type) %>%
  summarise(count = n_distinct(sandwichDist),
            .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
  mutate(human = DHF + DHV + DHP,
         dog = DDF + DDV + DDH,
         All = human + dog) %>%
  rename(humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) 
  
top_bySandwich[, cols_to_normalize] <- sweep(top_bySandwich[, cols_to_normalize], 2, sum_per_cols, `/`) *100

top_bySandwich <- top_bySandwich %>%
  arrange(desc(All))

plot_data <- top_bySandwich %>%
  mutate(B_genes=gsub("UniRef90_","",B_genes)) %>%
  pivot_longer(
    cols = all_of(cols_to_normalize),
    names_to = "Group",
    values_to = "pct_sandwich"
  ) %>%
  group_by(Group) %>%
  slice_max(order_by = pct_sandwich, n = 10, with_ties = FALSE) %>%
  arrange(Group, desc(pct_sandwich)) %>%
  mutate(Rank = row_number()) %>%
  ungroup() %>%
  mutate(joint_id = paste(Group, B_genes, sep = "__")) %>%
  mutate(joint_id = factor(joint_id, levels = unique(joint_id))) %>%
  mutate(Group = factor(Group, cols_to_normalize))

plot_pLGT <- plot_data %>%
  filter(Group %in% types) %>%
  mutate(Group = factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog"))) %>%
  ggplot(aes(x = pct_sandwich, y = joint_id, fill = Group)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Group, scales = "free", ncol=2) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x), limits = rev) +
  theme_minimal() +
  labs(
    title = "Top 5 Most Prevalent Transferred Genes in LGT events",
    x = "Proportion in Unique LGT Events (%)",
    y = "Cargo Genes"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 7.5),
    plot.title = element_text(size = 10, face="bold"),
    axis.title.y=element_blank(),
    axis.title.x=element_text(size=9)
  ) +
  scale_fill_manual(values = site_col)

plot_pLGT

save_plot(plot_pLGT, "cargo_relab", height = 16,width = 15, unit = "cm", scale = 1)


# TABLE OUT
top10_cargo <- plot_data

header_comment <- "# Note: Proportion  (%) of LGT identified with a specific gene."
file_path <- "summary_tables/top_10_cargo_bySandwich.tsv"

writeLines(header_comment, con = file_path)
write.table(top10_cargo, file_path, 
            quote = F,sep = "\t",row.names = F)

annotated_top10 <- top10_cargo %>%
  mutate(index = row_number()) %>%
  mutate(ANNOTATIONS.UNIREF90 = str_split(B_genes, pattern = "\\|")) %>%
  unnest(ANNOTATIONS.UNIREF90) %>%
  left_join(annot_byGenes %>% mutate(ANNOTATIONS.UNIREF90 = gsub("UniRef90_", "", ANNOTATIONS.UNIREF90)), by = "ANNOTATIONS.UNIREF90", keep = FALSE) #

out <- top_bySandwich %>%
mutate(B_genes = gsub("UniRef90_", "", B_genes)) %>% 
  right_join(annotated_top10 %>% select(-pct_sandwich, -joint_id),
             by = "B_genes", keep = FALSE) %>%
  filter(!Group %in% c("human", "dog")) %>%
  mutate(Group = factor(Group, levels = c("All", types))) %>%
  arrange(Group, Rank, index) %>%
  select(-index) %>%
  rename(Cargo = B_genes, Group_Rank = Rank) %>%
  select(Group, Group_Rank, ANNOTATIONS.UNIREF90, All,
         humanFecal, humanOral, humanSkin,
         dogFecal, dogOral, dogSkin, 
         ACCESSION, Consolidated_Function, 
         is_amr, is_cazy) %>%
  mutate(across(all_of(c(types, "All")), \(x) as.character(round(x, 3)))) %>%
  mutate(across(all_of(c("ACCESSION", "Consolidated_Function")), \(x) ifelse(is.na(x), "--", x)))

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/top_10_cargo_genes_by_Sandwich_full.tsv"

writeLines(header_comment, con = file_path)
write.table(out, file_path,
            quote = F, sep="\t", row.names = F)
  
short_out <-  out %>%
  filter(Group %in% types) %>%
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
  rename(pct_sandwich=All) %>%
  mutate(ANNOTATIONS.UNIREF90=gsub("UniRef90_","",ANNOTATIONS.UNIREF90)) %>%
  rename(UniRef90=ANNOTATIONS.UNIREF90) %>%
  mutate(DB = case_when(
    grepl("K",ACCESSION) ~ "KO",
    grepl("IP",ACCESSION) ~ "IP",
    grepl("PT",ACCESSION) ~ "PT",
    grepl("PF",ACCESSION) ~ "PF",
    grepl("SS",ACCESSION) ~ "SF",
    grepl("G3D",ACCESSION) ~ "GD",
    grepl("Var",ACCESSION) ~ "Various",
    TRUE ~ "--"
  ))  %>%
  select(Group, UniRef90, pct_sandwich, ACCESSION, DB,Consolidated_Function)


header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/top_10_cargo_genes_by_Sandwich_short.tsv"

writeLines(header_comment, con = file_path)
write.table(short_out, file_path,
            quote = F, sep="\t", row.names = F, append = T)


##############
# ANNOTATION #
##############
master_all_annot <- "input_file/annotations/master_all_annotations_filled_withGO.rds"
annot_bySAMPLE <- readRDS(master_all_annot)

df <- indata %>%
  filter(DIRECTION == "B>A") %>%
  mutate(sandwichDist=paste(sandwich, A_gapLen, sep="_"),
         host=case_when(
           type %in% c("DHF","DHV","DHP") ~ "human",
           TRUE ~ "dog")) %>%
  select(sandwichDist, CONTIG_NAME, B_genes, original_sid, type, host) %>% 
  mutate(index = row_number()) %>%
  separate_rows(B_genes, sep = "\\|") %>%
  group_by(index) %>%
  mutate(g_index = row_number()) %>%
  ungroup() %>%
  rename(SAMPLE=original_sid) %>%
  mutate(pair = paste(SAMPLE, CONTIG_NAME, B_genes, sep = "+"))

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

annotated_df <- annotated_df %>%
  filter(SAMPLE %in% final_samples) %>%
  mutate(short_sample=substr(SAMPLE, 1,6)) %>%
  left_join(pop.df %>% 
              select(short_sample, individual, population) %>% 
              distinct(),
            by="short_sample")
  

annotated_df %>% 
  group_by(type) %>%
  summarise(n=n_distinct(SAMPLE),
            n_indv = n_distinct(individual)
  )

annotated_df %>% 
  group_by(type, population) %>%
  summarise(n=n_distinct(SAMPLE),
            n_indv = n_distinct(individual)
  )

#View(annotated_df[is.na(annotated_df$population),])

unique(annotated_df$short_sample[is.na(annotated_df$population)])

table(annotated_df$short_sample %in% final_samples)

##############
# GO PROCESS #
##############
goProcess_byType <- annotated_df %>%
  mutate(GO_BIOLOGICAL_PROCESS = ifelse(is.na(GO_BIOLOGICAL_PROCESS), "No GO Entry", GO_BIOLOGICAL_PROCESS)) %>%
  group_by(type, GO_BIOLOGICAL_PROCESS) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

goProcess_byHost<- annotated_df %>%
  mutate(GO_BIOLOGICAL_PROCESS = ifelse(is.na(GO_BIOLOGICAL_PROCESS), "No GO Entry", GO_BIOLOGICAL_PROCESS)) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(host, GO_BIOLOGICAL_PROCESS) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

goProcess_All <- annotated_df %>%
  mutate(GO_BIOLOGICAL_PROCESS = ifelse(is.na(GO_BIOLOGICAL_PROCESS), "No GO Entry", GO_BIOLOGICAL_PROCESS)) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(GO_BIOLOGICAL_PROCESS) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

summary_GOProcess <- goProcess_All %>% select(GO_BIOLOGICAL_PROCESS, n_Sandwich) %>% rename(All=n_Sandwich) %>%
  left_join(goProcess_byType %>% select(GO_BIOLOGICAL_PROCESS, type, n_Sandwich) %>%
              pivot_wider(names_from = type, values_from = n_Sandwich, values_fill = 0),
            by="GO_BIOLOGICAL_PROCESS") %>%
  left_join(goProcess_byHost %>% select(GO_BIOLOGICAL_PROCESS, host, n_Sandwich) %>%
              pivot_wider(names_from = host, values_from = n_Sandwich, values_fill = 0),
            by="GO_BIOLOGICAL_PROCESS") %>%
  arrange(desc(All)) %>%
  rename(humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) 


sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Sandwich, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[cols_to_normalize]

summary_GOProcess[, cols_to_normalize] <- sweep(summary_GOProcess[, cols_to_normalize], 2, sum_per_cols, `/`) *100


plot_data <- summary_GOProcess %>%
  rename(Process = GO_BIOLOGICAL_PROCESS) %>%
  pivot_longer(cols = all_of(cols_to_normalize), names_to = "Group", values_to = "pct_sandwich") %>%
  mutate(Group = factor(Group, c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin","human","dog", "All"))) %>%
  group_by(Group) %>%
  slice_max(order_by = pct_sandwich, n = 11, with_ties = FALSE) %>%
  ungroup()

group_level <- c(types)

plot_list <- lapply(group_level, function(x){
  sub_df <- plot_data %>%
    filter(Group != "All") %>%
    mutate(Process = gsub("GO:00","GO:",Process)) %>%
    rename(percentage=pct_sandwich) %>%
    filter(Group==x) %>%
    select(Group, Process, percentage) %>%
    arrange(desc(percentage)) %>%
    mutate(percentage = round(percentage, 3)) %>%
    filter(percentage > 0) %>%
    mutate(Process = droplevels(factor(Process, levels = rev(unique(Process)))))
  
  p <- ggplot(sub_df, aes(x = percentage, y = Process, fill = Group)) +
    geom_col(show.legend = FALSE) +
    theme_minimal() +
    labs(title = x, x="Proportion in Unique LGT (%)") +
    theme(
      legend.position = "inside",
      plot.title = element_text(face = "bold", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 7.5),
      axis.title.y = element_blank()
    ) +
    scale_fill_manual(values = site_col) +
    geom_text(
      data = subset(sub_df, Process == "No GO Entry"),
      aes(x = max(percentage)+1, label = "*"),
      color = "black",
      size = 5,
      vjust = 0.8
    ) +
    scale_x_continuous(limits = c(0,5))
  
  vals <- ceiling(sort(unique(round(sub_df$percentage,1))))
  
  if (length(vals) >= 2) {
    max_val <- round(max(sub_df$percentage))
    second_max <- round(vals[length(vals) - 1]*1.1)
    
    y_names <- levels(sub_df$Process)
    
    p1 <- sub_df %>%
      mutate(percentage = ifelse(percentage>second_max,second_max,percentage)) %>%
      ggplot(aes(x = percentage, y = Process, fill = Group)) +
      geom_col(show.legend = FALSE) +
      theme_minimal(base_size=10) +
      labs(title = x) +
      theme(
        legend.position = "inside",
        plot.title = element_text(face = "bold", size = 10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 7.5),
        axis.title = element_blank()
      ) +
      scale_fill_manual(values = site_col) +
      scale_x_continuous(limits = c(0,second_max), expand = c(0,0))+
      coord_cartesian(xlim = c(0, second_max)) +
      scale_y_discrete(drop = FALSE)
    
    p2 <- sub_df %>%
      mutate(percentage = ifelse(percentage > second_max, percentage, 0)) %>%
      ggplot(aes(x = percentage, y = Process, fill = Group)) +
      geom_col(show.legend = FALSE) +
      theme_minimal(base_size=10) +
      labs(title = "") +
      theme(
        legend.position = "inside",
        plot.title = element_text(face = "bold", size = 10),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_blank(),
        axis.title = element_blank()
      ) +
      scale_fill_manual(values = site_col) +
      geom_text(
        data = subset(sub_df, Process == "No GO Entry"),
        aes(x = max(percentage) + 0.5, label = "*"),
        color = "black",
        size = 5,
        vjust = 0.8
      ) +
      scale_x_continuous(breaks = c(max_val-2, max_val, max_val+2),
                         expand = c(0,0),
                         labels = c(paste0("||  ",max_val-2), as.character(max_val), as.character(max_val+2)))+
      coord_cartesian(xlim = c(max_val-2, max_val+2)) +
      scale_y_discrete(drop = FALSE)
    
    p <- (p1 | plot_spacer() | p2) +
      plot_layout(widths = c(0.75,0.005,0.25))
  } 
  
  return(p)
})

plot_list[[6]]

saveRDS(plot_list,file.path(figpath,"GO_process_relab_list.rds"))

names(plot_list) <- group_level

combined_plot <- wrap_plots(plot_list[grepl("human", group_level)], ncol = 1) +
  plot_annotation(
    title = "Most Prevalent Biological Process in LGT by Cohort",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  ) 

  
combined_plot

save_plot(
  combined_plot,
  "GO_process_relab_human",
  height = 18,
  width = 15,
  unit = "cm",
  scale = 1,
  dpi=1000
)

combined_plot <- wrap_plots(plot_list[grepl("dog", group_level)], ncol = 1) +
  plot_annotation(
    title = "Most Prevalent Biological Process in LGT by Cohort",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  ) 

combined_plot

save_plot(
  combined_plot,
  "GO_process_relab_dog",
  height = 18,
  width = 15,
  unit = "cm",
  scale = 1,
  dpi=1000
)


#######
# AMR #
#######
# Fix: Add to AMR LIST by function
amr_prefixes <- unique(c(
  annotated_df %>%
    filter(is_amr == TRUE) %>%
    pull(what_amr) %>%
    na.omit() %>%
    tolower() %>%
    unique(), 
  annotated_df %>%
    filter(is_amr == TRUE | grepl("mycin", Consolidated_Function, ignore.case = TRUE)) %>%
    pull(Consolidated_Function) %>%
    na.omit() %>%
    sub("[ _].*$", "", .) %>%
    tolower() %>%
    unique(),
  "streptogramin"
)) 

amr_prefixes <- c(amr_prefixes[!amr_prefixes %in% c("amr", "mfs", "malonyl-coa")], "rifamycin")
amr_pattern  <- paste0("(?i)\\b(", paste(amr_prefixes, collapse = "|"), ")")

annotated_df_fixed <- annotated_df %>%
  mutate(
    .matched = is_amr == TRUE | grepl(amr_pattern, Consolidated_Function, perl = TRUE),
    is_amr = ifelse(.matched, TRUE, is_amr),
    what_amr = ifelse(
      .matched & is.na(what_amr),
      str_extract(Consolidated_Function, regex(amr_pattern, ignore_case = TRUE)),
      what_amr
    ),
    what_amr = tolower(what_amr),
    what_amr = case_when(
      what_amr == "27-o-demethylrifamycin" ~ "rifamycin",
      what_amr == "virginiamycin" ~ "streptogramin",
      TRUE ~ what_amr
    )
  ) %>%
  select(-.matched)

AMR_byType <- annotated_df_fixed %>%
  filter(is_amr == TRUE) %>%
  mutate(B_genes = paste0(B_genes, "\n (", tolower(what_amr),")")) %>%
  group_by(type, B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

AMR_byHost<- annotated_df_fixed %>%
  filter(is_amr == TRUE) %>%
  mutate(B_genes = paste0(B_genes, "\n (", tolower(what_amr),")")) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(host, B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

AMR_All <- annotated_df_fixed %>%
  filter(is_amr == TRUE) %>%
  mutate(B_genes = paste0(B_genes, "\n (", tolower(what_amr),")")) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(relab=100*n_B_genes/n_Sandwich)

summary_AMR <- AMR_All %>% select(B_genes, n_Sandwich) %>% rename(All=n_Sandwich) %>%
  left_join(AMR_byType %>% select(B_genes, type, n_Sandwich) %>%
              pivot_wider(names_from = type, values_from = n_Sandwich, values_fill = 0),
            by="B_genes") %>%
  left_join(AMR_byHost %>% select(B_genes, host, n_Sandwich) %>%
              pivot_wider(names_from = host, values_from = n_Sandwich, values_fill = 0),
            by="B_genes") %>%
  arrange(desc(All)) %>%
  rename(humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF,
         dogOral = DDV,
         dogSkin = DDH) 


select_group <- cols_to_normalize[cols_to_normalize %in% colnames(summary_AMR)]
sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Sandwich, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[select_group]

summary_AMR[, select_group] <- round(sweep(summary_AMR[, select_group], 2, sum_per_cols, `/`) *100,3)


plot_data <- summary_AMR %>%
  rename(feature = B_genes) %>%
  pivot_longer(cols = all_of(select_group), names_to = "Group", values_to = "pct_sandwich") %>%
  mutate(Group = factor(Group, select_group)) %>%
  group_by(Group) %>%
  slice_max(order_by = pct_sandwich, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  filter(Group %in% types) %>%
  filter(pct_sandwich > 0)
    

group_level <- as.character(unique(plot_data$Group))

plot_list <- lapply(group_level, function(x){
  sub_df <- plot_data %>%
    filter(Group == x) %>%
    mutate(
      feature = gsub("UniRef90_", "", feature),
      percentage = round(pct_sandwich, 3)
    ) %>%
    select(Group, feature, percentage) %>%
    arrange(desc(percentage)) %>%
    mutate(feature = factor(feature, levels = rev(unique(feature))))
  
  p <- ggplot(sub_df, aes(x = percentage, y = feature, fill = Group)) +
    geom_col(show.legend = FALSE) +
    theme_minimal() +
    labs(title = x, x="Proportion in Unique LGT (%)") +
    theme(
      legend.position = "inside",
      plot.title = element_text(face = "bold", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 7.5),
      axis.title.y = element_blank()
    ) +
    scale_fill_manual(values = site_col)
  
  return(p)
  })

names(plot_list) <- group_level
target_order <- c("humanFecal", "dogFecal",
                  "humanOral", "dogOral",
                  "humanSkin", "dogSkin")

ordered_plots <- lapply(target_order, function(x) {
  plot_list[[x]] + 
    ylab("UniRef90 (Class)") + 
    theme(axis.title.y = element_text(size = 10, angle = 90))
})
names(ordered_plots) <- target_order

combined_plot <- wrap_plots(ordered_plots) +
  plot_layout(
    ncol = 2,
    heights = c(4, 4, 1.5),
    axis_titles = "collect"
  ) +
  plot_annotation(
    title = "Most Transferred AMR by Cohort",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  )

combined_plot

save_plot(
  combined_plot,
  "AMR_relab_nSandwich",
  height = 8,
  width = 15,
  unit = "cm",
  scale = 1.2,
  dpi=600
)

out2 <- plot_data %>%
  mutate(
    class = sub(".*\n ", "", feature),
    feature = gsub("UniRef90_","",sub("\n.*", "", feature))
  ) %>%
  left_join(annot_byGenes %>% mutate(feature = gsub("UniRef90_","",ANNOTATIONS.UNIREF90)) %>% select(feature, ACCESSION, Consolidated_Function), by="feature") %>%
  mutate(Group = droplevels(factor(Group, levels = types))) %>%
  arrange(Group, desc(pct_sandwich)) %>%
  mutate(across(all_of(c("ACCESSION", "Consolidated_Function")), \(x) ifelse(is.na(x), "--", x))) %>%
  mutate(DB = case_when(
    grepl("K",ACCESSION) ~ "KO",
    grepl("IP",ACCESSION) ~ "IP",
    grepl("PT",ACCESSION) ~ "PT",
    grepl("PF",ACCESSION) ~ "PF",
    grepl("SS",ACCESSION) ~ "SF",
    grepl("G3D",ACCESSION) ~ "GD",
    grepl("Var",ACCESSION) ~ "Various",
    TRUE ~ "--"
  )) %>%
  mutate(feature = as.character(feature),
         feature = gsub("UniRef90_","", feature)) %>%
  select(Group, feature, pct_sandwich, ACCESSION, DB, Consolidated_Function)

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/AMR_genes_by_Sandwich_full.tsv"

writeLines(header_comment, con = file_path)
write.table(out2, file_path,
            quote = F, sep="\t", row.names = F)





#######
# CAZY #
#######
# FiX: Add to CAZY List
cazy_lookup <- c(
  "Six-hairpin glycosidase superfamily" = "GH",
  "maltose-6'-phosphate_glucosidase"   = "GH4"
)

test <- annotated_df_fixed

for (i in names(cazy_lookup)) {
  idx <- which(
    grepl(i, test$Consolidated_Function, ignore.case = TRUE) & 
      (is.na(test$is_cazy) | test$is_cazy == FALSE)
  )
  
  if (length(idx) > 0) {
    test[idx, "what_cazy"] <- cazy_lookup[i]
    test[idx, "is_cazy"]   <- TRUE
  }
}

test[!is.na(test$what_cazy), "is_cazy"] <- TRUE
annotated_df_fixed2 <- test


CAZY_byType <- annotated_df_fixed2 %>%
  filter(is_cazy == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", what_cazy, ")")) %>%
  group_by(type, B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )

CAZY_byHost<- annotated_df_fixed2 %>%
  filter(is_cazy == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", what_cazy, ")")) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(host, B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )

CAZY_All <- annotated_df_fixed2 %>%
  filter(is_cazy == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", what_cazy, ")")) %>%
  mutate(B_genes = ifelse(is.na(B_genes), "No GO Entry", B_genes)) %>%
  mutate(host=case_when(
    type %in% c("DHF","DHV","DHP") ~ "human",
    TRUE ~ "dog"
  )) %>%
  group_by(B_genes) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )

summary_CAZY <- CAZY_All %>% select(B_genes, n_Sandwich) %>% rename(All=n_Sandwich) %>%
  left_join(CAZY_byType %>% select(B_genes, type, n_Sandwich) %>%
              pivot_wider(names_from = type, values_from = n_Sandwich, values_fill = 0),
            by="B_genes") %>%
  left_join(CAZY_byHost %>% select(B_genes, host, n_Sandwich) %>%
              pivot_wider(names_from = host, values_from = n_Sandwich, values_fill = 0),
            by="B_genes") %>%
  rename(humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) %>%
  arrange(desc(All))

select_group <- cols_to_normalize[cols_to_normalize %in% colnames(summary_CAZY)]
sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Sandwich, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[select_group]

summary_CAZY[, select_group] <- round(sweep(summary_CAZY[, select_group], 2, sum_per_cols, `/`) *100,3)


plot_data <- summary_CAZY %>%
  rename(feature = B_genes) %>%
  pivot_longer(cols = all_of(select_group), names_to = "Group", values_to = "pct_sandwich") %>%
  mutate(Group = factor(Group, select_group)) %>%
  group_by(Group) %>%
  slice_max(order_by = pct_sandwich, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  filter(Group %in% types) %>%
  filter(pct_sandwich > 0)


group_level <- as.character(unique(plot_data$Group))

plot_list <- lapply(group_level, function(x){
  sub_df <- plot_data %>%
    filter(Group == x) %>%
    mutate(
      feature = gsub("UniRef90_", "", feature),
      percentage = round(pct_sandwich, 3)
    ) %>%
    select(Group, feature, percentage) %>%
    arrange(desc(percentage)) %>%
    mutate(feature = factor(feature, levels = rev(unique(feature))))
  
  p <- ggplot(sub_df, aes(x = percentage, y = feature, fill = Group)) +
    geom_col(show.legend = FALSE) +
    theme_minimal() +
    labs(title = x, x="Proportion in Unique LGT (%)") +
    theme(
      legend.position = "inside",
      plot.title = element_text(face = "bold", size = 10),
      strip.text = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 7.5),
      axis.title.y = element_blank()
    ) +
    scale_fill_manual(values = site_col)
  
  return(p)
})

names(plot_list) <- group_level
target_order <- c("humanFecal", "dogFecal",
                  "humanOral", "dogOral", 
                  "humanSkin", "dogSkin")

ordered_plots <- lapply(target_order, function(x) {
  p <- plot_list[[x]] + ylab("UniRef90 (Family)") +  theme(axis.title.y = element_text(size=10, angle=90))
  
  return(p)
})
names(ordered_plots) <- target_order


combined_plot <- wrap_plots(ordered_plots, ncol = 2, heights = c(1,1,0.5)) +
  plot_annotation(
    title = "Most Transferred CAZymes by Cohort",
    caption = "GT = glycosyltransferases; GH = glycoside hydrolases;\nPL = polysaccharide lyases; CE = carbohydrate esterases;\nAA = auxiliary activities; CBM = carbohydrate-binding modules",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0)
    )
  ) +
  plot_layout(axis_titles = "collect")

combined_plot

save_plot(
  combined_plot,
  "CAZY_relab_nSandwich",
  height = 18,
  width = 15,
  unit = "cm",
  scale = 1.4,
  dpi=600
)

out3 <- plot_data %>%
  mutate(
    class = sub(".*\n ", "", feature),
    feature = sub("\n.*", "", feature)
  ) %>%
  left_join(annot_byGenes %>% mutate(feature = gsub("UniRef90_","",ANNOTATIONS.UNIREF90)) %>% select(feature, ACCESSION, Consolidated_Function), by="feature") %>%
  mutate(Group = droplevels(factor(Group, levels = types))) %>%
  arrange(Group, desc(pct_sandwich)) %>%
  mutate(across(all_of(c("ACCESSION", "Consolidated_Function")), \(x) ifelse(is.na(x), "--", x))) %>%
  mutate(DB = case_when(
    grepl("K",ACCESSION) ~ "KO",
    grepl("IP",ACCESSION) ~ "IP",
    grepl("PT",ACCESSION) ~ "PT",
    grepl("PF",ACCESSION) ~ "PF",
    grepl("SS",ACCESSION) ~ "SF",
    grepl("G3D",ACCESSION) ~ "GD",
    grepl("Var",ACCESSION) ~ "Various",
    TRUE ~ "--"
  )) %>%
  mutate(feature = as.character(feature),
         feature = gsub("UniRef90_","", feature)) %>%
  select(Group, feature, pct_sandwich, ACCESSION, DB, Consolidated_Function)

header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/CAZY_genes_by_Sandwich_full.tsv"

writeLines(header_comment, con = file_path)
write.table(out3, file_path,
            quote = F, sep="\t", row.names = F)



###################
# Function Survey #
###################
annotated_df_fixed2 <- annotated_df_fixed2 %>%
  mutate(
    Consolidated_Function = ifelse(Consolidated_Function=="--" | Consolidated_Function == "", NA_character_, Consolidated_Function)
  )

AMR_byType <- annotated_df_fixed2 %>%
  filter(is_amr == TRUE)%>%
  mutate(B_genes = paste0(gsub("UniRef90_","",B_genes), "\n (", str_to_sentence(what_amr),")")) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(category="Antimicrobial")

CAZY_byType <- annotated_df_fixed2 %>%
  filter(is_cazy == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", what_cazy, ")")) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="CAZymes")

trp_byType <- annotated_df_fixed2 %>%
  filter(grepl("transposase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Transposase")

trf_byType <- annotated_df_fixed2 %>%
  filter(grepl("transferase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Transferase")

pol_byType <- annotated_df_fixed2 %>%
  filter(grepl("polymerase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA/RNA Polymerase")

rec_byType <- annotated_df_fixed2 %>%
  filter(grepl("recomb|integrase|invertase", tolower(Consolidated_Function)) | grepl("dna binding", tolower(GO_MOLECULAR_FUNCTION))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA Recombinase")

tox_byType <- annotated_df_fixed2 %>%
  filter(grepl("toxin", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Toxin/Antitoxin")

prm_byType <- annotated_df_fixed2 %>%
  filter(grepl("permease", tolower(Consolidated_Function)) | grepl("symporter", tolower(GO_MOLECULAR_FUNCTION))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Permease")

dmg_byType <- annotated_df_fixed2 %>%
  filter((grepl("damage|repair", tolower(Consolidated_Function)) & grepl("dna", tolower(Consolidated_Function)))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA Repair")

vir_byType <- annotated_df_fixed2 %>%
  filter(grepl("endonuclease|crispr|anti-phage|antiviral|abortive infection|prophage|capsid", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Phage|Virus-related")

mov_byType <- annotated_df_fixed2 %>%
  filter(grepl("chemotaxis|flagella|cilia", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Cell Motility")

unk_byType <- annotated_df_fixed2 %>%
  filter(grepl("uncharacterized_protein|unknown", tolower(Consolidated_Function)) | is.na(Consolidated_Function)) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Unknown Function")

functions_df <- bind_rows(AMR_byType, 
                         CAZY_byType,
                         trp_byType,
                         trf_byType,
                         prm_byType,
                         pol_byType,
                         rec_byType,
                         dmg_byType,
                         tox_byType,
                         vir_byType,
                         mov_byType,
                         unk_byType) %>%
  group_by(type, category) %>%
  summarise(n_Sandwich=sum(n_Sandwich), .groups = "drop") %>%
  pivot_wider(
    id_cols = c("category"),
    names_from = type,
    values_from = n_Sandwich,
    values_fill = 0
  ) %>%
  arrange(category) %>%
  select(category, DHF, DHV, DHP, DDF, DDV, DDH) %>%
  mutate(human = DHF + DHV + DHP,
         dog = DDF + DDV + DDH,
         All = human + dog) %>%
  rename(humanFecal = DHF,
         humanOral = DHV,
         humanSkin = DHP,
         dogFecal = DDF, 
         dogOral = DDV,
         dogSkin = DDH) %>%
  arrange(desc(All)) %>%
  mutate(category=factor(category, unique(category)))

sum_per_cols <- bind_rows(total_genes_byType, total_genes_byHost, total_genes_all)
sum_per_cols <- setNames(sum_per_cols$n_Sandwich, sum_per_cols$Group)
sum_per_cols <- sum_per_cols[cols_to_normalize]

functions_df[, cols_to_normalize] <- sweep(functions_df[, cols_to_normalize], 2, sum_per_cols, `/`) *100

write.table(functions_df, "summary_tables/functions_byType.tsv",
            quote = F, sep="\t", row.names = F)


plot_data <- functions_df %>%
  pivot_longer(cols=all_of(cols_to_normalize), names_to = "Group", values_to = "pct_sandwich")

p <- plot_data %>%
  filter(Group %in% types) %>%
  rename(type = Group, percentage = pct_sandwich) %>%
  mutate(
    host = ifelse(grepl("human", type, ignore.case = TRUE), "human", "dog"),
    body_site = gsub("human|dog", "", type, ignore.case = TRUE),
    type = droplevels(factor(type, types)),
    host = factor(host, c("human", "dog")),
    body_site = factor(body_site, c("Fecal","Oral", "Skin"))
  ) %>%
  ggplot(aes(y = category, x = percentage, fill = type)) +
  facet_grid(body_site~host, scale="free_x") +
  geom_col(position = "dodge") +
  geom_text(
    aes(x= percentage +1, label = {
      lab=sprintf("%.2f%%", percentage);
      ifelse(lab=="0.00%","NA",lab)}),
    size = 2.5, 
    angle=0,
    vjust = 0.4,
    hjust=0, 
  ) +
  labs(
    title = "Proportion of LGT with Associated Gene",
    y = "Proportion of Unique LGT (%)",
    fill = "Category",
    caption = "Only LGT with direction is considered; an LGT can have both ARG and CAZymes"
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.3))
  ) +
  scale_y_discrete(expand = c(0.1,0.1), limits=rev)+
  scale_fill_manual(values=site_col) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=0, hjust = 1, size=8),
        panel.border = element_rect(linewidth = 0.5, colour="black"),
        strip.text = element_text(face="bold", size=10),
        panel.spacing.y = unit(0.5, "cm"),
        plot.margin = margin(0.5,0.5,0.5,1, "cm")) +
  guides(fill = guide_legend(nrow = 1))

p


save_plot(p, "LGT_function_prop_byType", width = 15, height = 17, unit = "cm",dpi=600, scale=1.1)


# TABLE OUT
library(dplyr)

# 1. Total LGT events per type
total_lgt_per_type <- annotated_df_fixed2 %>%
  group_by(type) %>%
  summarise(total_sandwiches = n_distinct(sandwichDist), .groups = "drop")

# 2. General helper function using base R eval
get_top10_genes <- function(df, filter_expr, category_name) {
  df %>%
    filter(eval(filter_expr, envir = .)) %>%
    mutate(
      B_genes = gsub("UniRef90_", "", B_genes),
      host = ifelse(grepl("human", type, ignore.case = TRUE), "human", "dog"),
      body_site = gsub("^(human|dog)_?", "", type, ignore.case = TRUE)
    ) %>%
    group_by(body_site, type, B_genes, Consolidated_Function, ACCESSION) %>%
    summarise(n_Sandwich = n_distinct(sandwichDist), .groups = "drop") %>%
    left_join(total_lgt_per_type, by = "type") %>%
    mutate(
      percent_LGT = (n_Sandwich / total_sandwiches) * 100,
      category = category_name
    ) %>%
    filter(percent_LGT > 0) %>% 
    group_by(body_site) %>%
    slice_max(order_by = n_Sandwich, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(DB = case_when(
      grepl("K", ACCESSION)   ~ "KO",
      grepl("IP", ACCESSION)  ~ "IP",
      grepl("PT", ACCESSION)  ~ "PT",
      grepl("PF", ACCESSION)  ~ "PF",
      grepl("SS", ACCESSION)  ~ "SF",
      grepl("G3D", ACCESSION) ~ "GD",
      grepl("Var", ACCESSION) ~ "Various",
      TRUE                    ~ "--"
    )) %>%
    mutate(body_site=substr(body_site, 2,3),
           body_site=gsub("DH","DS",body_site),
           body_site=gsub("HP","HS",body_site),
           body_site=droplevels(factor(body_site, c("HF","HV","HS","DF","DV","DS")))) %>%
    mutate(
      Consolidated_Function = ifelse(
        str_count(ACCESSION, ";") > 1,
        "Various",
        Consolidated_Function
      )
    ) %>%
    arrange(body_site, desc(percent_LGT)) %>%
    select(
      site = body_site,
      UniRef90 = B_genes,
      percent_LGT,
      ACC_No = ACCESSION,
      DB,
      Consolidated_Function,
      category
    )
}

# 3. Define expression filters with base quote
filters <- list(
  Transposase     = quote(grepl("transposase", tolower(Consolidated_Function))),
  Transferase     = quote(grepl("transferase", tolower(Consolidated_Function))),
  Polymerase      = quote(grepl("polymerase", tolower(Consolidated_Function))),
  DNA_Recombinase = quote(grepl("recomb|integrase|invertase", tolower(Consolidated_Function)) | grepl("dna binding", tolower(GO_MOLECULAR_FUNCTION))),
  Viral_related   = quote(grepl("endonuclease|crispr|anti-phage|antiviral|abortive infection|prophage|capsid", tolower(Consolidated_Function))),
  Toxin_Antitoxin = quote(grepl("toxin", tolower(Consolidated_Function))),
  Permease        = quote(grepl("permease", tolower(Consolidated_Function)) | grepl("symporter", tolower(GO_MOLECULAR_FUNCTION))),
  DNA_Repair      = quote(grepl("damage|repair", tolower(Consolidated_Function)) & grepl("dna", tolower(Consolidated_Function))),
  Cell_Movement   = quote(grepl("chemotaxis|flagella|cilia", tolower(Consolidated_Function))),
  CAZymes         = quote(is_cazy == TRUE),
  Antimicrobial   = quote(is_amr == TRUE),
  Unknown         = quote(grepl("uncharacterized_protein|unknown", tolower(Consolidated_Function)) | is.na(Consolidated_Function))
)

# 4. Combine results across all categories
top10_all_categories <- bind_rows(
  lapply(names(filters), function(cat) {
    get_top10_genes(annotated_df_fixed2, filters[[cat]], cat)
  })
) %>%
  mutate(category=factor(category, names(filters))) %>%
  arrange(category, site, desc(percent_LGT))

View(top10_all_categories)
table(top10_all_categories$category)


header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- "summary_tables/All_functions_genes_by_Sandwich_full.tsv"

writeLines(header_comment, con = file_path)
write.table(top10_all_categories, file_path,
            quote = F, sep="\t", row.names = F)


write.table(annotated_df_fixed2, "input_file/annotated_cargo_df_fixed.tsv",
            quote = F, sep="\t", row.names = F)


# by contig
AMR_byType <- annotated_df_fixed2 %>%
  filter(is_amr == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", str_to_sentence(what_amr), ")")) %>%
  group_by(type) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(paste(CONTIG_NAME, SAMPLE, sep = "+")),
    .groups = "drop"
  ) %>%
  mutate(category = "Antimicrobial") %>%
  pivot_wider(
    id_cols = category,
    names_from = type,
    values_from = n_Contigs,
    values_fill = 0
  )

expected_cols <- c("DHF", "DHV", "DHP", "DDF", "DDV", "DDH")
missing_cols <- setdiff(expected_cols, colnames(AMR_byType))

for (col in missing_cols) {
  AMR_byType[[col]] <- 0
}

AMR_byType <- AMR_byType %>%
  select(category, all_of(expected_cols)) %>%
  mutate(
    human = DHF + DHV + DHP,
    dog = DDF + DDV + DDH,
    All = human + dog
  ) %>%
  rename(
    humanFecal = DHF,
    humanOral = DHV,
    humanSkin = DHP,
    dogFecal = DDF,
    dogOral = DDV,
    dogSkin = DDH
  ) %>%
  arrange(desc(All)) %>%
  mutate(category = factor(category, unique(category)))

nContigs_byType <- nContigs %>%
  mutate(type = substr(SAMPLE, 1, 3)) %>%
  mutate(type = case_when(
    type == "DHF" ~ "humanFecal",
    type == "DHV" ~ "humanOral",
    type == "DHP" ~ "humanSkin",
    type == "DDF" ~ "dogFecal",
    type == "DDV" ~ "dogOral",
    type == "DDH" ~ "dogSkin"
  )) %>%
  group_by(type) %>%
  summarise(total_contig_mio = sum(n_contig) / 1e6, .groups = "drop")

sum_per_cols <- setNames(nContigs_byType$total_contig_mio, nContigs_byType$type)

sum_per_cols <- c(
  sum_per_cols,
  All = sum(sum_per_cols, na.rm = TRUE),
  human = sum(sum_per_cols[grepl("^human", names(sum_per_cols))], na.rm = TRUE),
  dog = sum(sum_per_cols[grepl("^dog", names(sum_per_cols))], na.rm = TRUE)
)

cols_to_normalize <- intersect(names(sum_per_cols), colnames(AMR_byType))

sum_per_cols <- sum_per_cols[cols_to_normalize]

AMR_byType[, cols_to_normalize] <- sweep(AMR_byType[, cols_to_normalize], 2, sum_per_cols, `/`)
round(AMR_byType[,-1],2) 

################################
## Functional cat by samples ##
###############################
AMR_bySample <- annotated_df_fixed2 %>%
  filter(is_amr == TRUE)%>%
  mutate(B_genes = paste0(gsub("UniRef90_","",B_genes), "\n (", str_to_sentence(what_amr),")")) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )  %>%
  mutate(category="Antimicrobial")

CAZY_bySample <- annotated_df_fixed2 %>%
  filter(is_cazy == TRUE) %>%
  mutate(B_genes = paste0(gsub("UniRef90_", "", B_genes), "\n (", what_cazy, ")")) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="CAZymes")

trp_bySample <- annotated_df_fixed2 %>%
  filter(grepl("transposase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Transposase")

trf_bySample <- annotated_df_fixed2 %>%
  filter(grepl("transferase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Transferase")

pol_bySample <- annotated_df_fixed2 %>%
  filter(grepl("polymerase", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA/RNA Polymerase")

rec_bySample <- annotated_df_fixed2 %>%
  filter(grepl("recomb|integrase|invertase", tolower(Consolidated_Function)) | grepl("dna binding", tolower(GO_MOLECULAR_FUNCTION))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA Recombinase")

tox_bySample <- annotated_df_fixed2 %>%
  filter(grepl("toxin", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Toxin/Antitoxin")

prm_bySample <- annotated_df_fixed2 %>%
  filter(grepl("permease", tolower(Consolidated_Function)) | grepl("symporter", tolower(GO_MOLECULAR_FUNCTION))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Permease")

dmg_bySample <- annotated_df_fixed2 %>%
  filter((grepl("damage|repair", tolower(Consolidated_Function)) & grepl("dna", tolower(Consolidated_Function)))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="DNA Repair")

mov_bySample <- annotated_df_fixed2 %>%
  filter(grepl("chemotaxis|flagella|cilia", tolower(Consolidated_Function))) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Cell Motility")


unk_bySample <- annotated_df_fixed2 %>%
  filter(grepl("uncharacterized_protein|unknown", tolower(Consolidated_Function)) | is.na(Consolidated_Function)) %>%
  mutate(B_genes = gsub("UniRef90_","",B_genes)) %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  ) %>%
  mutate(category="Unknown Function")

n_distinct_sample <- annotated_df_fixed2 %>%
  group_by(SAMPLE) %>%
  summarise(
    n_Sample = n_distinct(SAMPLE),
    n_B_genes = n_distinct(B_genes),
    n_Sandwich = n_distinct(sandwichDist),
    n_Contigs = n_distinct(index),
    .groups = "drop"
  )
  
functions_df <- bind_rows(AMR_bySample, 
                          CAZY_bySample,
                          trp_bySample,
                          trf_bySample,
                          prm_bySample,
                          pol_bySample,
                          rec_bySample,
                          dmg_bySample,
                          tox_bySample, 
                          #mov_bySample,
                          unk_bySample) %>%
  group_by(SAMPLE, category) %>%
  summarise(n_Sandwich=sum(n_Sandwich), .groups = "drop") %>%
  left_join(n_distinct_sample %>%
              select(SAMPLE, sample_total=n_Sandwich), by="SAMPLE") %>%
  mutate(n_Sandwich = round(100*n_Sandwich/sample_total,1)) %>%
  pivot_wider(
    id_cols = c("category"),
    names_from = SAMPLE,
    values_from = n_Sandwich,
    values_fill = 0
  ) %>%
  mutate(All = rowSums(across(-category))) %>%
  arrange(desc(All)) %>%
  select(-All) %>%
  #bind_rows(
  #summarise(., across(-category, ~ sum(.x) / 100)) %>%
  #    mutate(category = "Others")
  #) %>%
  mutate(category=factor(category, unique(category)))

write.table(functions_df, "summary_tables/functions_bySample.tsv",
            quote = F, sep="\t", row.names = F)

bar_data <- functions_df %>%
  pivot_longer(
    cols = -category,
    names_to = "SAMPLE",
    values_to = "pct_sandwich"
  ) %>%
  mutate(short_sample = substr(SAMPLE, 1,6)) %>%
  left_join(pop.df %>% select(-SAMPLE) %>% distinct(), by="short_sample") %>%
  mutate(type = substr(SAMPLE, 1, 3)) %>%
  mutate(type = case_when(
    type == "DHF" ~ "humanFecal",
    type == "DHV" ~ "humanOral",
    type == "DHP" ~ "humanSkin",
    type == "DDF" ~ "dogFecal",
    type == "DDV" ~ "dogOral",
    type == "DDH" ~ "dogSkin"
  )) %>%
  mutate(type = droplevels(factor(type, types)))

sample_ord <- bar_data %>%
  group_by(type, SAMPLE) %>%
  summarise(
    sum_all = sum(pct_sandwich),
    sum_transp = sum(pct_sandwich[category == "Transposase"]),
    .groups = "drop"
  ) %>%
  arrange(type, desc(sum_all), desc(sum_transp)) %>%
  pull(SAMPLE)


library(RColorBrewer)

bar_p <- bar_data %>%
  mutate(SAMPLE = droplevels(factor(SAMPLE, levels = sample_ord))) %>%
  ggplot(aes(x = SAMPLE, y = pct_sandwich, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  facet_wrap(~type, scales = "free_x") +
  scale_y_continuous(limits = c(0, 180), expand = c(0,0)) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 0.8),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  )

bar_p

save_plot(bar_p, "functions_cat_bySample_barplot", width = 15, height = 12, unit = "cm", dpi=1200, scale=1.2)
