gc()
rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggrepel)

# basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

summary_dir=file.path(basedir, "summary_tables")
figpath=file.path(basedir, "figures")
out_dir=file.path(basedir, "summary_tables")

waafle_file <- "summary_tables/all_samples_internal_merged_indexed_filtered.tsv"

nContigs_file <- "input_file/input_contigs_compile.tsv"

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
  "humanFecal"="#0A5FCD", 
  "humanSkin"="#EE2C2C", 
  "humanOral"="#c59f00ff", 
  "dogFecal"="#0b0b8f", 
  "dogSkin"="#780247", 
  "dogOral"="#9A5324",
  "human" = "grey70",
  "dog" = "grey40",
  "All" = "grey10"
)


################
## Load Data  ##
################
cutoff = 0.05

#input_df = "Unique_LGT_glmer_output_test.tsv"
#figout = "lmer_volcano_lifestyle_LGT"
#this_model = "prevalence[0/1] ~ lifestyle_group + log10(total_LGT) + (1| community)"
#my_title = "Association: LGT prevalence vs Lifestyle group"

#input_df = "Unique_LGT_glmer_output_ordinalLS_twoset.tsv"
#figout = "lmer_volcano_lifestyle_LGT_ordinalLS"
#this_model = "prevalence[0/1] ~ lifestyle_cline2 + log10(total_LGT) + (1|community)"
#my_title = "Association: LGT prevalence vs Lifestyle (Ordinal)"

input_df = "Unique_LGT_glmer_output_ordinalLS.tsv"
figout = "lmer_volcano_lifestyle_LGT_ordinalLS"
this_model = "prevalence[0/1] ~ lifestyle_cline2 + log10(total_LGT) + (1|community)"
my_title = "Association: LGT prevalence vs Lifestyle (Ordinal)"

#input_df = "Unique_Cargo_glmer_output_test.tsv"
#figout = "lmer_volcano_lifestyle_Cargo"
#this_model = "prevalence[0/1] ~ lifestyle_group + log10(total_genes) + (1| community)"
#my_title = "Association: Gene prevalence vs Lifestyle"

#input_df = "Unique_Cargo_glmer_output_ordinalLS.tsv"
#figout = "lmer_volcano_lifestyle_Cargo_ordinal"
#this_model = "prevalence[0/1] ~ lifestyle_cline2 + log10(total_genes) + (1| community)"
#my_title = "Association: Gene prevalence vs Lifestyle (Ordinal)"


glmer_out <- read.delim(file.path(summary_dir, input_df)) %>%
  filter(grepl("lifestyle", metadata) & not0 > 5) %>%
  filter(
    p_val <= 0.05 | 
      between(coef, quantile(coef, 0.05, na.rm = TRUE), quantile(coef, 0.95, na.rm = TRUE))
  ) %>%
  mutate(feature = gsub("UniRef90_", "", feature)) %>%
  group_by(type) %>%
  mutate(q_val = p.adjust(p_val, method = "fdr", n = n())) %>%
  ungroup() %>%
  mutate(
    site = gsub("human|dog", "", type),
    host = gsub("Fecal|Oral|Skin", "", type),
    site = factor(site, levels = c("Fecal", "Oral", "Skin")),
    host = factor(host, levels = c("human", "dog")),
    type = droplevels(factor(type, levels = names(site_col))),
    significant = p_val < cutoff,
    colour_points = case_when(
      p_val < cutoff ~ as.character(type),
      TRUE ~ "not_significant"
    )) %>%
  mutate(colour_points=droplevels(factor(colour_points, c(names(site_col),"not_significant")))) %>%
  arrange(type,p_val)

sig_pvals <- glmer_out$p_val[glmer_out$q_val < 0.05 & !is.na(glmer_out$q_val)]
q_thresh_y <- if (length(sig_pvals) > 0) -log10(max(sig_pvals, na.rm = TRUE)) else NA_real_

hlines <- c(-log10(0.05), -log10(0.01), -log10(0.001)) #, q_thresh_y)
hlines <- hlines[!is.na(hlines)]

nfeatures <- length(unique(glmer_out$feature))

v_plot <- glmer_out %>%
  ggplot(aes(x = coef, y = -log10(p_val), colour = colour_points)) +
  geom_hline(yintercept = hlines, colour = "darkgrey", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "darkgrey", linewidth = 0.5) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(
    #data = subset(glmer_out, p_val < cutoff & !is.na(p_val)),
    data = glmer_out %>% 
      group_by(type) %>%
      slice_min(order_by = p_val, n=5) %>%
      mutate(label = paste(ifelse(q_val < 0.05, "*",""), feature, sep=""),
             nudge_x = coef*1.1),
    aes(label = label, colour = type),
    fontface="bold",
    size = 2.5,
    point.padding = 0,
    box.padding = 0,
    min.segment.length = 0,
    nudge_y = 0.15,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  #facet_grid(host ~ site) +
  #scale_colour_manual(values = c("FALSE" = "darkgrey", "TRUE" = "red")) +
  scale_colour_manual(values = c(site_col, not_significant="grey70"), drop = FALSE) +
  labs(
    title = my_title,
    x = "Coefficient",
    y = "-log10(p-value)",
    colour = "Sample_type",
    caption = paste("--- : p < 0.05, p < 0.01, p < 0.001;\n",
    "* = q < 0.05 with FDR, across", nfeatures, "features;\n",
    "Model:", this_model)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) + 
  scale_x_continuous(expand = c(0.3, 0.3))

v_plot

save_plot(v_plot, figout,
          width = 15, height=15, unit="cm", dpi=1000,
          scale=1.2)


################
# Annot Signif #
################


df_out <- glmer_out %>%
  filter(p_val < 0.05) %>%
  group_by(type) %>%
  slice_min(order_by = p_val, n = 10) %>%
  ungroup()


if (grepl("len_",df_out$feature)[1]){
  get_B_for_sandiwch <- read.delim(waafle_file) %>%
    select(A_gapLen, A1_gene, B_genes, A2_gene) %>%
    distinct() %>%
    rowwise() %>%
    mutate(sandwichDist = paste0("len_", A_gapLen, "_", 
                                 paste(na.omit(c(A1_gene, B_genes, A2_gene)), collapse = "|"))) %>%
    ungroup() %>%
    mutate(feature=gsub("UniRef90_","",sandwichDist)) %>%
    select(B_genes,feature) %>%
    distinct()
  
  df_out <- df_out %>%
    mutate(LGT_len = sapply(strsplit(feature, split = "_"), `[[`, 2)) %>%
    left_join(get_B_for_sandiwch) %>% 
    mutate(index = seq_along(feature))
  
  rm(get_B_for_sandiwch)
} else {
  df_out <- df_out %>%
    mutate(B_genes = feature)
}


library(stringr)
annot <- readRDS(master_all_annot) %>%
  rename(B_genes = ANNOTATIONS.UNIREF90_fixed) %>%
  select(B_genes, is_amr, is_cazy, ACCESSION, Consolidated_Function, GO_MOLECULAR_FUNCTION) %>%
  mutate(B_genes=gsub("UniRef90_","",B_genes)) %>%
  distinct() %>%
  group_by(B_genes) %>%
  summarise(
    is_amr = any(isTRUE(is_amr)),
    is_cazy= any(isTRUE(is_cazy)),
    ACCESSION={
      v <- unique(ACCESSION)
      v[v==""] <- NA
      v <- na.omit(v)
      str_flatten(v, collapse = "; ")
    },
    Consolidated_Function={
      v <- unique(Consolidated_Function)
      v[v==""] <- NA
      v <- na.omit(v)
      str_flatten(v, collapse = "; ")
    },
    GO_MOLECULAR_FUNCTION={
      v <- unique(GO_MOLECULAR_FUNCTION)
      v[v==""] <- NA
      v <- na.omit(v)
      str_flatten(v, collapse = "; ")
    },
    .groups = "drop"
  )


df_out <- df_out %>%
  separate_longer_delim(B_genes, delim = "|") %>%
  mutate(B_genes=gsub("UniRef90_","",B_genes)) %>%
  left_join(annot) %>%
  arrange(type, p_val)

amr_prefixes <- unique(c(
  df_out %>%
    filter(is_amr == TRUE | grepl("mycin", Consolidated_Function, ignore.case = TRUE)) %>%
    pull(Consolidated_Function) %>%
    na.omit() %>%
    unique(),
  "streptogramin"
)) 

amr_prefixes <- c(amr_prefixes[!amr_prefixes %in% c("amr", "mfs", "malonyl-coa")], "rifamycin")
amr_pattern  <- paste0("(?i)\\b(", paste(amr_prefixes, collapse = "|"), ")")
df_out$is_amr[grepl("mycin", df_out$Consolidated_Function)] <- "TRUE"

filters <- list(
  Transposase     = quote(grepl("transposase", tolower(Consolidated_Function))),
  Transferase     = quote(grepl("transferase", tolower(Consolidated_Function))),
  Polymerase      = quote(grepl("polymerase", tolower(Consolidated_Function))),
  DNA_Recombinase = quote(grepl("recomb|integrase|invertase", tolower(Consolidated_Function)) | grepl("dna binding", tolower(GO_MOLECULAR_FUNCTION))),
  Toxin_Antitoxin = quote(grepl("toxin", tolower(Consolidated_Function))),
  Permease        = quote(grepl("permease", tolower(Consolidated_Function)) | grepl("symporter", tolower(GO_MOLECULAR_FUNCTION))),
  DNA_Repair      = quote(grepl("damage|repair", tolower(Consolidated_Function)) & grepl("dna", tolower(Consolidated_Function))),
  Cell_Movement   = quote(grepl("chemotaxis|flagella|cilia", tolower(Consolidated_Function))),
  CAZymes         = quote(is_cazy == TRUE),
  Antimicrobial   = quote(is_amr == TRUE),
  Unknown         = quote(grepl("uncharacterized_protein|unknown", tolower(Consolidated_Function)) | is.na(Consolidated_Function))
)

filter_matrix <- sapply(filters, function(f) {
  res <- eval(f, envir = df_out)
  ifelse(is.na(res), FALSE, res)
})

df_out$Category <- apply(filter_matrix, 1, function(row) {
  matches <- names(row)[row]
  if (length(matches) == 0) "Other" else paste(matches, collapse = ";")
})

out <- df_out %>%
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
  select(type, B_genes, not0, coef, p_val, q_val, ACCESSION, DB, Consolidated_Function, Category) %>%
  rowwise()%>%
  mutate(
    ACCESSION=ifelse(str_count(ACCESSION, "; ") > 2, "Various", ACCESSION),
    Consolidated_Function=ifelse(str_count(Consolidated_Function, "; ") > 2, "Various", Consolidated_Function)
  ) %>%
  mutate(
    ACCESSION=ifelse(ACCESSION == "", "--", ACCESSION),
    Consolidated_Function=ifelse(Consolidated_Function == "", "--", Consolidated_Function)
  ) %>%
  ungroup() %>%
  arrange(type, p_val)


header_comment <- "# Note: HF = human faecal, HV = human oral, HS = human skin, DF = dog fecal, DV = dog oral, DS = dog skin, DB=Database, KO=Kegg Orthologs, IP=InterPro, PF=Pfam;,PT=PANTHER, SF=SUPFAM, GD=Gene3D."
file_path <- paste0("summary_tables/", figout, "_gene_table.tsv")

writeLines(header_comment, con = file_path)
write.table(out, file_path,
            quote = F, sep="\t", row.names = F)


