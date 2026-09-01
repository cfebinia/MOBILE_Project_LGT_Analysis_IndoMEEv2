rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)

# Paths
# my_lib <- "/home/caf77/R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))

basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
# setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir <- file.path(basedir, "index")

out_dir <- file.path(basedir, "summary_tables")
dir.create(out_dir,showWarnings = F,recursive = T)

figpath <- file.path(basedir, "figures")
dir.create(figpath,showWarnings = F,recursive = T)

waffle_file="raw_all_internal_merged_filtered.lgt.tsv"

################
## Load Data  ##
################
type <- c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")
df <- do.call(rbind, lapply(type, function(x) read.delim(file.path(input_dir,x,waffle_file))))

df$original_sid <- df$SAMPLE
df$SAMPLE <- substr(df$original_sid, 1, 6)
df$type <- substr(df$SAMPLE, 1, 3)

synteny_index_table <- read.delim(file.path(index_dir, "synteny_index_table_unAnnot.tsv"))
gene_index_table <- read.delim(file.path(index_dir, "gene_index_table_unAnnot.tsv"))
synteny_position_table <- read.delim(file.path(index_dir, "synteny_position_table_unAnnot.tsv"))
gene_position_table <- read.delim(file.path(index_dir, "gene_position_table_unAnnot.tsv"))

####################
## Synteny Stats  ##
####################
count_synteny <- df %>%
  group_by(SYNTENY) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Proportion = round(Freq / nrow(df), 5) * 100) %>%
  inner_join(synteny_index_table, by = "SYNTENY")

s_list <- strsplit(as.character(count_synteny$SYNTENY), "")
count_synteny$SYNTENY_LEN <- nchar(count_synteny$SYNTENY)
count_synteny$B_len <- sapply(s_list, function(x) sum(x == "B"))
count_synteny$star_len <- sapply(s_list, function(x) sum(x == "*"))
count_synteny$tilde_len <- sapply(s_list, function(x) sum(x == "~"))

calc_max_block <- function(s, char) {
  blocks <- unlist(strsplit(s, paste0("[^", char, "]")))
  if(length(blocks) == 0) return(0)
  max(nchar(blocks))
}
count_synteny$Max_A_len <- sapply(count_synteny$SYNTENY, calc_max_block, char = "A")
count_synteny$Max_B_len <- sapply(count_synteny$SYNTENY, calc_max_block, char = "B")

# Summarise Synteny 
# Calculate nestedA
count_synteny$nestedA <- sapply(count_synteny$SYNTENY, function(s) {
  b_pos <- which(strsplit(s, "")[[1]] == "B")
  if(length(b_pos) < 2) return(FALSE) 
  inner_region <- substr(s, min(b_pos), max(b_pos))
  grepl("A", inner_region)
})

# Export Synteny Summary
synteny_summary <- count_synteny %>%
  select(SYNTENY, Synteny_Index, DIRECTION, Freq, Proportion, 
         SYNTENY_LEN, B_len, star_len, tilde_len, 
         Max_A_len, Max_B_len, nestedA)

write.table(synteny_summary, 
            file = file.path(out_dir, "synteny_summary_unAnnot.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Threshold Visualization ---
b_dist <- count_synteny %>% group_by(Max_B_len) %>% summarise(Total_Events = sum(Freq))
pdf(file.path(figpath, "Max_B_len_threshold_check_unAnnot.pdf"), width = 7, height = 5)
barplot(b_dist$Total_Events, names.arg = b_dist$Max_B_len, col = "skyblue",
        main = "Distribution of LGT Block Lengths", xlab = "Max_B_len", ylab = "Frequency")
abline(v = 10, col = "red", lty = 2)
dev.off()

#####################
## Clade Counting  ##
#####################
count_clade <- df %>% 
  rename(Gene_SYNTENY = ANNOTATIONS.UNIREF90) %>%
  select(DIRECTION, SYNTENY, LOCI)

split_synteny <- strsplit(count_clade$SYNTENY, split = "")

count_clade$SYNTENY_LEN <- lengths(split_synteny)
count_clade$B_LEN       <- sapply(split_synteny, function(x) sum(x == "B"))
count_clade$AMB_LEN     <- sapply(split_synteny, function(x) sum(x %in% c("*", "~")))
count_clade$AMB_FRAC    <- round(count_clade$AMB_LEN / count_clade$SYNTENY_LEN,2)

Blen_tresh <- 10
Slen_tresh <- 35 
Ambg_tresh <- 0.1

table(count_clade$SYNTENY_LEN <= Slen_tresh)
table(count_clade$B_LEN <= Blen_tresh)
table(round(count_clade$AMB_FRAC,1))
table(count_clade$AMB_FRAC <= Ambg_tresh)

table(count_clade$SYNTENY_LEN <= Slen_tresh & count_clade$B_LEN <= Blen_tresh & count_clade$AMB_FRAC <= Ambg_tresh)

df <- df %>%
  mutate(to_exclude=ifelse(count_clade$SYNTENY_LEN <= Slen_tresh & count_clade$B_LEN <= Blen_tresh & count_clade$AMB_FRAC <= Ambg_tresh,
                           "NA","exclude"))
  
table(df$to_exclude)

df <- df %>% filter(to_exclude=="NA") %>% select(-to_exclude)

############################
## Sandwich Construction  ##
############################
unique_synteny <- unique(df$SYNTENY)

table(unique_synteny %in% synteny_index_table$SYNTENY)
table(unique_synteny %in% df$SYNTENY)
table(unique_synteny %in% gene_index_table$SYNTENY)

limit_Blen <- synteny_index_table$Synteny_Index[synteny_index_table$SYNTENY %in% unique_synteny]

# downsample to debug
#set.seed(7) #debug
#limit_Blen <- sample(unique(limit_Blen), 500) #debug

pos_sub <- synteny_position_table %>% filter(Synteny_Index %in% limit_Blen)
syn_sub <- synteny_index_table %>% filter(Synteny_Index %in% limit_Blen)
gene_sub <- gene_index_table %>% filter(Synteny_Index  %in% limit_Blen)
genepos_sub <- gene_position_table %>% filter(Gene_Synteny_Index %in% gene_sub$Gene_Synteny_Index)

# index the df
df_indexed <- df %>%
  mutate(no = row_number()) %>%
  rename(Gene_SYNTENY = ANNOTATIONS.UNIREF90)

# Join the mapping to the indexed dataframe
df_indexed <- df_indexed %>%
  left_join(syn_sub %>% select(-DIRECTION), by = "SYNTENY", keep = FALSE)

# helper function
find_flanker <- function(target_pos, direction, pos_data) {
  a_positions <- pos_data$Position[pos_data$Character == "A"]
  
  if (direction == -1) {
    valid_candidates <- a_positions[a_positions < target_pos]
    if (length(valid_candidates) > 0) {
      return(max(valid_candidates))
    }
  } else if (direction == 1) {
    valid_candidates <- a_positions[a_positions > target_pos]
    if (length(valid_candidates) > 0) {
      return(min(valid_candidates))
    }
  }
  
  return(NA_real_)
}

extract_gene <- function(idx, g_list) {
  if (is.na(idx)) return(rep(NA_character_, length(g_list)))
  vapply(g_list, function(x) if (length(x) >= idx) x[idx] else NA_character_, character(1))
}

extract_locus_coord <- function(idx, loc_list, coord_pos) {
  if (is.na(idx)) return(rep(NA_real_, length(loc_list)))
  vapply(loc_list, function(x) {
    if (length(x) >= idx && !is.na(x[idx])) {
      as.numeric(strsplit(x[idx], ":")[[1]][coord_pos])
    } else {
      NA_real_
    }
  }, numeric(1))
}

modify_unknown_genes <- function(data) {
  data$Gene_SYNTENY <- as.character(data$Gene_SYNTENY)
  data$LOCI <- as.character(data$LOCI)
  data$SYNTENY <- as.character(data$SYNTENY)
  data$CLADE_A <- as.character(data$CLADE_A)
  data$CLADE_B <- as.character(data$CLADE_B)
  
  for (i in 1:nrow(data)) {
    genes <- unlist(strsplit(data$Gene_SYNTENY[i], "\\|"))
    loci <- unlist(strsplit(data$LOCI[i], "\\|"))
    
    if (length(genes) == 0) {
      next
    }
    
    target_indices <- which(genes == "UniRef90_unknown")
    
    if (length(target_indices) > 0) {
      loci_split <- strsplit(loci, ":")
      starts <- as.numeric(sapply(loci_split, "[", 1))
      ends <- as.numeric(sapply(loci_split, "[", 2))
      
      clade <- paste0("CLADE_", unlist(strsplit(data$SYNTENY[i], "")))
      
      spA <- gsub("__", "_", data$CLADE_A[i]) %>% strsplit(split = "_") %>% unlist()
      spB <- gsub("__", "_", data$CLADE_B[i]) %>% strsplit(split = "_") %>% unlist()
      
      if(length(spA) < 3){spA[3] <- "sp"}
      if(length(spB) < 3){spB[3] <- "sp"}
      
      spA_clean <- paste0(str_to_title(substr(spA[-1], 1, 3)), collapse = "")
      spB_clean <- paste0(str_to_title(substr(spB[-1], 1, 3)), collapse = "")
      
      spec_map <- c("CLADE_A" = spA_clean, "CLADE_B" = spB_clean)
      spec <- spec_map[clade]
      
      if (length(spec) != length(genes)) {
        spec <- rep_len(spec, length(genes))
      }
      
      for (idx in target_indices) {
        length_suffix <- ends[idx] - starts[idx]
        genes[idx] <- paste0(genes[idx], "_", length_suffix, "_", spec[idx])
      }
    }
    
    data$Gene_SYNTENY[i] <- paste(genes, collapse = "|")
  }
  
  return(data)
}

modify_all_genes <- function(data) {
  data$Gene_SYNTENY <- as.character(data$Gene_SYNTENY)
  data$LOCI <- as.character(data$LOCI)
  data$SYNTENY <- as.character(data$SYNTENY)
  data$CLADE_A <- as.character(data$CLADE_A)
  data$CLADE_B <- as.character(data$CLADE_B)
  
  for (i in 1:nrow(data)) {
    genes <- unlist(strsplit(data$Gene_SYNTENY[i], "\\|"))
    loci <- unlist(strsplit(data$LOCI[i], "\\|"))
    
    if (length(genes) == 0) {
      next
    }
    
    loci_split <- strsplit(loci, ":")
    starts <- as.numeric(sapply(loci_split, "[", 1))
    ends <- as.numeric(sapply(loci_split, "[", 2))
    
    for (idx in 1:length(genes)) {
      length_suffix <- ends[idx] - starts[idx]
      genes[idx] <- paste0(genes[idx], "_", length_suffix)
    }
    
    clade <- paste0("CLADE_", unlist(strsplit(data$SYNTENY[i], "")))
    
    spA <- gsub("__", "_", data$CLADE_A[i]) %>% strsplit(split = "_") %>% unlist()
    spB <- gsub("__", "_", data$CLADE_B[i]) %>% strsplit(split = "_") %>% unlist()
    
    if(length(spA) < 3){spA[3] <- "sp"}
    if(length(spB) < 3){spB[3] <- "sp"}
    
    spA_clean <- paste0(str_to_title(substr(spA, 1, 3)), collapse = "")
    spB_clean <- paste0(str_to_title(substr(spB, 1, 3)), collapse = "")
    
    spec_map <- c("CLADE_A" = spA_clean, "CLADE_B" = spB_clean)
    spec <- spec_map[clade]
    
    if (length(spec) != length(genes)) {
      spec <- rep_len(spec, length(genes))
    }
    
    genes <- paste(genes, spec, sep = "_")
    
    data$Gene_SYNTENY[i] <- paste(genes, collapse = "|")
  }
  
  return(data)
}
# si=567

sandwich_df <- lapply(unique(df_indexed$Synteny_Index), function(si) {
  this_pos <- pos_sub[pos_sub$Synteny_Index == si, ]
  b_indices <- this_pos$Position[this_pos$Character == "B"]
  
  direction_condition <- syn_sub$DIRECTION[syn_sub$Synteny_Index == si]
  direction_BtoA <- length(direction_condition) > 0 && direction_condition[1] == "B>A"
  
  a1_idx <- find_flanker(min(b_indices), -1, this_pos)
  a2_idx <- if (direction_BtoA) find_flanker(max(b_indices), 1, this_pos) else NA
  
  if (direction_BtoA) {
    sandwich_synteny <- paste0(c("A", rep("B", length(b_indices)), "A"), collapse = "")
  } else {
    sandwich_synteny <- paste0(c("A", rep("B", length(b_indices))), collapse = "")
  }
  
  genes_sub <- df_indexed[df_indexed$Synteny_Index == si, c("no", "SYNTENY","CLADE_A","CLADE_B","Gene_SYNTENY", "LOCI")]
  
  genes_sub <- modify_unknown_genes(genes_sub)
  # genes_sub <- modify_all_genes(genes_sub)
  
  g_list <- strsplit(genes_sub$Gene_SYNTENY, "\\|")
  loc_list <- strsplit(genes_sub$LOCI, "\\|")

  
  a1_genes <- extract_gene(a1_idx, g_list)
  a2_genes <- extract_gene(a2_idx, g_list)
  
  b_genes_list <- vapply(g_list, function(g) paste(g[b_indices], collapse = "|"), character(1))
  
  a1_ends <- extract_locus_coord(a1_idx, loc_list, 2)
  target_start_idx <- if (direction_BtoA) a2_idx else a1_idx + 1
  a2_starts <- extract_locus_coord(target_start_idx, loc_list, 1)
  
  b_ends <- extract_locus_coord(max(b_indices), loc_list, 2)
  b_starts <- extract_locus_coord(min(b_indices), loc_list, 1)
  
  A_gapLen <- a2_starts - a1_ends
  B_len <- b_ends - b_starts
  
  data.frame(
    no = genes_sub$no,
    sandwich_synteny = sandwich_synteny,
    B_count = length(b_indices),
    A1_gene = a1_genes,
    A2_gene = a2_genes,
    B_genes = b_genes_list,
    A_gapLen = A_gapLen,
    B_Len = B_len,
    stringsAsFactors = FALSE
  )
  
}) %>% 
  bind_rows() %>% 
  rowwise() %>% 
  mutate(sandwich = paste0(na.omit(c(A1_gene, B_genes, A2_gene)), collapse = "|")) %>%
  ungroup() %>% 
  arrange(no)

df_indexed <- df_indexed %>%
  right_join(sandwich_df, by = "no") %>%
  arrange(Synteny_Index, B_count, sandwich) %>%
  mutate(Sandwich_Index = match(sandwich, unique(sandwich))) %>%
  arrange(no)

# eliminate reversed
df_indexed$no_reversed_sandwich <- sapply(strsplit(df_indexed$sandwich, split = "\\|"), function(x) {
  original <- paste(x, collapse = "_")
  reversed <- paste(rev(x), collapse = "_")
  if (original >= reversed) {
    return(original)
  } else {
    return(reversed)
  }
})

df_indexed$no_reversed_Bgenes <- sapply(strsplit(df_indexed$B_genes, split = "\\|"), function(x) {
  original <- paste(x, collapse = "_")
  reversed <- paste(rev(x), collapse = "_")
  if (original >= reversed) {
    return(original)
  } else {
    return(reversed)
  }
})

table(df_indexed$no_reversed_sandwich == df_indexed$no_reversed_sandwich)

df_indexed <- df_indexed %>%
  distinct(SAMPLE,CONTIG_NAME,no_reversed_sandwich, A_gapLen, .keep_all = T) %>%
  select(-no_reversed_sandwich, -no_reversed_Bgenes)

# check completeness
table(is.na(df_indexed$sandwich))
table(is.na(df_indexed$A1_gene))
table(is.na(df_indexed$B_genes))
table(df_indexed$DIRECTION,is.na(df_indexed$A2_gene))


# write indexed output
write.table(df_indexed, file.path(out_dir, "all_samples_internal_merged_indexed_filtered.tsv"), 
            sep="\t", quote=F, row.names=F)


