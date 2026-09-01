library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)


basedir <- "Imaging_Lab_PC1/WAAFLE_Extra/"
setwd(basedir)
annotation_file <- "input_file/annotations/master_all_annotations.tsv"
cazy_info <- "input_file/annotations/CAZY_db/fam-substrate-mapping.tsv"
cazy_subfam_info <- "input_file/annotations/CAZY_db/CAZyID_subfam_mapping.tsv"
cazy_subfam_EC <- "input_file/annotations/CAZY_db/subfam_EC_mapping.tsv"

waafle_file <- "summary_tables/all_samples_internal_merged_indexed_filtered.tsv"
gene_unknown <- "input_file/annotations/geneID_map_unknown_only.tsv"
gene_map <- "input_file/annotations/geneID_map.tsv"

# load data
annot_df <- read.delim(file.path(basedir, annotation_file), quote = "") 
gene_map <- read.delim(file.path(basedir, gene_map))


####
cazy_ref <- read.delim(file.path(basedir, cazy_info)) %>%
  group_by(Family) %>%
  summarise(across(
    everything(),
    ~ {
      valid <- .[!is.na(.) & trimws(.) != ""]
      str_c(unique(valid), collapse = ",")
    }
  ), .groups = "drop") %>%
  distinct() %>%
  mutate(EC_Number=gsub(",","|",EC_Number))

cazy_subfam_info <- read.delim(file.path(basedir, cazy_subfam_info), header = F) %>%
  distinct() %>%
  mutate(V2_split = strsplit(V2, split = "\\|")) %>%
  unnest(V2_split) %>%
  group_by(V1) %>%
  summarise(
    V2_unique = str_c(unique(V2_split), collapse = "|"),
    .groups = "drop"
  ) %>%
  select(V1, V2_unique) %>%
  distinct() %>%
  rename(Subfamily=V1, Elements=V2_unique) %>%
  mutate(Family = str_split_i(Subfamily, pattern = "_", i = 1))

cazy_subfam_to_EC <- read.delim(file.path(basedir, cazy_subfam_EC), header=F)
colnames(cazy_subfam_to_EC) <- c("Subfamily", "EC_Number","n","Accessions","substrate")

# map CAZY to genes
cazy_map <- annot_df %>%
  mutate(ANNOTATIONS.UNIREF90=ANNOTATIONS.UNIREF90_fixed) %>%
  select(ANNOTATIONS.UNIREF90, CAZY_Domain, CAZY_validity) %>%
  filter(!is.na(CAZY_Domain)) %>%
  separate_longer_delim(cols = CAZY_Domain, delim = ",") %>%
  mutate(Family = str_split_i(CAZY_Domain, pattern = "_", i = 1)) %>%
  distinct() %>%
  left_join(cazy_ref, by = "Family") %>%
  mutate(
    EC_Number = ifelse(grepl("^\\d", CAZY_Domain), CAZY_Domain, EC_Number),
    EC_Number = ifelse(is.na(EC_Number) | trimws(EC_Number) == "", NA_character_, EC_Number)
  ) %>%
  rename(Substrate = Substrate_curated) %>%
  mutate(Substrate = ifelse(is.na(Substrate) | trimws(Substrate) == "", NA_character_, Substrate)) %>%
  mutate(Class = gsub("[0-9]", "", Family))
  
missing_substrate <- cazy_map %>%
  filter((is.na(Substrate_high_level)|is.na(Substrate)) & is.na(EC_Number)) %>% 
  select(CAZY_Domain, Family) %>%
  separate_longer_delim(cols = c(CAZY_Domain, Family), delim = "; ") %>%
  mutate(Class = gsub("[0-9]", "", Family)) %>%
  mutate(ID = as.numeric(gsub("[^0-9]", "", Family))) %>%
  filter(ID>0) %>%
  distinct() %>%
  arrange(Class, ID)

missing_EC <- cazy_map %>%
  filter(is.na(EC_Number)) %>% 
  filter(!CAZY_Domain %in% missing_substrate$CAZY_Domain) %>%
  select(CAZY_Domain, Family) %>%
  separate_longer_delim(cols = c(CAZY_Domain, Family), delim = "; ") %>%
  mutate(Class = gsub("[0-9]", "", Family)) %>%
  mutate(ID = as.numeric(gsub("[^0-9]", "", Family))) %>%
  filter(ID>0) %>%
  distinct() %>%
  arrange(Class, ID)

missing_cazy <- union(missing_substrate$CAZY_Domain, missing_EC$CAZY_Domain)

match_subfam <- cazy_subfam_info %>%
  select(Subfamily, Elements) %>%
  distinct() %>%
  separate_longer_delim(Elements, delim = "|") %>%
  distinct() %>%
  mutate(
    Family = str_split_i(Subfamily, pattern = "_", i = 1)) %>%
  filter(Elements %in% missing_cazy) %>%
  left_join(cazy_ref, by="Family") %>%
  select(-Subfamily) %>%
  distinct() %>%
  filter(Elements == Family) %>%
  group_by(Elements) %>%
  summarise(
    Family = {
      all_fam <- unlist(str_split(Family, "\\|"))
      v <- all_fam[!is.na(all_fam) & all_fam != "" & all_fam != "-"]
      str_c(unique(v), collapse = "; ")
    },
    Substrate_high_level = {
      all_sub <- unlist(str_split(Substrate_high_level, ";|,|\\|"))
      v <- all_sub[!is.na(all_sub) & all_sub != "" & all_sub != "-"]
      str_c(unique(v), collapse = ",")
    },
    Substrate = {
      all_sub <- unlist(str_split(Substrate_curated, ";|,|\\|"))
      v <- all_sub[!is.na(all_sub) & all_sub != "" & all_sub != "-"]
      str_c(unique(v), collapse = ",")
    },
    Name = {
      v <- Name[!is.na(Name) & Name != "" & Name != "-"]
      str_c(unique(v), collapse = "; ")
    },
    EC_Number = {
      all_ec <- unlist(str_split(EC_Number, "\\|"))
      v <- all_ec[!is.na(all_ec) & all_ec != "" & all_ec != "-"]
      str_c(unique(v), collapse = "; ")
    },
    .groups = "drop"
  ) %>%
  filter(!if_all(
    c(Substrate_high_level, Substrate, Name, EC_Number), 
    ~ is.na(.) | . == ""
  )) %>%
  select(-Elements,-EC_Number)

where_family_is_EC <- cazy_map %>%
  filter(EC_Number == CAZY_Domain & (Substrate_high_level == "" | is.na(Substrate_high_level))) %>%
  pull(EC_Number) %>%
  unique()

match_EC <- cazy_subfam_to_EC %>%
  mutate(Family = str_split_i(Subfamily, pattern = "_", i = 1)) %>%
  filter((Family %in% missing_cazy | EC_Number %in% where_family_is_EC) & !Family %in% match_subfam$Family) %>%
  group_by(Family) %>%
  summarise(
    Substrate = {
      all_sub <- unlist(str_split(substrate, ",|\\|"))
      v <- all_sub[!is.na(all_sub) & all_sub != "" & all_sub != "-"]
      str_c(unique(v), collapse = ",")
    },
    EC_Number = {
      all_ec <- unlist(str_split(EC_Number, "\\|"))
      v <- all_ec[!is.na(all_ec) & all_ec != "" & all_ec != "-"]
      str_c(unique(v), collapse = "|")
    },
    .groups = "drop"
  ) %>%
  mutate(Substrate = ifelse(Substrate == "", NA_character_, Substrate)) 

cazy_map_missing_subfam <- cazy_map %>%
  filter(CAZY_Domain %in% match_subfam$Family & is.na(EC_Number)) %>%
  select(-Substrate_high_level, -Substrate, -Name, -EC_Number) %>%
  left_join(match_subfam, by = "Family")

cazy_map <- cazy_map %>%
  filter(!(CAZY_Domain %in% match_subfam$Family & is.na(EC_Number))) %>%
  bind_rows(cazy_map_missing_subfam)

cazy_map_missing_ec <- cazy_map %>%
  filter((Family %in% missing_cazy | EC_Number %in% where_family_is_EC) & !Family %in% match_subfam$Family) %>%
  select(-EC_Number, -Substrate) %>%
  left_join(match_EC, by = "Family")

cazy_map <- cazy_map %>%
  filter(!((Family %in% missing_cazy | EC_Number %in% where_family_is_EC) & !Family %in% match_subfam$Family)) %>%
  bind_rows(cazy_map_missing_ec)


cazy_map_missing <- cazy_map %>%
  filter(if_all(c(Substrate_high_level, Substrate, EC_Number), ~ is.na(.x))) %>%
  select(-Name)

cazy_class_lookup <- tibble(
  Class = c("GH", "GT", "PL", "CE", "AA", "CBM"),
  CAZY_Class_Name = c(
    "glycoside hydrolases",
    "glycosyl transferases",
    "polysaccharide lyases",
    "carbohydrate esterases",
    "auxiliary activities",
    "carbohydrate-binding modules"
  )
)

cazy_map_final <- cazy_map %>%
  mutate(Class = na_if(Class, "...")) %>% 
  left_join(cazy_class_lookup, by="Class", keep=FALSE) %>%
  select(-Name)
  
get_KO <- annot_df %>%
  mutate(ANNOTATIONS.UNIREF90 = ANNOTATIONS.UNIREF90_fixed) %>%
  filter(ANNOTATIONS.UNIREF90 %in% cazy_map_final$ANNOTATIONS.UNIREF90) %>%
  select(ANNOTATIONS.UNIREF90, KO_ID, Consolidated_Function) %>%
  distinct() %>%
  group_by(ANNOTATIONS.UNIREF90) %>%
  summarise(
    KO_ID = {
      all <- unlist(str_split(KO_ID, ",|;|\\|"))
      v <- all[!is.na(all) & all != "" & all != "-"]
      str_c(unique(v), collapse = "; ")
    },
    Consolidated_Function = {
      all <- unlist(str_split(Consolidated_Function, ",|;|\\|"))
      v <- all[!is.na(all) & all != "" & all != "-"]
      str_c(unique(v), collapse = "; ")
    },
    .groups = "drop"
  )

cazy_map_final <- cazy_map_final %>%
  left_join(get_KO, by = "ANNOTATIONS.UNIREF90", keep = FALSE) %>%
  mutate(KO_ID=ifelse(KO_ID=="",NA,KO_ID)) %>%
  mutate(
    Consolidated_Function_new = case_when(
      !is.na(KO_ID) ~ Consolidated_Function,
      TRUE ~ paste(CAZY_Class_Name, Substrate_high_level, paste("(", EC_Number, ")", sep = ""))
    )
  ) %>%
  select(-Consolidated_Function) %>%
  rename(Consolidated_Function = Consolidated_Function_new) %>%
  distinct() %>%
  arrange(KO_ID)

# MAP AMR to Genes
amr_map <- annot_df %>%
  mutate(ANNOTATIONS.UNIREF90=ANNOTATIONS.UNIREF90_fixed) %>%
  select(ANNOTATIONS.UNIREF90,AMR_Element.symbol, AMR_Class, AMR_Type, KO_ID, Consolidated_Function) %>%
  filter(!is.na(AMR_Type)) %>%
  distinct() %>%
  group_by(ANNOTATIONS.UNIREF90) %>%
  summarise(across(
    everything(),
    ~ {
      valid <- .[!is.na(.) & trimws(.) != ""]
      str_c(unique(valid), collapse = "; ")
    }
  ), 
  .groups = "drop") %>%
  mutate(KO_ID=ifelse(KO_ID=="",NA,KO_ID)) %>%
  mutate(
    Consolidated_Function_new = case_when(
      !is.na(KO_ID) ~ Consolidated_Function,
      TRUE ~ paste(AMR_Type, AMR_Class, paste("(", AMR_Element.symbol, ")", sep = ""))
    )
  ) %>%
  select(-Consolidated_Function) %>%
  rename(Consolidated_Function = Consolidated_Function_new) %>%
  distinct() %>%
  arrange(KO_ID)
  

other_genes <- annot_df %>%
  mutate(ANNOTATIONS.UNIREF90=ANNOTATIONS.UNIREF90_fixed) %>%
  filter(!(ANNOTATIONS.UNIREF90 %in% amr_map$ANNOTATIONS.UNIREF90 |
           ANNOTATIONS.UNIREF90 %in% cazy_map_final$ANNOTATIONS.UNIREF90)) %>%
  select(ANNOTATIONS.UNIREF90,KO_ID, Consolidated_Function) %>%
  distinct() %>%
  group_by(ANNOTATIONS.UNIREF90) %>%
  summarise(across(
    everything(),
    ~ {
      valid <- .[!is.na(.) & trimws(.) != ""]
      str_c(unique(valid), collapse = "; ")
    }
  ), 
  .groups = "drop") %>%
  mutate(is_amr=FALSE, 
         what_amr=NA_character_,
         is_cazy=FALSE,
         what_cazy=NA_character_) %>%
  arrange(KO_ID)

consolidated_map <- other_genes %>%
  bind_rows(amr_map %>%
              mutate(is_amr=TRUE, 
                     what_amr=AMR_Class,
                     is_cazy=FALSE,
                     what_cazy=NA_character_) %>%
              select(all_of(colnames(other_genes)))) %>%
  bind_rows(cazy_map_final %>%
              select(ANNOTATIONS.UNIREF90, Family, KO_ID, Consolidated_Function) %>%
              group_by(ANNOTATIONS.UNIREF90) %>%
              summarise(across(
                everything(),
                ~ {
                  valid <- .[!is.na(.) & trimws(.) != ""]
                  str_c(unique(valid), collapse = "; ")
                }
              ), 
              .groups = "drop") %>%
              mutate(is_amr=FALSE, 
                     what_amr=NA_character_,
                     is_cazy=TRUE,
                     what_cazy=Family) %>%
              select(all_of(colnames(other_genes)))) %>%
  mutate(KO_ID=ifelse(KO_ID=="",NA,KO_ID)) %>%
  arrange(ANNOTATIONS.UNIREF90)

consolidated_map %>% filter(duplicated(ANNOTATIONS.UNIREF90))

# deduplicate
# if ANNOTATIONS.UNIREF90 is duplicated, keep only the line with the least element in KO_ID
consolidated_map <- consolidated_map %>%
  group_by(ANNOTATIONS.UNIREF90) %>%
  mutate(
    min_ko = map_chr(str_split(KO_ID, "; "), ~ {
      valid_elements <- .x[!is.na(.x) & .x != ""]
      if (length(valid_elements) == 0) return(NA_character_)
      min(valid_elements)
    })
  ) %>%
  arrange(min_ko, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-min_ko)

# validate
length(unique(consolidated_map$ANNOTATIONS.UNIREF90)) == length(unique(annot_df$ANNOTATIONS.UNIREF90_fixed))

# writeout
gene_annotation_mapping <- list(
  cazy=cazy_map_final,
  amr=amr_map)

saveRDS(gene_annotation_mapping,"input_file/annotations/gene_annotation_mapping.rds")
saveRDS(consolidated_map,"input_file/annotations/consolidated_gene_annotation_mapping.rds")

