type_col <- site_col[1:3]
names(type_col) <- c("Fecal","Skin","Oral")

hh_dogs <- read.delim(file.path(input_dir,"Household_counts_filledv2.tsv")) %>%
  filter(!is.na(sampleid))

sid_with_dogs <- hh_dogs$sampleid[which(hh_dogs$Dogs_HH > 0)]

library(dplyr)
library(stringr) # or stringi

dog_socials <- read.delim("input_file/socials/dog-owner_pairs.tsv") %>%
  mutate(owner_id = gsub("-","",owner_id),
         dog_id = gsub("-","", dog_id)) %>%
  separate_longer_delim(owner_id, delim = ";") %>%
  separate_longer_delim(dog_id, delim = ";") %>%
  filter(dog_id != "" & owner_id != "") %>%
  mutate(s1=as.character(owner_id), s2=as.character(dog_id)) %>%
  mutate(sample1_2 = paste(pmin(s1, s2), pmax(s1, s2), sep = "_"))

hh_pairs <- unique(dog_socials$sample1_2)
length(hh_pairs)

dog_human_pairs <- sLGT_long %>%
  filter(host_pair=="dog_human") %>%
  mutate(type =  case_when(
    (grepl("DDF|DHF", substr(seqID_1, 1,3))) & (grepl("DDF|DHF", substr(seqID_2, 1,3))) ~ "Fecal",
    (grepl("DDV|DHV", substr(seqID_1, 1,3))) & (grepl("DDV|DHV", substr(seqID_2, 1,3))) ~ "Oral",
    (grepl("DDH|DHP", substr(seqID_1, 1,3))) & (grepl("DDH|DHP", substr(seqID_2, 1,3))) ~ "Skin",
    TRUE ~ "MIXED")) %>%
  filter(type != "MIXED") %>%
  mutate(with_dogs = ifelse(
    ((sample_1 %in% sid_with_dogs) | (sample_2 %in% sid_with_dogs)),
    "with_dogs", "no_dogs")) %>%
  mutate(
    s1 = as.character(sub("_.*$", "", sample1_2)),
    s2 = as.character(sub("^.*_", "", sample1_2)),
    sample1_2_rev = paste(pmin(s1, s2), pmax(s1, s2), sep = "_")
  ) %>%
  mutate(pop1 = substr(sample_1,1,3),
         pop2 = substr(sample_2,1,3)) %>%
  mutate(relation = case_when(
    sample1_2_rev %in% hh_pairs ~ "same_household",
    sample1_2 %in% hh_pairs ~ "same_household",
    pop1 == pop2 ~ "same_village",
    TRUE ~ "different_village"
  )) %>%
  mutate(relation_withDogs = paste(relation, with_dogs, sep="_")) %>%
  select(-s1, -s2, -pop1, -pop2)

table(dog_human_pairs$type)
table(dog_human_pairs$relation_withDogs)
table(dog_human_pairs$sample1_2_rev %in% hh_pairs)
table(dog_human_pairs$sample1_2 %in% hh_pairs)
table(dog_human_pairs$relation == "same_village")


dog_human_pairs <- dog_human_pairs %>%
  mutate(type=factor(type, c("Fecal","Oral","Skin")),
         with_dogs=factor(with_dogs, c("with_dogs","no_dogs")),
         relation=factor(relation, c("different_village","same_village","same_household"))) %>%
  arrange(relation,with_dogs, type) %>%
  mutate(relation_withDogs=factor(relation_withDogs, unique(relation_withDogs)))

relation_levels <- levels(dog_human_pairs$relation_withDogs)

raw_plot_list <- dog_human_pairs %>%
  group_split(type) %>%
  map(function(sub_df) {
    current_type <- as.character(sub_df$type[1])
    current_host <- "human-dog pairs"
    
    setlist <- union(sub_df$seqID_1, sub_df$seqID_2)
    n_unique_human <- sum(grepl("DH",setlist))
    n_unique_dog <- sum(grepl("DD",setlist))
    title_text <- paste0(current_type, " (human = ", n_unique_human, ", dog = ", n_unique_dog, ")")
    
    # Ensure factor ordering (1: different_village, 2: same_village, 3: same_household)
    sub_df$relation <- droplevels(factor(sub_df$relation_withDogs, levels = relation_levels))
    nPairs <- table(sub_df$relation_withDogs)
    
    iLevels <- levels(sub_df$relation_withDogs)
    
    sub_df <- sub_df %>% arrange(relation_withDogs)
    
    # ----------------------------------------------------
    # 1. PAIRWISE WILCOXON TESTS
    # ----------------------------------------------------
    # Filter non-empty groups with sufficient data
    valid_groups <- sub_df %>%
      group_by(relation_withDogs) %>%
      summarise(n = sum(!is.na(LGT_similarity)), .groups = "drop") %>%
      filter(n > 0) %>%
      pull(relation_withDogs)
    
    df_brackets <- tibble()
    
    if (length(valid_groups) >= 2 && length(unique(na.omit(sub_df$LGT_similarity))) >= 2) {
      
      # Pairwise Wilcoxon rank-sum test with p-adjustment
      pw_res <- pairwise.wilcox.test(
        x = sub_df$LGT_similarity,
        g = sub_df$relation_withDogs,
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
            x1 = as.numeric(factor(group1, levels = iLevels)),
            x2 = as.numeric(factor(group2, levels = iLevels)),
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
      group_by(relation_withDogs) %>%
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
    p <- ggplot(sub_df, aes(x = relation_withDogs, y = LGT_similarity, fill = type)) +
      stat_summary(fun = "mean", geom = "bar", width = 0.9) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, colour = "black") +
      scale_fill_manual(values = type_col) +
      scale_x_discrete(
        expand = expansion(mult = c(0.1, 0.1)),
        labels = function(x) {
          paste0(x, "\n(nPairs = ", as.numeric(nPairs[x]), ")")
        }
      ) +
      labs(title = title_text, x="LGT Similarity") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_line(colour = "grey80", linewidth = 0.5),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 9.5, colour = "grey40")
      )
    
    # Add dynamic brackets if comparisons exist
    df_brackets <- df_brackets %>%
      filter(p_stars != "ns")
    
    if (!is.null(df_brackets) & nrow(df_brackets) > 0) {
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

raw_plot_list