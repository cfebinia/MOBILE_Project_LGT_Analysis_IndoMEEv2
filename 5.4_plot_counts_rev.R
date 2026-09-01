rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(patchwork)
library(purrr)

# Paths
# my_lib <- "/home/caf77/R/x86_64-pc-linux-gnu-library/4.2"
# .libPaths(c(my_lib, .libPaths()))

# basedir="C:/Users/caf77_Local/Documents/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
basedir="C:/Users/caf77/OneDrive - University of Cambridge/Imaging_Lab_PC1/WAAFLE_Extra"
setwd(basedir)

input_dir=file.path(basedir,"input_file")
index_dir <- file.path(basedir, "index")

summary_dir <- file.path(basedir, "summary_tables")

figpath <- file.path(basedir, "figures")
dir.create(figpath,showWarnings = F,recursive = T)

waffle_file="all_samples_internal_merged_indexed_filtered.tsv"
LGT_count="nLGT_by_sample.tsv" #"nLGT_by_seqID.tsv"
outliers_list="outliers_waafle.txt"
alpha_rds <- "adiv_metrics.rds"
n_contigs <- "input_file/input_contigs_compile.tsv"

site_col <- c(
  "humanFecal"="#56B4E9", 
  "humanSkin"="#EE2C2C", 
  "humanOral"="#ef9f00", 
  "dogFecal"="#2850c8", 
  "dogSkin"="#780247", 
  "dogOral"="#9A5324"  
)

pop.ord <- c("BTU", "ORT", "ORS", "APT", "TBU", "BSP", "LDY")
type.ord <- c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")

pLGT <- "pLGT_bySample_withRep.tsv"
nContigs <- "nContig_bySample_withRep.tsv"

#######################
## Helper Functions  ##
#######################
save_plot <- function(p, name, width, height, unit="in", scale=1, dpi = 600) {
  formats <- c("rds", "png", "pdf")
  for(fmt in formats){
    file <- file.path(figpath, paste0(name, ".", fmt))
    if(fmt == "rds") saveRDS(p, file)
    else ggsave(file, plot = p, width = width, height = height, dpi = dpi, unit=unit, scale=scale) 
  }
}

################
## Load Data  ##
################
outliers <- read.delim(file.path(input_dir,"outliers_waafle.txt"))

pop.df <- read.delim(file.path(basedir, "input_file/Consolidated_SampleID_combined_final_SeqIDmatched.tsv")) %>%
  mutate(SAMPLE = substr(seqID,1,6)) %>%
  select(SAMPLE, population, individual) %>%
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


#total_contigs <- read.delim(n_contigs) %>%
#  rename(SAMPLE=file_name) %>%
#  filter(!SAMPLE %in% outliers$sampleid) %>%
#  mutate(type=substr(SAMPLE,1,3)) %>%
#  mutate(type=case_when(
#    type=="DHF" ~ "humanFecal",
#    type=="DHV" ~ "humanOral",
#    type=="DHP" ~ "humanSkin",
#    type=="DDF" ~ "dogFecal",
#    type=="DDV" ~ "dogOral",
#    type=="DDH" ~ "dogSkin"
#  )) 

#df <- read.delim(file.path(summary_dir,LGT_count)) %>%
#  filter(!SAMPLE %in% outliers$sampleid) %>%
#  mutate(type=case_when(
#    type=="DHF" ~ "humanFecal",
#    type=="DHV" ~ "humanOral",
#    type=="DHP" ~ "humanSkin",
#    type=="DDF" ~ "dogFecal",
#    type=="DDV" ~ "dogOral",
#    type=="DDH" ~ "dogSkin"
#  )) %>%
#  mutate(population=factor(population, pop.ord),
 #        group=factor(group,c("Early-transition","Late-transition","Agriculture")),
 #        type=factor(type, type.ord),
  #       host=factor(host, c("human","dog")),
   #      site=factor(site, c("Fecal","Oral","Skin"))) 


df <- read.delim(file.path(summary_dir,pLGT)) %>%
  select(original_sid, count=raw_count_LGT, nSandWichDist, n_contig, pLGT) %>%
  filter(!original_sid %in% outliers$sampleid) %>%
  mutate(SAMPLE=substr(original_sid, 1,6),
         type=substr(SAMPLE,1,3)) %>%
  mutate(type=case_when(
    type=="DHF" ~ "humanFecal",
    type=="DHV" ~ "humanOral",
    type=="DHP" ~ "humanSkin",
    type=="DDF" ~ "dogFecal",
    type=="DDV" ~ "dogOral",
    type=="DDH" ~ "dogSkin"
  )) %>%
  mutate(site = gsub("human|dog", "", type),
         host = ifelse(grepl("human", type), "human", "dog")) %>%
  left_join(pop.df %>% select(SAMPLE, individual, population, group), by="SAMPLE") %>%
  mutate(population=factor(population, pop.ord),
         group=factor(group,c("Early-transition","Late-transition","Agriculture")),
         type=factor(type, type.ord),
         host=factor(host, c("human","dog")),
         site=factor(site, c("Fecal","Oral","Skin"))) 

total_pLGT <- df %>%
  rename(count_unique=count) %>%
  group_by(type) %>%
  summarise(count_unique=sum(count_unique), 
            total_nContig=sum(n_contig),
            .groups = "drop") %>%
  mutate(pLGT=(count_unique/total_nContig)*100)

df %>% 
  group_by(type) %>%
  summarise(nSample=n_distinct(SAMPLE))

total_contigs <-  df %>%
  select(original_sid, SAMPLE, n_contig, type)

#####################
## Summary Tables  ##
#####################
get_meanSE <- function(x, digits = 1) {
  x <- na.omit(x)
  i <- round(mean(x), digits)
  j <- round(sd(x) / sqrt(length(x)), digits)
  k <- paste(i, "±", j, sep = "")
  return(k)
}

get_range <- function(x, digits = 1) {
  x <- na.omit(x)
  i <- round(quantile(x, c(0, 1)), digits)
  k <- paste(i[1], "--", i[2], sep = "")
  return(k)
}

############
## Plots  ##
############
size_annotate <- df %>%
  group_by(type) %>%
  mutate(x_lab = max(pLGT)*1.15) %>%
  group_by(population, type) %>%
  summarise(
    n = paste("(n=", n_distinct(SAMPLE), ")", sep = ""),
    x_lab = first(x_lab),
    .groups = "drop"
  )

plot_nLGT <- df %>%
  ggplot(aes(x = pLGT, y = population, colour = type, fill = type)) +
  geom_density_ridges(alpha = 0.3) +
  facet_grid(~type, scales = "free") +
  theme_minimal(base_size = 20) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(),
    axis.text.x = element_text(size=16)
  ) +
  geom_text(
    data = size_annotate,
    aes(x = x_lab, label = n, y=population),
    inherit.aes = FALSE,
    vjust = 1,
    hjust = 1, 
    size = 5
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = seq(0, 100, by = 0.2)) +
  scale_fill_manual(values = site_col) +
  scale_colour_manual(values = site_col) +
  labs(
    y = "Community",
    x = "pLGT (%)"
  )

plot_nLGT

# boxplot
counts <- df %>% 
  group_by(type) %>%
  summarise(nSample=n_distinct(SAMPLE))
type_labels <- setNames(paste0(counts$type, "\n(", counts$nSample, ")"), counts$type)

my_comparisons <- list(c("dogFecal", "humanFecal"),
                       c("dogOral","humanOral"),
                       c("dogSkin","humanSkin"))



nLGT_box <- ggplot(df, 
                   aes(x = type, y = pLGT, fill = type)) +
  geom_point(
    aes(fill = type, colour = type), shape = 21, alpha = 0.5, size = 2.5,
    position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
  ) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_point(data=total_pLGT, aes(x=type, y=pLGT), col="black", shape=16, size=5)+
  scale_x_discrete(labels = type_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  scale_fill_manual(values = site_col) +
  scale_colour_manual(values = site_col) +
  theme_minimal(base_size = 20) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    legend.title = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    x = "Sample Type",
    y = "pLGT (%)"
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 8
  )

nLGT_box


comp_1 <- nLGT_box / plot_nLGT + 
  plot_layout(heights = c(0.6, 1)) + 
  plot_annotation(tag_levels = "A") & 
  theme(
    axis.text = element_text(size = 18),
    plot.tag = element_text(size = 35, face = "bold"),
    legend.title = element_text()
    )

plot(comp_1)

save_plot(comp_1, "nLGT_composite", width = 15, height = 18, scale=1)


######

my_comparisons <- list(
  c("Late-transition", "Early-transition"),
  c("Early-transition", "Agriculture")
)

plot_list <- df %>%
  split(list(.$host, .$site), drop = TRUE) %>%
  map(function(sub_df) {
    
    # Extract titles based on the current subset
    current_host <- unique(sub_df$host)
    current_site <- unique(sub_df$site)
    
    counts <- sub_df %>% 
      group_by(group) %>%
      summarise(nSample=n_distinct(SAMPLE))
    type_labels <- setNames(paste0(counts$group, "\n(", counts$nSample, ")"), counts$type)
    
    my_comparisons <- combn(as.character(unique(sub_df$group)),2,simplify = F)
    
    ggplot(sub_df, aes(x = group, y = pLGT, fill = type)) +
      geom_point(
        aes(colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, fill=NA) +
      scale_x_discrete(labels = type_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1, size=14),
        axis.title.x = element_blank(),
        plot.title = element_text(face="bold", size=14)
      ) +
      labs(title = paste(current_host, current_site, sep = ""), y="pLGT (%)")+
      stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        size = 5
      )
  })

combined_plot <- wrap_plots(plot_list[c(1,3,5,2,4,6)]) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    #title = "# LGT event / # Total Contig",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_Ls", height = 10, width = 15, unit = "cm", scale=2.2)

#######################
# Assoc with richness #
#######################
library(psych)

adf <- readRDS(file.path(input_dir, alpha_rds)) %>%
  select(-faith) %>%
  mutate(SAMPLE = substr(sampleid, 1, 6)) %>%
  filter(SAMPLE %in% unique(df$SAMPLE)) %>%
  select(-sampleid) %>%
  group_by(SAMPLE) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")


corr_adiv <- apply(adf[, -1], 2, function(metric) {
  cordf <- data.frame(SAMPLE = adf$SAMPLE, adiv = metric) %>%
    left_join(df %>%
                select(SAMPLE, pLGT, n_contig, type) %>%
                group_by(type, SAMPLE) %>%
                summarise(pLGT = mean(pLGT, na.rm=T),
                          n_contig = mean(n_contig, na.rm=T),
                          .groups = "drop"),
              by = "SAMPLE") 
  
  res <- corr.test(x = cordf$adiv, y = cordf$pLGT, method = "spearman", adjust = "none")
  
  # Extract the correlation coefficient and p-value
  c(r = as.numeric(res$r), p = as.numeric(res$p))
})

corr_adiv <- data.frame(t(corr_adiv))
corr_adiv$p.adj <- p.adjust(corr_adiv$p)


df_long <-  adf %>%
  tidyr::pivot_longer(
    cols = -SAMPLE, 
    names_to = "metric_name", 
    values_to = "adiv_value"
  ) %>%
  left_join(df %>%
               select(SAMPLE, pLGT, n_contig, type) %>%
               group_by(type, SAMPLE) %>%
               summarise(pLGT = mean(pLGT, na.rm=T),
                         n_contig = mean(n_contig, na.rm=T),
                         .groups = "drop"),
                 by = "SAMPLE") %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  mutate(metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance")))

corr_raw <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  split(list(.$metric_name), drop = TRUE) %>%
  map_dfr(function(sub_df) {
    if (nrow(sub_df) < 3) return(NULL)
    res <- corr.test(x = sub_df$adiv_value, y = sub_df$pLGT, method = "spearman", adjust = "none")
    data.frame(
      metric_name = unique(sub_df$metric_name),
      r = as.numeric(res$r),
      p = as.numeric(res$p),
      stringsAsFactors = FALSE
    )
  })

x_positions <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  group_by(metric_name) %>%
  summarise(x_pos = min(adiv_value, na.rm = TRUE) * 1.05, .groups = "drop")

corr_annotations <- corr_raw %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p_stars = as.character(symnum(p.adj, 
                                  cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                  symbols = c("***", "**", "*", "ns"))),
    label = paste0("rho == ", round(r, 2), "~'", p_stars, "'"),
    metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance"))
  ) %>%
  left_join(x_positions, by = "metric_name")

plot_cor <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  ggplot(aes(x = adiv_value, y = pLGT, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 1.2) +
  geom_text(
    data = corr_annotations,
    aes(x = x_pos, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.5,
    hjust = 0,
    size = 4.5,
    parse = TRUE
  ) +
  facet_grid(. ~ metric_name, scales = "free") +
  theme_minimal(base_size = 20) +
  scale_color_manual(values = site_col) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(t = 2, b = 2, r = 2, l = 2), size = 14, face="bold"),
    axis.title = element_text(size=14, face="bold")
  ) +
  labs(
    x = "Alpha Diversity",
    y = "pLGT (%)"
  )

plot_cor

save_plot(plot_cor, "nLGT_vs_Richness_ALL", height = 5, width = 15, unit = "cm", scale=1.8)



#By Sample Type
df_long <- adf %>%
  tidyr::pivot_longer(
    cols = -SAMPLE, 
    names_to = "metric_name", 
    values_to = "adiv_value"
  ) %>%
  left_join(df %>%
              select(SAMPLE, pLGT, n_contig, type) %>%
              group_by(type, SAMPLE) %>%
              summarise(pLGT = mean(pLGT, na.rm=T),
                        n_contig = mean(n_contig, na.rm=T),
                        .groups = "drop"),
            by = "SAMPLE")

df_long <- read.delim(file.path(input_dir, "allUnified_draft3_weakFilter_mapStats.tsv"), stringsAsFactors = FALSE) %>%
  rename(nReads=bwa_counts_total_pass) %>%
  mutate(SAMPLE=substr(sample,1,6)) %>%
  select(SAMPLE, nReads) %>%
  filter(SAMPLE %in% df_long$SAMPLE) %>%
  group_by(SAMPLE) %>%
  summarise(nReads=mean(nReads)/10^6, .groups = "drop") %>%
  right_join(df_long,
            by = "SAMPLE")

corr_adiv_by_type <- df_long %>%
  split(list(.$type, .$metric_name), drop = TRUE) %>%
  map_dfr(function(sub_df) {
    
    # Check for sufficient data points to run a correlation
    if (nrow(sub_df) < 3) return(NULL)
    
    res <- corr.test(
      x = sub_df$adiv_value, 
      y = sub_df$pLGT, 
      method = "spearman", 
      adjust = "none"
    )
    
    data.frame(
      type = unique(sub_df$type),
      metric = unique(sub_df$metric_name),
      r = as.numeric(res$r),
      p = as.numeric(res$p),
      stringsAsFactors = FALSE
    )
  })

corr_adiv_by_type$p.adj <- p.adjust(corr_adiv_by_type$p)
corr_adiv_by_type[corr_adiv_by_type$p.adj<0.05,] %>% arrange(type)


df_long <- df_long %>%
  mutate(metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance")))

corr_raw <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  split(list(.$type, .$metric_name), drop = TRUE) %>%
  map_dfr(function(sub_df) {
    if (nrow(sub_df) < 3) return(NULL)
    res <- corr.test(x = sub_df$adiv_value, y = sub_df$pLGT, method = "spearman", adjust = "none")
    data.frame(
      type = unique(sub_df$type),
      metric_name = unique(sub_df$metric_name),
      r = as.numeric(res$r),
      p = as.numeric(res$p),
      stringsAsFactors = FALSE
    )
  })

x_positions <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  group_by(metric_name) %>%
  summarise(x_pos = min(adiv_value, na.rm = TRUE) * 1.05, .groups = "drop")

xy_positions <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  group_by(type, metric_name) %>%
  summarise(
    x_pos = {
      center_vals <- adiv_value[adiv_value < quantile(adiv_value, 0.25, na.rm = TRUE) | 
                                  adiv_value > quantile(adiv_value, 0.75, na.rm = TRUE)]
      vals <- table(round(center_vals, 0))
      min_count <- min(vals)
      as.numeric(names(vals[vals == min_count])[1])
    },
    y_pos = {
      center_vals <- pLGT[pLGT < quantile(pLGT, 0.25, na.rm = TRUE) | 
                            pLGT > quantile(pLGT, 0.75, na.rm = TRUE)]
      vals <- table(round(center_vals, 1))
      min_count <- min(vals)
      as.numeric(names(vals[vals == min_count])[1])
    },
    hjust = ifelse(x_pos < median(adiv_value, na.rm = TRUE), 0, 1),
    vjust = ifelse(y_pos < median(pLGT, na.rm = TRUE), 0, 1),
    x_pos = ifelse(x_pos < median(x_pos, na.rm = TRUE), x_pos*0.8, x_pos*1.2),
    y_pos = ifelse(y_pos < median(pLGT, na.rm = TRUE), y_pos*0.8, y_pos*1.2),
    .groups = "drop"
  )

corr_annotations <- corr_raw %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p_stars = as.character(symnum(p.adj, 
                                  cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                  symbols = c("***", "**", "*", "ns"))),
    label = paste0("rho == ", sprintf("%.2f", r), " ~ '", p_stars, "'"),
    metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance"))
  ) %>%
  left_join(xy_positions, by = c("type", "metric_name"))

plot_cor <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  ggplot(aes(x = adiv_value, y = pLGT, colour = type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, col = "black", linewidth = 1.2) +
  geom_text(
    data = corr_annotations,
    aes(x = x_pos, y = y_pos, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = 4.5,
    parse = TRUE
  ) +
  facet_grid(type ~ metric_name, scales = "free") +
  theme_minimal(base_size = 20) +
  scale_color_manual(values = site_col) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(t = 2, b = 2, r = 2, l = 2), size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "Alpha Diversity",
    y = expression(paste(frac("#LGT events", "# Total Contigs"), " (%)"))
  )

plot_cor

save_plot(plot_cor, "nLGT_vs_Richness", height = 18, width = 15, unit = "cm", scale=1.5)


## lm
library(purrr)
library(MuMIn)

df_long <- df_long %>%
  mutate(pLGT_scaled = pLGT) %>%
  na.omit()
map_type <- c(
  "DHF" = "humanFecal",
  "DDF" = "dogFecal",
  "DHV" = "humanOral",
  "DDV" = "dogOral",
  "DHP" = "humanSkin",
  "DDH" = "dogSkin"
)

df_long$sType <- map_type[substr(df_long$SAMPLE,1,3)]
table(df_long$sType)

xy_positions <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  group_by(type, metric_name) %>%
  summarise(
    x_pos = {
      center_vals <- adiv_value[adiv_value < quantile(adiv_value, 0.1, na.rm = TRUE) | 
                                  adiv_value > quantile(adiv_value, 0.9, na.rm = TRUE)]
      vals <- table(round(center_vals, 0))
      min_count <- if (length(vals) > 0) min(vals) else 0
      if (length(vals) > 0) as.numeric(names(vals[vals == min_count])[1]) else mean(adiv_value, na.rm = TRUE)
    },
    y_pos = {
      center_vals <- pLGT_scaled[pLGT_scaled < quantile(pLGT_scaled, 0.1, na.rm = TRUE) | 
                                   pLGT_scaled > quantile(pLGT_scaled, 0.9, na.rm = TRUE)]
      vals <- table(round(center_vals, 1))
      min_count <- if (length(vals) > 0) min(vals) else 0
      if (length(vals) > 0) as.numeric(names(vals[vals == min_count])[1]) else mean(pLGT_scaled, na.rm = TRUE)
    },
    hjust = ifelse(x_pos < median(adiv_value, na.rm = TRUE), 0, 1),
    vjust = ifelse(y_pos < median(pLGT_scaled, na.rm = TRUE), 0, 1),
    x_pos = ifelse(x_pos < median(x_pos, na.rm = TRUE), x_pos * 1.1, x_pos * 0.8),
    y_pos = ifelse(y_pos < median(pLGT_scaled, na.rm = TRUE), y_pos * 1.1, y_pos * 0.8),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = ifelse(grepl("Fecal", type), y_pos + 0.15, y_pos),
    y_pos = ifelse(grepl("dogSkin", type), y_pos + 0.4, y_pos),
    x_pos = ifelse(grepl("dominance", metric_name) & grepl("Skin|humanOral", type), x_pos + 0.5, x_pos),
    x_pos = ifelse(grepl("richness", metric_name), x_pos + 600, x_pos),
    x_pos = ifelse(grepl("shannon", metric_name) & grepl("dogOral", type), x_pos + 1, x_pos),
    x_pos = ifelse(grepl("shannon", metric_name) & grepl("dogSkin", type), x_pos - 1, x_pos)
  )

model_outputs <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  mutate(metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance"))) %>%
  group_by(metric_name, type) %>%
  filter(n() >= 3) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~ lm(pLGT ~ scale(adiv_value) + log10(nReads), data = .x)),
    summary_mod = purrr::map(model, summary),
    
    grid_df = purrr::map2(data, model, function(sub_df, mod) {
      grid <- expand.grid(
        adiv_value = seq(min(sub_df$adiv_value, na.rm = TRUE), max(sub_df$adiv_value, na.rm = TRUE), length.out = 100),
        nReads = mean(sub_df$nReads, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      grid$pred <- predict(mod, newdata = grid)
      grid
    }),
    
    stats_df = purrr::map2(data, summary_mod, function(sub_df, summary_obj) {
      r_sq <- summary_obj$r.squared
      coef_matrix <- summary_obj$coefficients
      
      pred_row <- grep("adiv_value", rownames(coef_matrix), value = TRUE)[1]
      
      beta <- if (!is.na(pred_row)) coef_matrix[pred_row, "Estimate"] else NA
      p_val <- if (!is.na(pred_row)) coef_matrix[pred_row, "Pr(>|t|)"] else NA
      
      p_label <- ifelse(is.na(p_val), "NA",
                        ifelse(p_val < 0.001, "***",
                               ifelse(p_val < 0.01, "**",
                                      ifelse(p_val < 0.05, "*", "ns"))))
      
      data.frame(
        label = paste0("B = ", round(beta, 2), p_label, ", R2 = ", round(r_sq, 2)),
        stringsAsFactors = FALSE
      )
    })
  ) %>%
  ungroup() %>%
  left_join(xy_positions, by = c("metric_name", "type")) %>%
  mutate(
    stats_df = purrr::map2(stats_df, seq_len(n()), function(sdf, idx) {
      sdf$x_pos <- x_pos[idx]
      sdf$y_pos <- y_pos[idx]
      sdf$hjust <- hjust[idx]
      sdf$vjust <- vjust[idx]
      sdf
    })
  )

pred_df <- model_outputs %>%
  select(metric_name, type, grid_df) %>%
  tidyr::unnest(cols = c(grid_df))

stats_df <- model_outputs %>%
  select(metric_name, type, stats_df) %>%
  tidyr::unnest(cols = c(stats_df))

plot_lm <- df_long %>%
  filter(metric_name %in% c("richness", "shannon", "dominance")) %>%
  mutate(metric_name = factor(metric_name, levels = c("richness", "shannon", "dominance"))) %>%
  ggplot(aes(x = adiv_value, y = pLGT_scaled, colour = type)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_line(
    data = pred_df, 
    aes(y = pred), 
    linewidth = 1.2,
    colour = "black"
  ) +
  geom_text(
    data = stats_df,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    parse = FALSE,
    size = 4,
    hjust = 0.5,
    vjust = 1
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.1, 0.1)),
    breaks = function(x) {
      b <- scales::extended_breaks()(x)
      b[b >= 0]
    }
  ) +
  facet_grid(type ~ metric_name, scales = "free") +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = site_col) +
  labs(
    x = "Alpha Diversity",
    y = expression(paste(frac("#LGT events", "# Total Contigs"), " (%)")),
    caption = "Model: pLGT_scaled ~ adiv_value + log10(nReads)\n(fitted independently for each metric_name and type;\nlines represent predictions at mean nReads)"
  ) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(margin = margin(t = 2, b = 2, r = 2, l = 2), size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.caption = element_text(size = 14, face = "plain", hjust = 1)
  )

plot_lm
save_plot(plot_lm, "nLGT_vs_Richness_lmer_log10", height = 18, width = 15, unit = "cm", scale = 1.8)


#######################
# by Rel_young
######################
metadata <- read.delim(file.path(input_dir,"combined_ind_metadata.tsv")) %>%
  rename(sampleid=Sample) %>%
  select(sampleid, age, sex) %>%
  distinct()

hh_id <- read.delim("input_file/socials/addresses.txt")

get_mode <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NA)
  ux <- unique(x_clean)
  ux[which.max(tabulate(match(x_clean, ux)))]
}

hh_member <- read.delim(file.path(input_dir, "Household_counts_filledv2.tsv")) %>%
  filter(!is.na(Total_HH)) %>%
  left_join(hh_id, by = "sampleid") %>%
  mutate(Rel_young = Rel_children + Rel_childrenstep + Rel_neicenephewinlaw + Rel_grandchildren) %>%
  select(sampleid, population, HH_ID, Total_HH, Rel_young) %>%
  left_join(metadata, by = "sampleid") %>%
  group_by(population, HH_ID) %>%
  filter(Total_HH == get_mode(Total_HH)) %>%
  ungroup() %>%
  rename(individual=sampleid) %>%
  select(-population)


df_relYoung <- df %>%
  filter(host == "human") %>%
  filter(group != "Agriculture") %>%
  inner_join(hh_member, by="individual") %>%
  mutate(with_young = ifelse(Rel_young>0, "with_young", "no_young")) %>%
  mutate(group2 = paste(group, with_young, sep = "_")) %>%
  mutate(group = droplevels(group),
         population = droplevels(population),
         with_young = factor(with_young, c("with_young", "no_young")),
         SAMPLE = as.character(SAMPLE)) %>%
  arrange(group, population, with_young, SAMPLE) %>%
  mutate(group2 = factor(group2, unique(group2))) %>%
  na.omit()

# HH size
library(lmerTest)
library(MuMIn)

plot_list <- df_relYoung %>%
  #split(list(.$site, .$group)) %>%
  split(.$site) %>%
  map(function(sub_df) {
    #sub_df <- df_relYoung %>% filter(site=="Fecal")
    lower_bound <- quantile(sub_df$pLGT  , 0.005, na.rm = TRUE)
    upper_bound <- quantile(sub_df$pLGT  , 0.995, na.rm = TRUE)
    sub_df <- sub_df[sub_df$pLGT   >= lower_bound & sub_df$pLGT  <= upper_bound, ]
    
    n_ind <- length(unique(sub_df$individual))
    if (nrow(sub_df) == 0) return(NULL)
    
    current_host  <- as.character(unique(sub_df$host))[1]
    current_site  <- as.character(unique(sub_df$site))[1]
    
    fit <- lmer(pLGT ~ log10(Total_HH+0.1) + age + log10(n_contig) + (1|population), data = sub_df)
    fit_summary <- summary(fit)
    
    pred_grid <- data.frame(
      Total_HH  = seq(min(sub_df$Total_HH , na.rm = TRUE), max(sub_df$Total_HH , na.rm = TRUE), length.out = 100),
      age = mean(sub_df$age, na.rm=TRUE),
      n_contig = mean(sub_df$n_contig, na.rm = TRUE)
    )
    pred_grid$pred_pLGT <- predict(fit, newdata = pred_grid, re.form = NA)
    
    r2_val <- r.squaredGLMM(fit)[1, "R2m"]
    dat <- as.matrix((fit_summary$coefficients))
    get <- grepl("Total_HH", rownames(dat))
    
    if (nrow(fit_summary$coefficients) >= 2) {
      pval <- round(as.matrix(fit_summary$coefficients)[get, 5],3)
      pval_str <- ifelse(pval > 0.05, "ns", 
                         ifelse(pval < 0, "***",
                                ifelse(pval < 0.01, "**",
                                       ifelse(pval < 0.05, "*","m"))))
    } else {
      pval_str <- "NA"
    }
    
    if (nrow(fit_summary$coefficients) >= 2) {
      beta <- round(as.matrix(fit_summary$coefficients)[get, 1],3)
      beta_str <- paste("B =",beta)
    } else {
      beta <- NA
    }
    
    stat_label <- paste0(
      beta_str,
      pval_str,", ",
      "R² = ", format(r2_val, digits = 3)
    )
    
    n_ind <- length(unique(sub_df$individual))
    
    p <- ggplot(sub_df, aes(x = Total_HH , y = pLGT)) +
      geom_point(
        aes(fill = type, colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_line(
        data = pred_grid,
        aes(x = Total_HH , y = pred_pLGT),
        colour = "black",
        linewidth = 1
      ) +
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = stat_label, 
        hjust= 1.1,
        vjust = 1.5, 
        size = 4.5, 
        fontface = "italic"
      ) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
      ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)", x = "Total_HH  (count)")
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    return(p)
  })

plot_list <- compact(plot_list)

combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "LMM Model: pLGT ~ log10(Total_HH+0.1) + age + log10(n_contig) + (1|community); R2 is the marginal R2",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_Total_HH_lm", height = 6, width = 15, unit = "cm", scale=2)


### by children
plot_list <- df_relYoung %>%
  filter(age < 50) %>%
  split(.$site) %>%
  map(function(sub_df) {
    
    current_host  <- as.character(unique(sub_df$host))
    current_site  <- as.character(unique(sub_df$site))
    
    counts <- sub_df %>%
      group_by(with_young) %>%
      summarise(
        nSample=n_distinct(SAMPLE)
      )
    
    type_labels <- setNames(paste0(counts$with_young, "\n(", counts$nSample, ")"), counts$with_young)
    
    sub_df$x <- type_labels[sub_df$with_young]
    
    counts_young <- table(sub_df$with_young)
    valid_groups <- names(counts_young)[counts_young >= 2]
    
    n_ind <- length(unique(sub_df$individual))
    
    p <- ggplot(sub_df, aes(x = with_young, y = pLGT, fill = type)) +
      geom_point(
        aes(colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, fill = NA) +
      scale_x_discrete(labels = type_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
       ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)")

    if (length(valid_groups) >= 2) {
      p <- p + geom_signif(
        comparisons = list(c("with_young", "no_young")),
        test = "wilcox.test",
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        textsize = 5,
        vjust = -0.7,
        tip_length = 0.03
      )
    }
    
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  })

combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "Age < 50",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_with_young", height = 8, width = 15, unit = "cm", scale=2)

# Rel_young (continuous)
library(lmerTest)
library(MuMIn)

plot_list <- df_relYoung %>%
  #split(list(.$site, .$group)) %>%
  split(.$site) %>%
  map(function(sub_df) {
    # sub_df <- subset(df_relYoung, site=="Fecal")
    lower_bound <- quantile(sub_df$pLGT  , 0.005, na.rm = TRUE)
    upper_bound <- quantile(sub_df$pLGT  , 0.995, na.rm = TRUE)
    sub_df <- sub_df[sub_df$pLGT   >= lower_bound & sub_df$pLGT  <= upper_bound, ] %>%
      filter(age < 50)
    
    if (nrow(sub_df) == 0) return(NULL)
    
    current_host  <- as.character(unique(sub_df$host))[1]
    current_site  <- as.character(unique(sub_df$site))[1]
    
    fit <- lmer(pLGT ~ Rel_young + age + log10(n_contig) + (1|population), data = sub_df)
    fit_summary <- summary(fit)

    pred_grid <- data.frame(
      Rel_young  = seq(min(sub_df$Rel_young , na.rm = TRUE), max(sub_df$Rel_young , na.rm = TRUE), length.out = 100),
      age = mean(sub_df$age, na.rm=TRUE),
      n_contig = mean(sub_df$n_contig, na.rm = TRUE)
    )
    pred_grid$pred_pLGT <- predict(fit, newdata = pred_grid, re.form = NA)
    
    r2_val <- r.squaredGLMM(fit)[1, "R2m"]
    dat <- as.matrix((fit_summary$coefficients))
    get <- grepl("Rel_young", rownames(dat))
    
    if (nrow(fit_summary$coefficients) >= 2) {
      pval <- round(as.matrix(fit_summary$coefficients)[get, 5],3)
      pval_str <- ifelse(pval > 0.05, "ns", 
                         ifelse(pval < 0, "***",
                                ifelse(pval < 0.01, "**",
                                       ifelse(pval < 0.05, "*","m"))))
    } else {
      pval_str <- "NA"
    }
    
    if (nrow(fit_summary$coefficients) >= 2) {
      beta <- round(as.matrix(fit_summary$coefficients)[get, 1],3)
      beta_str <- paste("B =",beta)
    } else {
      beta <- NA
    }
    
    stat_label <- paste0(
      beta_str,
      pval_str,", ",
      "R² = ", format(r2_val, digits = 3)
    )
    
    n_ind <- length(unique(sub_df$individual))
    
    p <- ggplot(sub_df, aes(x = Rel_young , y = pLGT)) +
      geom_point(
        aes(fill = type, colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_line(
        data = pred_grid,
        aes(x = Rel_young , y = pred_pLGT),
        colour = "black",
        linewidth = 1
      ) +
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = stat_label, 
        hjust= 1.1,
        vjust = 1.5, 
        size = 4.5, 
        fontface = "italic"
      ) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
      ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)", x = "Rel_young  (count)")
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  })

plot_list <- compact(plot_list)

combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "LMM Model: pLGT ~ Rel_young + age + log10(n_contig) + (1|community); R2 is the marginal R2",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_Rel_young_lm", height = 6, width = 15, unit = "cm", scale=2)


#######################
# by Dog Ownership
######################
hh_dogs <- read.delim(file.path(input_dir, "Household_counts_filledv2.tsv")) %>%
  filter(!is.na(sampleid)) %>%
  select(individual = sampleid, Dogs_HH)

df_humans_wDogs <- df %>%
  filter(host == "human") %>%
  filter(group != "Agriculture") %>%
  inner_join(hh_dogs, by="individual") %>%
  mutate(with_dogs = ifelse(Dogs_HH>0, "with_dogs", "no_dogs")) %>%
  mutate(group2 = paste(group, with_dogs, sep = "_")) %>%
  mutate(group = droplevels(group),
         population = droplevels(population),
         with_dogs = factor(with_dogs, c("with_dogs", "no_dogs")),
         SAMPLE = as.character(SAMPLE)) %>%
  arrange(group, population, with_dogs, SAMPLE) %>%
  mutate(group2 = factor(group2, unique(group2))) %>%
  na.omit()

### by ownership only
plot_list <- df_humans_wDogs %>%
  split(.$site) %>%
  map(function(sub_df) {
    
    current_host  <- as.character(unique(sub_df$host))
    current_site  <- as.character(unique(sub_df$site))
    
    counts <- sub_df %>%
      group_by(with_dogs) %>%
      summarise(
        nSample=n_distinct(individual)
      )
    type_labels <- setNames(paste0(counts$with_dogs, "\n(", counts$nSample, ")"), counts$with_dogs)

    sub_df$x <- type_labels[sub_df$with_dogs]
    
    counts_dogs <- table(sub_df$with_dogs)
    valid_groups <- names(counts_dogs)[counts_dogs >= 2]
    
    n_ind <- length(unique(sub_df$individual))
    
    p <- ggplot(sub_df, aes(x = with_dogs, y = pLGT, fill = type)) +
      geom_point(
        aes(colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, fill = NA) +
      scale_x_discrete(labels = type_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
      ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)")
    
    if (length(valid_groups) >= 2) {
      p <- p + geom_signif(
        comparisons = list(c("with_dogs", "no_dogs")),
        test = "wilcox.test",
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
        textsize = 5,
        vjust = -0.7,
        tip_length = 0.03
      )
    }
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  })

combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_withDogs", height = 8, width = 15, unit = "cm", scale=2)

# nDogs
library(lmerTest)
library(MuMIn)

plot_list <- df_humans_wDogs %>%
  split(.$site) %>%
  #split(list(.$site, .$group)) %>%
  map(function(sub_df) {
    #sub_df <- subset(df_humans_wDogs, site=="Fecal")
    lower_bound <- quantile(sub_df$pLGT, 0.005, na.rm = TRUE)
    upper_bound <- quantile(sub_df$pLGT, 0.995, na.rm = TRUE)
    sub_df <- sub_df[sub_df$pLGT >= lower_bound & sub_df$pLGT <= upper_bound, ] %>%
      rename(nDogs = Dogs_HH)
    
    if (nrow(sub_df) == 0) return(NULL)
    
    current_host  <- as.character(unique(sub_df$host))[1]
    current_site  <- as.character(unique(sub_df$site))[1]
    
    fit <- lmer(pLGT ~ nDogs + log10(n_contig) + (1|population), data = sub_df)
    fit_summary <- summary(fit)
    
    pred_grid <- data.frame(
      nDogs = seq(min(sub_df$nDogs, na.rm = TRUE), max(sub_df$nDogs, na.rm = TRUE), length.out = 100),
      n_contig = mean(sub_df$n_contig, na.rm = TRUE)
    )
    pred_grid$pred_pLGT <- predict(fit, newdata = pred_grid, re.form = NA)
    
    r2_val <- r.squaredGLMM(fit)[1, "R2m"]
    dat <- as.matrix((fit_summary$coefficients))
    get <- grepl("nDogs", rownames(dat))
    
    if (nrow(fit_summary$coefficients) >= 2) {
      pval <- round(as.matrix(fit_summary$coefficients)[get, 5],3)
      pval_str <- ifelse(pval > 0.05, "ns", 
                         ifelse(pval < 0, "***",
                                ifelse(pval < 0.01, "**",
                                       ifelse(pval < 0.05, "*","m"))))
    } else {
      pval_str <- "NA"
    }
    
    if (nrow(fit_summary$coefficients) >= 2) {
      beta <- round(as.matrix(fit_summary$coefficients)[get, 1],3)
      beta_str <- paste("B =",beta)
    } else {
      beta <- NA
    }
    
    stat_label <- paste0(
      beta_str,
      pval_str,", ",
      "R² = ", format(r2_val, digits = 3)
    )
    
    n_ind <- length(unique(sub_df$individual))
    
    p <- ggplot(sub_df, aes(x = nDogs, y = pLGT)) +
      geom_point(
        aes(fill = type, colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_line(
        data = pred_grid,
        aes(x = nDogs, y = pred_pLGT),
        colour = "black",
        linewidth = 1
      ) +
      annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = stat_label, 
        hjust = 1.1, 
        vjust = 1.5, 
        size = 5, 
        fontface = "italic"
      ) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
      ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)", x = "nDogs (count)")
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  })

plot_list <- compact(plot_list)

combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = "LMM Model: pLGT ~ nDogs + log10(n_contig) + (1|community); R2 is the marginal R2)",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot
save_plot(combined_plot, "nLGT_by_nDogs_lm", height = 6, width = 15, unit = "cm", scale=2)

###################
#### by lifestyle
###################
plot_list <- df_humans_wDogs %>%
  split(list(.$site, .$with_dogs), drop = TRUE) %>%
  map(function(sub_df) {
    
    if (nrow(sub_df) == 0) return(NULL)
    
    current_host   <- as.character(unique(sub_df$host))[1]
    current_site   <- as.character(unique(sub_df$site))[1]
    current_subset <- as.character(unique(sub_df$with_dogs))[1]
    
    counts <- sub_df %>% 
      group_by(group2) %>%
      summarise(nSample = n_distinct(SAMPLE), .groups = "drop")
    
    type_labels <- setNames(
      paste0(counts$group2, "\n(", counts$nSample, ")"), 
      counts$group2
    )
    
    sub_df$x <- factor(type_labels[as.character(sub_df$group2)])
    
    unique_x <- unique(sub_df$x)
    
    if (length(unique_x) >= 2) {
      my_comparisons <- combn(as.character(unique_x), m = 2, simplify = FALSE)
    } else {
      my_comparisons <- NULL
    }
    
    n_ind <- length(unique(sub_df$individual))
    p <- ggplot(sub_df, aes(x = x, y = pLGT, fill = type)) +
      geom_point(
        aes(colour = type),
        shape = 21,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, seed = 1)
      ) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, fill = NA) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      scale_fill_manual(values = site_col) +
      scale_colour_manual(values = site_col) +
      theme_minimal(base_size = 18) +
      theme(
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14, face="bold"),
        plot.subtitle = element_text(size=12, face="italic")
      ) +
      labs(title = paste(current_host, current_site, sep=""), subtitle=paste("(n=", n_ind, ")", sep = ""), 
           y = "pLGT (%)")
    
    if (!is.null(my_comparisons)) {
      p <- p + stat_compare_means(
        comparisons = my_comparisons,
        method = "wilcox.test",
        label = "p.signif",
        size = 5
      )
    }
    
    if(current_site!="Oral"){
      p <- p + theme(axis.title.x = element_blank())
    }
    
    if(current_site != "Fecal"){
      p <- p + theme(axis.title.y = element_blank())
    }
    
    return(p)
  })

plot_list <- compact(plot_list)


ordered_plot_list <- plot_list[c(1, 4, 2, 5, 3, 6)]
ordered_plot_list <- plot_list
combined_plot <- wrap_plots(ordered_plot_list, ncol = 3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.margin = margin(0.5,0.5,0.5,2, "cm"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_plot

save_plot(combined_plot, "nLGT_by_Ls_withDogs", height = 15, width = 15, unit = "cm", scale=2)

