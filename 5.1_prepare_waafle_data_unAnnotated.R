rm(list=ls())

basedir="Documents/Analysis/MOBILE/WAAFLE_Extra"
setwd(basedir)

type=c("humanFecal","dogFecal","humanOral","dogOral","humanSkin","dogSkin")
cType=c("DHF","DDF","DHV","DDV","DHP","DDH")
output_dir=file.path(basedir,"index")

input_dir="Documents/Analysis/MOBILE/MOBILE_microbiome/input_files/waafle"
waffle_file="raw_all_internal_merged_filtered.lgt.tsv"

dir.create(output_dir,showWarnings = F,recursive = T)

# load WAAFLE output

df <- do.call(rbind, lapply(type, function(x) read.delim(file.path(input_dir,x,waffle_file))))

rownames(df) <- NULL
df$original_sid <- df$SAMPLE
df$SAMPLE <- substr(df$original_sid, 1,6)
df$type <- substr(df$SAMPLE, 1,3)

table(substr(unique(substr(df$SAMPLE,1,6)),1,3))

##############
## Indexing ##
##############
# Index unique synteny and output the gene positional table
synteny <- data.frame(SYNTENY=unique(df$SYNTENY), LEN=nchar(unique(df$SYNTENY)))
synteny <- synteny[order(synteny$SYNTENY,decreasing = T),]
synteny <- synteny[order(synteny$LEN,decreasing = F),]
synteny <- unique(synteny$SYNTENY)


# Main table: Position details
results <- list()
for (i in seq_along(synteny)) {
  chars <- strsplit(synteny[i], "")[[1]]
  for (j in seq_along(chars)) {
    if (chars[j] %in% c("A", "B", "*", "~")) {
      results[[length(results) + 1]] <- data.frame(
        Synteny_Index = i,
        Position = j,
        Character = chars[j]
      )
    }
  }
}
synteny_position_table <- do.call(rbind, results)

# Separate table: String index definitions
synteny_index_table <- data.frame(
  Synteny_Index = seq_along(synteny),
  SYNTENY = synteny
)
directions <- data.frame(SYNTENY=df$SYNTENY, DIRECTION=df$DIRECTION)
directions$SYNTENY <- factor(directions$SYNTENY, synteny)
directions <- directions[order(directions$SYNTENY),]
directions <- unique(directions)
synteny_index_table <- merge.data.frame(synteny_index_table, directions, by="SYNTENY")
synteny_index_table <- synteny_index_table[order(synteny_index_table$Synteny_Index),]
synteny_index_table <- unique(synteny_index_table)
rownames(synteny_index_table) <- NULL


# Index and get the positions of genes
# Input: Vector of gene strings
gene_syn <- data.frame(SYNTENY=df$SYNTENY, Gene_SYNTENY=df$ANNOTATIONS.UNIREF90, LOCI_STRING=df$LOCI)
gene_syn <- unique(gene_syn)
gene_syn$SYNTENY <- factor(gene_syn$SYNTENY,levels = synteny_index_table$SYNTENY)
gene_syn <- gene_syn[order(gene_syn$Gene_SYNTENY),]
gene_syn <- gene_syn[order(gene_syn$SYNTENY),]
rownames(gene_syn) <- NULL

gene_strings <- unique(gene_syn$Gene_SYNTENY)

position_list <- list()
for (i in seq_along(gene_strings)) {
  genes <- strsplit(gene_strings[i], "\\|")[[1]]
  loci <- strsplit(gene_syn$LOCI_STRING[gene_syn$Gene_SYNTENY == gene_strings[i]][1], "\\|")[[1]]
  clade <- strsplit(as.character(gene_syn$SYNTENY[gene_syn$Gene_SYNTENY == gene_strings[i]][1]), "")[[1]]
  
  for (j in seq_along(genes)) {
    locus_parts <- strsplit(loci[j], ":")[[1]]
    position_list[[length(position_list) + 1]] <- data.frame(
      Gene_Synteny_Index = i,
      Position = j,
      Gene_Name = genes[j],
      Clade = clade[j],
      Start = as.numeric(locus_parts[1]),
      Stop = as.numeric(locus_parts[2])
    )
  }
}
gene_position_table <- do.call(rbind, position_list)

gene_index_table <- unique(gene_syn)
gene_index_table$Gene_Synteny_Index <- match(gene_index_table$Gene_SYNTENY, gene_strings)

gene_index_table <- merge.data.frame(gene_index_table, synteny_index_table, by="SYNTENY")
gene_index_table <- gene_index_table[order(gene_index_table$Gene_Synteny_Index),]
gene_index_table <- gene_index_table[order(gene_index_table$Synteny_Index),]
gene_index_table <- unique(gene_index_table)
rownames(gene_index_table) <- NULL

# Output
head(synteny_position_table)
head(synteny_index_table)
head(gene_position_table)
head(gene_index_table)

# Write DB
write.table(file=file.path(output_dir,"synteny_position_table_unAnnot.tsv"), synteny_position_table,
            quote = F, sep = "\t",row.names =F)
write.table(file=file.path(output_dir,"synteny_index_table_unAnnot.tsv"), synteny_index_table,
            quote = F, sep = "\t",row.names =F)
write.table(file=file.path(output_dir,"gene_position_table_unAnnot.tsv"), gene_position_table,
            quote = F, sep = "\t",row.names =F)
write.table(file=file.path(output_dir,"gene_index_table_unAnnot.tsv"), gene_index_table,
            quote = F, sep = "\t",row.names =F)
