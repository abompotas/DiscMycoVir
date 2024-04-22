library(ggplot2)
library(gplots)
library(dplyr)
library(magrittr)
library(stringr)
library(topGO)
library(KEGGREST)
library(biomaRt)
library(org.Hs.eg.db)
library(Rgraphviz)
library(GOxploreR)
library(GO.db)
library(igraph)

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).", call = FALSE)
  quit(save = 'no')
}
input_file <- args[1]
type <- args[2]
uniprotids <- args[3]
output_path <- args[4]

data1 <- read.csv(input_file, header = TRUE, sep = '\t')
data1$GeneRatio <- log(data1$GeneRatio)
datanew <- data1 %>% dplyr::filter(p.adjust < 0.1)

# Extract the relevant columns
go_terms <- datanew$Term
go_terms <- gsub(".*~", "", go_terms)

datanew$Term <- go_terms

plot_height <- max(10, min(2 * nrow(datanew), 25))

jpeg(filename = paste(output_path, paste(type, "enrichment_plot.jpg", sep = '_'), sep = ''), units = "cm", width = 40, height = plot_height, res = 600)
par(cex = 2.0)
ggplot(datanew,
       aes(x = reorder(Term, GeneRatio), y = GeneRatio)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 18) +
  scale_colour_gradient(limits = c(0, 0.1), low = "red", high = "blue") +
  coord_flip() +
  labs(y = 'GeneRatio', x = 'Term') +
  xlab(NULL) +
  theme(
    axis.text.x = element_text(size = 19, vjust = 0.5),
    axis.text.y = element_text(size = 19)
  )

dev.off()


# Read the UniProt IDs from the file
uniprot_ids <- read.table(uniprotids, header = FALSE, stringsAsFactors = FALSE)
uniprot_ids <- uniprot_ids$V1
# Map UniProt IDs to Entrez Gene IDs using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapped_ids <- getBM(attributes = c("uniprotswissprot", "entrezgene_id"), filters = "uniprotswissprot", values = uniprot_ids, mart = mart)

# Create a list of Entrez Gene IDs for the genes of interest
gene_list <- as.character(mapped_ids$entrezgene)

# Create a list of all human genes
all_genes <- unique(mapped_ids$entrezgene_id)

# Create a named numeric vector indicating whether a gene is in the list of interest
gene_in_list <- ifelse(all_genes %in% gene_list, 1, 0)
names(gene_in_list) <- all_genes

# Convert the named numeric vector to a factor with 2 levels
gene_in_list <- factor(gene_in_list, levels = c(0, 1), labels = c("uninteresting", "interesting"))

if (type == 'GOTERM_BP_DIRECT') {
  ontology = "BP"
}else if (type == 'GOTERM_CC_DIRECT') {
  ontology = 'CC'
}else if (type == 'GOTERM_MF_DIRECT') {
  ontology = 'MF'
}else {
  ontology = FALSE
}

if (ontology != FALSE) {
  # Create the topGO data object
  GOdata <- new("topGOdata", ontology = ontology, allGenes = gene_in_list, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'entrez')
  test <- new("elimScore", testStatistic = GOKSTest, name = 'Fisher test', cutOff = 0.01)
  resultfisher <- getSigGroups(GOdata, test)
  tab <- GenTable(GOdata, raw.p.value = resultfisher, topNodes = length(resultfisher@score), numChar = 120)

  jpeg(filename = paste(output_path, paste(ontology, "gotree_plot.jpg", sep = '_'), sep = ''), units = 'cm', width = 10, height = 10, res = 1000)
  par(cex = 0.8)
  plot <- showSigOfNodes(GOdata, score(resultfisher), firstSigNodes = 3, useInfo = "def")
  print(plot)
  dev.off()

  if (ontology == "BP") {
    terms_level3 <- Level2GOTermBP(level = 2, organism = 'Homo sapiens')
  } else if (ontology == "MF") {
    terms_level3 <- Level2GOTermMF(level = 2, organism = 'Homo sapiens')
  } else {
    terms_level3 <- Level2GOTermCC(level = 2, organism = 'Homo sapiens')
  }


  common_mf <- intersect(tab$GO.ID, terms_level3)
  GO_table_subset <- tab[tab$GO.ID %in% intersect(tab$GO.ID, terms_level3),]

  # Sort the rows of the table by the Annotated column in descending order
  GO_table_sorted <- GO_table_subset[order(-GO_table_subset$Annotated),]


  # Keep only the top 10 rows
  GO_table_top15 <- head(GO_table_sorted, 15)

  jpeg(filename = paste(output_path, paste(ontology, "barplot_third_level_plot.jpg", sep = '_'), sep = ''), units = 'cm', width = 40, height = 20, res = 600)
  par(cex = 2.0)


  # Create the bar plot
  plot <- ggplot(data = GO_table_top15, aes(x = Annotated, y = Term)) +
    geom_bar(stat = "identity", fill = 'steelblue') +
    labs(x = "Number of genes", y = "GO Term") +
    geom_text(aes(label = Annotated), hjust = -0.2, size = 5) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x = element_text(size = 17, vjust = 0.5),
      axis.text.y = element_text(size = 17)
    )
  print(plot)
  dev.off()

}




