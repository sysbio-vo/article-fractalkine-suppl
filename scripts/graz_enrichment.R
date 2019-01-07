library(enrichR)
library(dplyr)
library(biomaRt)
library(karyoploteR) 
library(regioneR)
source("plots_utils.R")

# Read data
degs <- read.table("../degs/degs_graz.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")
exprs <- read.table("../exprs/exprs_graz.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")
enrich <- read.table("../degs/graz_enrichment_matrix_shinygo_enrichr.tsv",sep="\t", header = T,
                             stringsAsFactors = F, quote = "", check.names = F)
proteins <- read.table("../degs/degs_graz_known_proteins.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")
genes <- read.table("../degs/degs_graz_known_genes.tsv",sep="\t", header = T,
                    stringsAsFactors = F, quote = "")
chloc <- read.table("../degs/degs_graz_chrom_loc.tsv",sep="\t", header = T,
                    stringsAsFactors = F, quote = "")


# EnrichR
dbs <- listEnrichrDbs()
dbs <- c("MGI_Mammalian_Phenotype_2017")
enriched <- enrichr(genes$SYMBOL, dbs)
#printEnrich(enriched, "output.txt" , sep = "\t", columns = c(1:9))

bp <- enriched[["MGI_Mammalian_Phenotype_2017"]]
mpo <- bp[which(bp$Combined.Score>10),]

mpo$Category <- "MPO"
l <- t(sapply(mpo$Term, function(x) {
                    v <- unlist(strsplit(x, "_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
                    id <- v[1]
                    name <- paste(v[-1], collapse=" ")
                    return(c(id, name))
                 }))

mpo$ID <- l[,1]
mpo$Term <- l[,2]

gen <- unlist(strsplit(paste(mpo$Genes, collapse=";"), ";"))
gen <- unique(gen)
m <- matrix(0, ncol = nrow(mpo), nrow = length(gen))
df <- data.frame(m, row.names = gen)
colnames(df) <- mpo$Term

g <- strsplit(mpo$Genes, ";")
for (i in 1:length(g)) {
  df[which(rownames(df) %in% g[[i]]), i] <- 1
}

df$FC <- degs$FC[match(rownames(df), degs$SYMBOL)]  

#GOHeat(df[,-8], nlfc = 0)
pl <- GOHeat(df, nlfc = 1, fill.col = c('red', 'gray95', 'blue'))
save_plot(paste("../plots/mpo_heatmap_enrichr.pdf", sep=""),
          base_height=5, base_width=10, pl, ncol=1)
save_plot(paste("../plots/mpo_heatmap_enrichr.svg", sep=""),
          base_height=5, base_width=10, pl, ncol=1)

# ShinyGO & enrichR
enrich$FC <- degs$FC[match(rownames(enrich), degs$SYMBOL)]
enrich[is.na(enrich)] <- 0
pl <- GOHeat(enrich, nlfc = 1, fill.col = c('red', 'gray95', 'blue'))
save_plot(paste("../plots/enrichment_heatmap.pdf", sep=""),
          base_height=5, base_width=9, pl, ncol=1)
save_plot(paste("../plots/enrichment_heatmap.svg", sep=""),
          base_height=5, base_width=9, pl, ncol=1)

# Chromosomal locations
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
## Filter on HGNC symbol, retrieve genomic location and band
my.symbols <- degs$SYMBOL
my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol"),
                    values = list(hgnc_symbol=my.symbols),
                    mart = ensembl)

joint <- merge(degs, my.regions, by.x="SYMBOL", by.y="hgnc_symbol", all=T)
joint <- joint[order(abs(joint$logFC), decreasing = T),]
## Not all the locations were present, I had to edit the file manually afterwards.
#write.table(joint, "../degs/degs_graz_chrom_loc.tsv",sep="\t", quote=F, row.names=F, col.names = T)

## Download and unpack this file - ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gene.info <- read.table("../misc/Homo_sapiens.gene_info",sep="\t", header = T,
                    stringsAsFactors = F, quote = "", fill=F, check.names = F, allowEscapes = T,
                    comment.char = "")

gene.info <- gene.info[,c("GeneID", "type_of_gene")]
degs.full <- merge(chloc, gene.info, by.x="ENTREZID", by.y="GeneID", all.x=T)
degs.full[which(degs.full$SYMBOL=="LOC105369845"),]$type_of_gene <- "ncRNA"
write.table(degs.full, "../degs/degs_graz_full_info.tsv",sep="\t", quote=F, row.names=F, col.names = T)

ch.size <- read.table("../misc/hg38_ch_sizes.tsv",sep="\t", header = T,
                        stringsAsFactors = F, quote = "", fill=F, check.names = F, allowEscapes = T,
                        comment.char = "")

## Plot with karyotypeR package
data <- toGRanges(degs.full[,c("chromosome_name", "start_position", "end_position", "SYMBOL", "type_of_gene")])
seqlevelsStyle(data) <- "UCSC"

typeCol <- brewer.pal(length(levels(factor(data$type_of_gene))),"Set2")
type <- factor(data$type_of_gene, labels=typeCol)

ncod <- data[which(data$type_of_gene=="ncRNA"),]
cod <- data[which(data$type_of_gene!="ncRNA"),]
seqlevelsStyle(ncod) <- "UCSC"
seqlevelsStyle(cod) <- "UCSC"

pdf("../plots/ch_location.pdf", width=17, height=13)
kp <- plotKaryotype(genome="hg38", plot.type = 2)
kpPlotMarkers(kp, data=cod, labels=cod$SYMBOL,
              text.orientation = "horizontal", label.color = "blue",
              marker.parts = c(0.2, 0.7, 0.1), r1=0.7, cex=0.7,
              adjust.label.position = T, clipping = F, data.panel=1)
kpPlotMarkers(kp, data=ncod, labels=ncod$SYMBOL,
              text.orientation = "horizontal", label.color = "red",
              marker.parts = c(0.2, 0.7, 0.1), r1=0.7, cex=0.7,
              adjust.label.position = T, clipping = F, data.panel=2)
dev.off()
