library(enrichR)
library(dplyr)
library(biomaRt)
library(karyoploteR) 
library(regioneR)
library(EnhancedVolcano)
library(TissueEnrich)
library(tidyr)
library(pheatmap)
library(org.Hs.eg.db)
library(cowplot)
source("plots_utils.R")

# Read data
degs <- read.table("../degs/degs_graz.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")
degs.all <- read.table("../degs/degs_graz_allgenes.tsv",sep="\t", header = T,
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





## Tissue enrichment
gs <- GeneSet(geneIds=degs$SYMBOL,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
bs <- GeneSet(geneIds=degs.all$SYMBOL,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
diff <- setdiff(geneIds(gs), intersect(geneIds(gs),geneIds(bs)))

gs <- GeneSet(geneIds=setdiff(degs$SYMBOL, diff),organism="Homo Sapiens",geneIdType=SymbolIdentifier())

output<-teEnrichment(inputGenes = gs, backgroundGenes = bs)

seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
head(enrichmentOutput)
enrichmentOutput <- enrichmentOutput[which(enrichmentOutput$Tissue.Specific.Genes>0),]
enrichmentOutput <- enrichmentOutput[which(enrichmentOutput$Log10PValue>0),]

### pval plot
ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

### FC plot
ggplot(enrichmentOutput,aes(x=reorder(Tissue,-fold.change),y=fold.change,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = 'Fold change')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

### Heatmap for tissue specific genes
seExp<-output[[2]][["Bone Marrow"]]
exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
exp$Gene<-row.names(exp)
exp<-exp %>% gather(Tissue=1:(ncol(exp)-1))

ggplot(exp, aes(key, Gene)) + geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")+
  labs(x='', y = '') +
  theme_bw() +
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())

### Genes for specific tissue
seGroupInf<-output[[3]][["Placenta"]]
groupInf<-data.frame(assay(seGroupInf))
print(head(groupInf))

### Heatmap for all genes acrosstissue specific genes
exps <- c()
hgenes <- c()
grInf <- data.frame(assay(output[[3]][[i]]))

for (i in 1:length(output[[2]])) {
  seExp<-output[[2]][[i]]
  exp<-setNames(data.frame(assay(seExp), row.names = rowData(seExp)[,1]), colData(seExp)[,1])
  hgenes <- c(hgenes, rownames(exp))
  #hgenes <- unique(hgenes)
  #print(hgenes)
  seGroupInf<-output[[3]][[i]]
  groupInf<-data.frame(assay(seGroupInf))
  grInf <- rbind(grInf, groupInf)
}


### Unmapped genes
unmapped <- print(geneIds(output[[4]]))

all <- c(as.character(unique(grInf$Gene)), unmapped)
degs[-which(degs$SYMBOL %in% all),]$SYMBOL

### Heatmap for all the genes
load("../misc/combine-expression.rda")
hpa <- dataset[[1]]$expressionData
anno<-select(org.Hs.eg.db, rownames(hpa), c("GENENAME", "SYMBOL"), keytype="ENSEMBL")
anno <- anno[which(anno$SYMBOL %in% degs$SYMBOL),]
hpa <- hpa[which(rownames(hpa) %in% anno$ENSEMBL), ]

ind <- match(anno$ENSEMBL, rownames(hpa))
rownames(hpa) <- anno$SYMBOL
hpa <- hpa[apply(hpa, MARGIN = 1, FUN = function(x) sd(x) != 0),]
exprs.hpa <- exprs[which(exprs$SYMBOL %in% rownames(hpa)),]
rownames(exprs.hpa) <- exprs.hpa$SYMBOL
exprs.hpa <- exprs.hpa[,-c(1,2)]
ind <- match(rownames(hpa), rownames(exprs.hpa))
exprs.hpa <- exprs.hpa[ind,]
exprs.m <- rowSums(exprs.hpa[,c(1, 3, 5, 7)])
exprs.p <- rowSums(exprs.hpa[,c(2, 4, 6, 8)])

hpa$CD16 <-exprs.m
hpa$CD16.stimulated <-exprs.p

hpa <- hpa[,c("Bone.Marrow", "Heart.Muscle", "Kidney",
              "Placenta","Smooth.Muscle","Spleen",
              "CD16", "CD16.stimulated")]

pl <- pheatmap(hpa, cluster_cols = T, cluster_rows = T, drop_levels = F,
               scale = "row", treeheight_row = 40, treeheight_col = 20,
               cellheight = 4, cellwidth=17, border_color = NA,
               fontsize_row = 4, fontsize = 8, silent=T)

save_plot(paste("../plots/tissue/degs_graz_heatmap_tissue.pdf", sep=""),
          base_height=5, base_width=5, pl, ncol=1)
save_plot(paste("../plots/tissue/degs_graz_heatmap_tissue.svg", sep=""),
          base_height=5, base_width=5, pl, ncol=1)



### Volcano plot
pl <- EnhancedVolcano(degs.all,
                      lab = degs.all$SYMBOL,
                      x = "logFC",
                      y = "P.Value",
                      FCcutoff = 0.5,
                      pCutoff = 0.01,
                      ylim=c(0, 6),
                      selectLab = rownames(enrich)
                      )

save_plot(paste("../plots/volcano.pdf", sep=""),
          base_height=9, base_width=9, pl, ncol=1)
save_plot(paste("../plots/volcano.svg", sep=""),
          base_height=9, base_width=9, pl, ncol=1)

exprs_degs <- exprs[which(exprs$ENTREZID %in% degs$ENTREZID),]
write.table(exprs_degs, "../exprs/exprs_degs_graz.tsv",sep="\t", quote=F, row.names=F)
