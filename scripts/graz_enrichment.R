library(enrichR)
library(GOplot)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)


# Read data

#degs <- read.table("../degs/degs_graz_all.tsv",sep="\t", header = T,
#                   stringsAsFactors = F, quote = "")
degs <- read.table("../degs/degs_graz.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")
exprs <- read.table("../exprs/exprs_graz.tsv",sep="\t", header = T,
                   stringsAsFactors = F, quote = "")


### EnrichR
dbs <- listEnrichrDbs()
dbs <- c("MGI_Mammalian_Phenotype_2017")
enriched <- enrichr(degs$SYMBOL, dbs)
#printEnrich(enriched, "output.txt" , sep = "\t", columns = c(1:9))

bp <- enriched[["MGI_Mammalian_Phenotype_2017"]]
head(bp)

bp <- bp[which(bp$P.value<0.01),]

mpo <- bp
mpo$Category <- "MPO"

l <- t(sapply(mpo$Term, function(x) {
                    v <- unlist(strsplit(x, "_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
                    id <- v[1]
                    name <- paste(v[-1], collapse=" ")
                    return(c(id, name))
                 }))

mpo$ID <- l[,1]
mpo$Term <- l[,2]
mpo$Genes <- gsub(";", ", ", mpo$Genes)
mpo$adj_pval <- mpo$P.value
mpo <- mpo[,c("Category", "ID", "Term", "Genes", "adj_pval")]

### topGO

genesList <- degs$P.Value
names(genesList) <- degs$ENTREZID
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = genesList, geneSel = function(x){x<=0.05},
                    nodeSize = 5,
                    annot = annFUN.org,mapping="org.Hs.eg.db", ID = "entrez")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultFisher
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
resultKS.elim
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 15)

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 10, useInfo = 'all')
numSigGenes(sampleGOdata)

### GOplot

data(EC)
EC$genes
head(EC$david)

go_obj <- EC
go_obj$eset <- exprs
go_obj$enrichr <- mpo
go_obj$process <- mpo$Term
go_obj$genes <- degs[,c("SYMBOL", "logFC")]
colnames(go_obj$genes) <- c("ID", "logFC")
go_obj$geneslist <- degs[,c("SYMBOL", "logFC")]

circ <- circle_dat(go_obj$enrichr, EC$genes)
chord <- chord_dat(circ, go_obj$genes, go_obj$process)
GOHeat(chord[,-8], nlfc = 0)

GOHeat(chord, nlfc = 1, fill.col = c('red', 'gray95', 'green'))

