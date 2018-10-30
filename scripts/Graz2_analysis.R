library(oligo)
library(hugene20sthsentrezg.db)
library(arrayQualityMetrics)
library(sva)
library(stringr) 
library(ggplot2)
library(ggfortify)
library(limma)
source("plots_utils.R")


##install.brainarray("hugene20st", version = "22.0.0")
#install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.hugene20st.hs.entrezg_22.0.0.tar.gz",
#                 repos = NULL, type = "source")

pd <- read.table("../pdata/pdata_FHR_Graz2.tsv", header = TRUE, sep = "\t", stringsAsFactors = F)
pd <- pd[which(pd$GroupB!="SABGEI"),]
pd <- pd[which(pd$GroupB!="KRIK"),]

celFiles <- paste("../raws/20160711-FHR-GRAZ2/", rownames(pd), sep="")

rawData <- oligo::read.celfiles(filenames=celFiles, pkgname="pd.hugene20st.hs.entrezg")

row.names(pd) <- sampleNames(rawData)
pdata <- AnnotatedDataFrame(data=pd)
phenoData(rawData) <- pdata

# Add ScanDate to pdata
#pd.plain$ScanDate <- str_replace_all(arawData@protocolData@data$dates, "T", " ")
#pd.plain$ScanDate <- sapply(strsplit(pd.plain$ScanDate, split=' ', fixed=TRUE), function(x) (x[1]))

eset = oligo::rma(rawData)

arrayQualityMetrics(expressionset = eset,
                    outdir = "../plots/graz2/AQM_report_graz2_nooutlier",
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("GroupA"))


exprs <- exprs(eset)
colnames(exprs) <- pd$SampleName
anno <- AnnotationDbi::select(hugene20sthsentrezg.db, rownames(exprs), c("ENTREZID", "SYMBOL", "GENENAME"))
anno <- anno[which(anno$ENTREZID!="NA"),]
n_occur <- data.frame(table(anno$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
anno <- anno[which(anno$PROBEID %in% uniques),]

exprs <- exprs[which(rownames(exprs) %in% anno$PROBEID),]
rownames(exprs) <- anno[match(rownames(exprs), anno$PROBEID),]$SYMBOL


## Simple PCA

pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, c("GroupA", "GroupB"))
save_plot(paste("../plots/graz2/Graz2_PCA_nooutliers.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)
save_plot(paste("../plots/graz2/Graz2_PCA_nooutliers.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)

## Fancy PCA

library(factoextra)
fviz_eig(pca)

pl <- fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
pl
save_plot(paste("../plots/graz2/Graz2_PCA_cos_nooutliers.pdf", sep=""),
          base_height=5, pl, ncol=1)
save_plot(paste("../plots/graz2/Graz2_PCA_cos_nooutliers.svg", sep=""),
          base_height=5, pl, ncol=1)


indCol <- brewer.pal(length(levels(factor(pd$GroupB))),"Set3")
individuals <- as.character(factor(pd$GroupB, labels=indCol))

pl <- fviz_pca_biplot(pca,
                geom=c("point", "text"),
                label="all",
                col.var = "contrib",
                fill.ind = factor(pd$GroupB),
                #palette = indCol,
                col.ind = "black",
                pointshape=21, pointsize = 3,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,     # Avoid text overlapping
                select.var = list(contrib=5),
                invisible = "quali"
) + labs(fill = "Donors")
pl
save_plot(paste("../plots/graz2/Graz2_PCA_loadings_nooutliers.pdf", sep=""),
          base_height=6, pl, ncol=1)
save_plot(paste("../plots/graz2/Graz2_PCA_loadings_nooutliers.svg", sep=""),
          base_height=6, pl, ncol=1)

# DEGs
filterDEGS <- function(degs, pval, fc, adj) {
  # Filter by p-values
  if (missing(adj)) {
    degs <- degs[degs$adj.P.Val < pval,]
  } else if (adj==FALSE) {
    degs <- degs[degs$P.Value < pval,]
  } else if (adj==TRUE) {
    degs <- degs[degs$adj.P.Val < pval,]
  }
  
  # Sort by logFC
  degs <- degs[order(abs(degs$logFC), decreasing = TRUE),]
  # Filter by logFC
  degs <- degs[abs(degs$logFC) > fc,]  
  return (degs)
}


design = model.matrix(~GroupA+GroupB, data=pd)
#design = model.matrix(~GroupA, data=pd)
design
fit <- lmFit(exprs, design)
fit <- eBayes(fit)
# Get all the genes with logFC, p-values, no filtering
degs <- topTable(fit, adjust.method="fdr", number=nrow(fit), coef="GroupAp")

degs <- filterDEGS(degs, 0.01, 0.5, adj=F)

anno_name<-select(org.Hs.eg.db, rownames(degs), c("GENENAME"), keytype="SYMBOL")
degs$description <- anno_name$GENENAME
degs$SYMBOL <- rownames(degs)

degs$FC <- sign(degs$logFC)*(2^abs(degs$logFC))

degs <- degs[,c("SYMBOL", "logFC", "FC", "P.Value", "adj.P.Val", "description")]
  
write.table(degs, "../degs/degs_graz2.tsv",sep="\t", quote=F, row.names=FALSE)

