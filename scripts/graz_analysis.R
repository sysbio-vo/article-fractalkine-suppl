library(oligo)
library(hugene20sthsentrezg.db)
library(limma)
library(arrayQualityMetrics)
library(ggplot2)
library(ggfortify)
library(factoextra)
library(pheatmap)
source("plots_utils.R")
source("degs_utils.R")

# Install custom CDF package if not installed
#install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/pd.hugene20st.hs.entrezg_22.0.0.tar.gz",
#                 repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hugene20sthsentrezg.db_22.0.0.tar.gz",
#                 repos = NULL, type = "source")

# Read pheno and raw data
pd <- read.table("../pdata/pdata_FHR_graz.tsv", header = TRUE, sep = "\t", stringsAsFactors = F)

celFiles <- paste("../raws/", pd$FileName, sep="")
rawData <- oligo::read.celfiles(filenames=celFiles, pkgname="pd.hugene20st.hs.entrezg")

rownames(pd) <- sampleNames(rawData)
pdata <- AnnotatedDataFrame(data=pd)
phenoData(rawData) <- pdata

# Normalize using RMA function
eset = oligo::rma(rawData)

# Create quality control report
arrayQualityMetrics(expressionset = eset,
                    outdir = paste("../plots/AQM_report_graz", sep=""),
                    force = TRUE,
                    do.logtransform = FALSE,
                    intgroup = c("Donor"))

# Extract expression matrix
exprs <- exprs(eset)
colnames(exprs) <- pd$SampleID
rownames(pd) <- pd$SampleID

# Get annotations
anno <- AnnotationDbi::select(hugene20sthsentrezg.db, rownames(exprs), c("ENTREZID", "SYMBOL", "GENENAME"))
anno <- anno[which(anno$ENTREZID!="NA"),]
n_occur <- data.frame(table(anno$PROBEID))
uniques <- n_occur[n_occur$Freq == 1,]$Var1
anno <- anno[which(anno$PROBEID %in% uniques),]

# Replace expression matrix row names with gene symbol
exprs <- exprs[which(rownames(exprs) %in% anno$PROBEID),]
ind <- match(rownames(exprs), anno$PROBEID)
rownames(exprs) <- anno[ind,]$SYMBOL

e <- cbind(ENTREZID=anno[ind,]$ENTREZID, SYMBOL=anno[ind,]$SYMBOL, exprs)
write.table(e, "../exprs/exprs_graz.tsv",sep="\t", quote=F, row.names=F)

# Simple PCA
pca = prcomp(t(exprs))
pl <- pcaPlots(pca, pd, c("Donor", "Stimulus"))
save_plot(paste("../plots/graz_PCA.pdf", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)
save_plot(paste("../plots/graz_PCA.svg", sep=""),
          base_height=3, base_aspect_ratio = pl[[2]]/2, pl[[1]], ncol=2)

# PCA percentage of variance
pl <- fviz_eig(pca)
save_plot(paste("../plots/graz_PCA_var_percentage.pdf", sep=""),
          base_height=3, base_width=6, pl, ncol=1)
save_plot(paste("../plots/graz_PCA_var_percentage.svg", sep=""),
          base_height=3, base_width=6, pl, ncol=1)

# PCA loadings
pl <- ggloadings(pca, num_pca=1)
save_plot(paste("../plots/graz_PCA_loadings_1.pdf", sep=""),
          base_height=3, base_width=6, pl, ncol=1)
save_plot(paste("../plots/graz_PCA_loadings_1.svg", sep=""),
          base_height=3, base_width=6, pl, ncol=1)
pl <- ggloadings(pca, num_pca=2)
save_plot(paste("../plots/graz_PCA_loadings_2.pdf", sep=""),
          base_height=3, base_width=6, pl, ncol=1)
save_plot(paste("../plots/graz_PCA_loadings_2.svg", sep=""),
          base_height=3, base_width=6, pl, ncol=1)

# PCA quality of representation
pl <- fviz_pca_ind(pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE)
save_plot(paste("../plots/graz_PCA_cos.pdf", sep=""),
          base_height=5, pl, ncol=1)
save_plot(paste("../plots/graz_PCA_cos.svg", sep=""),
          base_height=5, pl, ncol=1)

# PCA loadings on 2D
indCol <- brewer.pal(length(levels(factor(pd$Donor))),"Set3")
individuals <- as.character(factor(pd$Donor, labels=indCol))

pl <- fviz_pca_biplot(pca,
                geom=c("point", "text"),
                label="all",
                col.var = "contrib",
                fill.ind = factor(pd$Donor),
                #palette = indCol,
                col.ind = "black",
                pointshape=21, pointsize = 3,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,     # Avoid text overlapping
                select.var = list(contrib=5),
                invisible = "quali"
) + labs(fill = "Donor")
save_plot(paste("../plots/graz_PCA_loadings.pdf", sep=""),
          base_height=6, pl, ncol=1)
save_plot(paste("../plots/graz_PCA_loadings.svg", sep=""),
          base_height=6, pl, ncol=1)

# Get differentially expressed genes
pd$Stimulus <- factor(pd$Stimulus, levels=c("none", "CX3CL1"))
design = model.matrix(~Stimulus+Donor, data=pd)
rownames(exprs) <- anno[ind,]$ENTREZID
fit <- lmFit(exprs, design)
fit <- eBayes(fit)

# Get all the genes with logFC, p-values, no filtering
degs <- topTable(fit, adjust.method="fdr", number=nrow(fit), coef="StimulusCX3CL1")
# Filter
degs <- filterDEGS(degs, 0.01, 0.5, adj=F)

# Add description
anno_name<-select(org.Hs.eg.db, rownames(degs), c("GENENAME", "SYMBOL"), keytype="ENTREZID")
degs$description <- anno_name$GENENAME
degs$SYMBOL <- anno_name$SYMBOL
degs$ENTREZID <- rownames(degs)

# Add FC and rearrange columns
degs$FC <- sign(degs$logFC)*(2^abs(degs$logFC))
degs <- degs[,c("SYMBOL", "ENTREZID","logFC", "FC", "P.Value", "adj.P.Val", "description")]
  
write.table(degs, "../degs/degs_graz.tsv",sep="\t", quote=F, row.names=FALSE)

# Heatmap on significant genes
rownames(exprs) <- anno[ind,]$SYMBOL
exprs.degs <- exprs[which(rownames(exprs) %in% degs$SYMBOL),]
ht_matrix <- exprs.degs

anno_col <- pd[,c("Donor", "Stimulus")]

pl <- pheatmap(ht_matrix, cluster_cols = T, cluster_rows = T, drop_levels = F,
         annotation_col = anno_col, scale = "row", treeheight_row = 40, treeheight_col = 20,
         cellheight = 4, cellwidth=17, border_color = NA,
         fontsize_row = 4, fontsize = 8, silent=T)
save_plot(paste("../plots/degs_graz_heatmap.pdf", sep=""),
          base_height=9.5, base_width=5, pl, ncol=1)
save_plot(paste("../plots/degs_graz_heatmap.svg", sep=""),
          base_height=9.5, base_width=5, pl, ncol=1)