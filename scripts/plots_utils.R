require(ggplot2)
require(cowplot)
require(ggfortify)
library(grid)
library(reshape2)
library(RColorBrewer)

getAspectRatio <- function(p){
  gb <- ggplot_build(p)
  g <- ggplot_gtable(gb)
  
  nullw <- sapply(g$widths, attr, "unit")
  nullh <- sapply(g$heights, attr, "unit")
  
  # ar of plot
  if(any(nullw == "null"))
    ar <- unlist(g$widths[nullw == "null"]) / unlist(g$heights[nullh == "null"])

  # ar of plot + legend
  g$fullwidth <- convertWidth(sum(g$widths), "in", valueOnly=TRUE)
  g$fullheight <- convertHeight(sum(g$heights), "in", valueOnly=TRUE)
  ar <- g$fullwidth / g$fullheight

  return(ar)
}

pcaPlots <- function(pca.data, pheno.data, meta.vars, title, ncol) {
  pheno.data[] <- lapply(pheno.data, as.character)
  plots <- c()
  ar <- -100
  for (i in meta.vars) {
    pl <- autoplot(pca.data, data = pheno.data, colour=i) +
                  coord_fixed()
    newar <- getAspectRatio(pl)
    if (ar<newar) {
      ar <- newar
    }
    plots <- c(plots, list(pl))
  }
  if(missing(ncol)) {
    pl <- plot_grid(plotlist = plots, ncol=length(meta.vars), align="hv")
  } else {
    pl <- plot_grid(plotlist = plots, ncol=ncol, align="hv")
  }
  if (!missing(title)) {
    title <- ggdraw() + draw_label(title, fontface='bold')  
    pl <- plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }
  pl <- pl + theme(plot.margin=margin(t=10, r=10, b=10, l=10))
  return(list(pl, ar))
}

illuPRatio <- function(df) {
  df$Ratio <- df$P95Grn / df$P05Grn
  hl <- df %>% group_by(Matrix) %>% summarize_at(vars(Ratio), mean)
  
  pl <- ggplot(df, aes(x=Section, y=Ratio)) + 
        geom_dotplot(binaxis='y') + ylab("P95 / P05") +
        geom_hline(aes(yintercept = Ratio), data=hl, linetype = "dashed", color='steelblue') +
        facet_grid(. ~ Matrix)
  
  return(pl)
}

geneRegr <- function(df, gene, var, group, x=median(df[,var], na.rm=T), coord=FALSE, maxy=0, miny=0) {
  df <- df[which(!outlier(df[,gene], logical=T)),]
  m <- lm(as.formula(paste(var, " ~ ", gene, sep="")), df);
  pl <- ggplot(df,aes_string(x=var,y=gene, colour=group)) + 
        geom_point() + geom_smooth(method=lm, linetype=1, colour="black", se=T, size=0.7) +
        #scale_colour_manual(values=c("green", "red", "blue")) +
        annotate("text", x=min(df[,var]+2.5, na.rm=T), y=max(df[,gene])-0.3, label=paste("R=",round(summary(m)$r.squared, digits = 3)), color="black", size=4) +
        geom_vline(aes(xintercept = x), linetype = "dashed", color='grey22') +
        theme_grey()
  
  if (coord) {
    pl <- pl + coord_cartesian(ylim=c(maxy,miny))
  }
  
  return(pl)
}

corPlot <- function(df, rep1, rep2, title="") {
  m <- lm(as.formula(paste(rep1, " ~ ", rep2, sep="")), df);
  pl <- ggplot(df,aes_string(x=rep1,y=rep2)) + 
    geom_point() + geom_smooth(method=lm, linetype=1, colour="black", se=F, size=0.7) +
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon', colour='black', alpha=0.5) + 
    scale_fill_continuous(low="#2FFF71", high="#FF4C48") +
    annotate("text", x=min(df[,rep1]+2.5, na.rm=T), y=max(df[,rep2])-0.3, label=paste("R=",round(summary(m)$r.squared, digits = 3)), color="black", size=4) +
    geom_vline(aes(xintercept = median(df[,rep1], na.rm=T)), linetype = "dashed", color='grey22') +
    geom_hline(aes(yintercept = median(df[,rep2], na.rm=T)), linetype = "dashed", color='grey22') +
    theme_grey() + ggtitle(title) + theme(plot.title = element_text(size=12))
  
  return(pl)
}

