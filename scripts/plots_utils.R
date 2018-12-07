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

ggloadings <- function(pca, num_pca=1, num_loading=5, xlab="genes") {
  loadings <- pca$rotation[,num_pca]
  loadings <- abs(loadings)
  loadings <- sort(loadings, decreasing = TRUE)
  label_names <- names(loadings)
  loadings <- loadings[1:num_loading]
  label_names <- label_names[1:num_loading]
  dat <- data.frame(pc = label_names, loadings = loadings)
  dat$pc <- factor(dat$pc, levels = unique(dat$pc))
  
  p <- ggplot(dat, aes(x = pc, y = loadings))
  p <- p + geom_bar(stat = "identity")
  p <- p + xlab(xlab) + ylab("contribution scores")
  return(p)
}

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

GOHeat <- function (data, nlfc, fill.col) 
{
  x <- y <- z <- NULL
  if (missing(nlfc)) 
    nlfc <- 0
  else nlfc <- nlfc
  if (missing(fill.col)) 
    fill.col <- c("firebrick", "white", "dodgerblue")
  else fill.col <- fill.col
  distance <- dist(data)
  cluster <- hclust(distance)
  M <- dim(data)[2]
  nterm <- M - nlfc
  if (nlfc == 0) {
    s <- rowSums(data[, 1:nterm])
    tmp <- NULL
    for (r in 1:nrow(data)) {
      tmp <- c(tmp, as.numeric(gsub(1, s[r], data[r, 1:nterm])))
    }
  }
  else {
    tmp <- NULL
    for (r in 1:nrow(data)) {
      tmp <- c(tmp, as.numeric(gsub(1, data[r, (nterm + 
                                                  1)], data[r, 1:nterm])))
    }
  }
  df <- data.frame(x = factor(rep(cluster$order, each = nterm)), y = rep(colnames(data[, 
                                                                                       1:nterm]), length(rownames(data))), z = tmp, lab = rep(rownames(data), 
                                                                                                                                              each = nterm))
  df_o <- df[order(df$x), ]
  g <- ggplot() +
    geom_tile(data = df_o, aes(x = x, y = y, fill = z), colour="white", size=0.25) +
    scale_x_discrete(breaks = 1:length(unique(df_o$x)), labels = unique(df_o$lab)) +
    theme(axis.text.x = element_text(size = 11, angle=60, hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 12),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed()
  if (nlfc == 0) {
    g + scale_fill_gradient2("Count", space = "Lab", low = fill.col[2], 
                             mid = fill.col[3], high = fill.col[1])
  }
  else {
    g + scale_fill_gradient2("FC", space = "Lab", low = fill.col[3], 
                             mid = fill.col[2], high = fill.col[1],
                             limits=c(round(min(df_o$z), 2)-0.1, round(max(df_o$z), 2)+0.1),
                             breaks = c(round(min(df_o$z), 2), 0, round(max(df_o$z), 2)))
  }
}
