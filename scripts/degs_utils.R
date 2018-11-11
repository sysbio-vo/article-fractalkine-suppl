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

