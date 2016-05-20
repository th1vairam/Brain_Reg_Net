##' generate a heatmap + dendrograms, ggplot2 style
##'
##' @param x data matrix
##' @param hm.colours vector of colours (optional)
##' @return invisibly returns a list of ggplot2 objects. Display them with ggheatmap.show()
##' @author Chris Wallace
##' @export
##' @examples
##' ## test run
##' ## simulate data
##' library(mvtnorm)
##' sigma=matrix(0,10,10)
##' sigma[1:4,1:4] <- 0.6
##' sigma[6:10,6:10] <- 0.8
##' diag(sigma) <- 1
##' X <- rmvnorm(n=100,mean=rep(0,10),sigma=sigma)
##'  
##' ## make plot
##' p <- ggheatmap(X)
##'  
##' ## display plot
##' ggheatmap.show(p)
require(ggdendro)
require(ggplot2)
ggheatmap <- function(x,
                      hm.colours=my.colours) {
  # Sourced from https://raw.githubusercontent.com/chr1swallace/random-functions/master/R/ggplot-heatmap.R
  if(is.null(colnames(x)))
    colnames(x) <- sprintf("col%s",1:ncol(x))
  if(is.null(rownames(x)))
    rownames(x) <- sprintf("row%s",1:nrow(x))
  ## plot a heatmap
  ## x is an expression matrix
  row.hc <- hclust(dist(x), "ward")
  col.hc <- hclust(dist(t(x)), "ward")
  row.dendro <- ggdendro::dendro_data(as.dendrogram(row.hc),type="rectangle")
  col.dendro <- ggdendro::dendro_data(as.dendrogram(col.hc),type="rectangle")
  
  ## dendro plots
  col.plot <- mydplot(col.dendro, col=TRUE, labels=TRUE) +
    ggplot2::scale_x_continuous(breaks = 1:ncol(x),labels=col.hc$labels[col.hc$order]) +
    ggplot2::theme(plot.margin = unit(c(0,0,0,0), "lines"))
  row.plot <- mydplot(row.dendro, row=TRUE, labels=FALSE) +
    ggplot2::theme(plot.margin = unit(rep(0, 4), "lines"))
  
  ## order of the dendros
  col.ord <- match(col.dendro$labels$label, colnames(x))
  row.ord <- match(row.dendro$labels$label, rownames(x))
  xx <- x[row.ord,col.ord]
  dimnames(xx) <- NULL
  xx <- reshape2::melt(xx) %>%
    dplyr::rename(association = value)
  
  centre.plot <- ggplot2::ggplot(xx, aes(X2,X1)) + geom_tile(aes(fill=association), colour="white") +
    ggplot2::scale_fill_gradientn(colours = hm.colours) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0),breaks = NULL) +
    ggplot2::theme(plot.margin = unit(rep(0, 4), "lines"))
  ret <- list(col=col.plot,row=row.plot,centre=centre.plot)
  invisible(ret)
}