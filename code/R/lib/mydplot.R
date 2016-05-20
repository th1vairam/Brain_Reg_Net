require(ggplot2)
require(ggdendro)
mydplot <- function(ddata, row=!col, col=!row, labels=col) {
  ## plot a dendrogram
  yrange <- range(ddata$segments$y)
  yd <- yrange[2] - yrange[1]
  nc <- max(nchar(as.character(ddata$labels$label)))
  tangle <- if(row) { 0 } else { 90 }
  tshow <- col
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::labs(x = NULL, y = NULL) + ggdendro::theme_dendro()
  if(row) {
    p <- p +
      ggplot2::scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
      ggplot2::coord_flip()
  } else {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  return(p)
}