ggphylomorpho_size_HTC <- function (tree, tipinfo, 
                                    taxa = taxon, 
                                    xvar = PC1, 
                                    yvar = PC2, 
                                    factorvar = group, 
                                    labelvar = splabel, 
                                    point.size = size, 
                                    title.plot = "Phylomorphospace", 
                                    label.x.axis = "PC1", 
                                    label.y.axis = "PC2", 
                                    repel = TRUE, 
                                    edge.width = 1, 
                                    fontface = "italic", 
                                    tree.alpha = 0.7) 
{
  require(ggplot2)
  require(phytools)
  require(ggrepel)
  mat <- cbind(eval(substitute(xvar), tipinfo), 
               eval(substitute(yvar), tipinfo))
  rownames(mat) <- eval(substitute(taxa), tipinfo)
  stopifnot(length(setdiff(tree$tip.label, rownames(mat))) == 0)
  
  xAnc <- fastAnc(tree, mat[, 1])
  yAnc <- fastAnc(tree, mat[, 2])
  all_node_coords <- data.frame(x = c(mat[tree$tip.label, 1], xAnc), 
                                y = c(mat[tree$tip.label, 2], yAnc), 
                                nodeid = 1:(tree$Nnode + length(tree$tip.label)))
  edges <- data.frame(tree$edge)
  names(edges) <- c("node1", "node2")
  
  edgecoords <- merge(merge(edges, all_node_coords, by.x = "node1", 
                            by.y = "nodeid"), 
                      all_node_coords, by.x = "node2", by.y = "nodeid")
  
  pointsForPlot <- data.frame(x = eval(substitute(xvar), tipinfo), 
                              y = eval(substitute(yvar), tipinfo), 
                              color = eval(substitute(factorvar), tipinfo), 
                              label = eval(substitute(labelvar), tipinfo), 
                              size = eval(substitute(size), tipinfo))
  
  
  theplot <- ggplot() + 
    geom_segment(data = edgecoords, 
                 aes(x = x.x, xend = x.y, y = y.x, yend = y.y), 
                 linewidth = edge.width, 
                 alpha = tree.alpha) + 
    geom_point(data = pointsForPlot, 
               aes(x = x, y = y, color = color, size = size)) + 
    labs(title = title.plot, x = label.x.axis, y = label.y.axis) + 
    theme_bw(5) + 
    theme(legend.position = "bottom")
  
  if (repel){
    theplot <- theplot + 
      geom_text_repel(data = pointsForPlot, 
                      aes(x = x, y = y, label = label), 
                      segment.alpha = 0.5, 
                      fontface = fontface)
  }
  else {
    theplot <- theplot + 
      geom_text(data = pointsForPlot, 
                aes(x = x, y = y, label = label), 
                fontface = fontface)
  }
  return(theplot)
}
