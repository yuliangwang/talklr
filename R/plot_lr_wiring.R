#' Plot ligand-receptor wiring diagram
#'
#' This function uses the iGraph package to plot the ligand-receptor wiring diagram among interacting cell types
#' @param ligand_exprs ligand expression matrix
#' @param receptor_exprs receptor expression matrix
#' @param cell_labels character vector of cell type names
#' @return a plot of the ligand-receptor wiring diagram
#' @export
plot_lr_wiring<-function(ligand_exprs,receptor_exprs,cell_labels){
  n_cell<-length(cell_labels)
  norm_mat<-matrix(0,n_cell,n_cell)
  ligand_exprs[ligand_exprs<1]<- 0 #if not expressed, set to 0
  receptor_exprs[receptor_exprs<1]<- 0
  for (i in 1:n_cell){
    norm_mat[i,]<-ligand_exprs[i]*receptor_exprs
  }
  # min_val<-min(norm_mat)
  # max_val<-max(norm_mat)
  # final_mat<-(norm_mat - min_val)/(max_val-min_val)
  total<-sum(norm_mat)
  if (total>0){
    final_mat<-norm_mat/total
    final_mat[final_mat<0.1]<-0
  }
  else {
    final_mat<-norm_mat
  }
  rownames(final_mat)<-cell_labels
  colnames(final_mat)<-cell_labels
  require(igraph)
  net<-graph_from_adjacency_matrix(final_mat,mode="directed",weighted=TRUE)
  E(net)$width <- 4*(E(net)$weight)
  co<-layout_in_circle(net)

  edgeList <- do.call(rbind, igraph::get.adjedgelist(net))
  if (sum(E(net)$weight)>0){
    par(mar = c(0.3, 0.3, 0.3, 0.3))
    plot(net,edge.label=round(E(net)$weight,digits = 2),edge.curved=TRUE,vertex.label.font=2,edge.arrow.size=0.6,layout=co,margin=c(0.4,0.4,0.4,0.4),vertex.size=25,vertex.label.cex=2,vertex.label.dist=4,vertex.label.degree=-pi/2,edge.color='black',vertex.color='white',edge.label.color='brown',edge.label.cex=1.5,edge.label.font=2)
  }
  else{
    par(mar = c(0.3, 0.3, 0.3, 0.3))
    plot(net,lty=0,vertex.label.font=2,layout=co,margin=c(0.4,0.4,0.4,0.4),vertex.size=25,vertex.label.cex=2,vertex.label.dist=4,vertex.label.degree=-pi/2,vertex.color='white',edge.label.color='brown',edge.label.cex=1.5,edge.label.font=2)
  }
}
