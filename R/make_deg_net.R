#'Prioritize ligand-receptor interactions based on differential expression
#'
#'A ligand-receptor pair is considerd when both the ligand and receptor genes are enriched in one cell type (e.g., expressed 3 fold higher in cell type A than each of the other cell types)
#'The resulting ligand-receptor network is then compared against talklr
#' @param expression_data expression matrix for ligand and receptor genes. 1st column "genes" (gene symbol, upper case); remaining columns: gene expression level
#' @param expressed_net the output from calling make_expressed_net.R function. To ensure the input network to both algorithms are the same.
#' @param fc_thresh fold change threshold above which a gene is considered enriched
#' @param pseudo_count a small value added to the expression matrices to avoid dividing by 0.
#' @return ligand-receptor network where both ligand and receptor genes are uniquely enriched in one cell type (i.e, marker genes for each cell type)
#' @export
make_deg_net<-function(expression_data,expressed_net,fc_thresh,pseudo_count){
  expression_data$genes<-toupper(expression_data$genes)
  n_cell<-ncol(expression_data)-1
  expressed_stat<-matrix(NA,nrow=nrow(expression_data),ncol=n_cell)
  for (i in 1:n_cell){
    mat<-matrix(rep(expression_data[,i+1],n_cell),nrow=nrow(expression_data),ncol=n_cell,byrow = F)
    expressed_stat[,i]<-rowSums(((mat+pseudo_count)/(expression_data[,2:ncol(expression_data)]+pseudo_count))>fc_thresh)>=(n_cell-1)
  }

  degs<- rowSums(expressed_stat)==1
  degs<-expression_data$genes[degs]
  deg_pairs<-expressed_net[expressed_net$Ligand.ApprovedSymbol %in% degs & expressed_net$Receptor.ApprovedSymbol %in% degs,]
  return(deg_pairs)
}
