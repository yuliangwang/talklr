#'Prioritize ligand-receptor interactions in a single condition
#'
#'Instead of using an expression threshold, this function identifies interesting ligand-receptor pairs among a specific list of expressed genes provided by the user.
#'Run this function separately for perturbed and baseline conditions before running disease_vs_normal_KL.R, so that ligand and receptor matrices are in the same order.
#' @param expression_data expression matrix for ligand and receptor genes. 1st column "genes" (gene symbol, upper case); remaining columns: gene expression level
#' @param expressed_genes a character vector of expressed genes
#' @param receptor_ligand receptor-ligand interaction database from Ramilowski JA et al, Nat. Comm. 2015. Provided with the package
#' @param KL_method "product" or "min": ligand-receptor interaction score is either defined by multiplying ligand and receptor expression values or take the minimum
#' @return KL divergence for all expressed ligand-receptor pairs, as well as ligand and receptor expression values across cell types. You may want to sort from highest to lowest to examine the top ones.
#' @export
make_expressed_net_specify_expressed_genes<-function(expression_data,expressed_genes,receptor_ligand,KL_method){
  #expression_data: 1st column "genes" (gene symbol, upper case); remaining columns: gene expression level
  expression_data$genes<-toupper(expression_data$genes)
  n_cell<-ncol(expression_data)-1
  expression_data<-expression_data[expression_data$genes %in% expressed_genes,]

  expressed_net<- receptor_ligand[receptor_ligand$Ligand.ApprovedSymbol %in% expressed_genes & receptor_ligand$Receptor.ApprovedSymbol %in% expressed_genes,]
  ind<-match(expressed_net$Ligand.ApprovedSymbol,expression_data$genes)
  expressed_net<-cbind(expressed_net,expression_data[ind,2:ncol(expression_data)])
  colnames(expressed_net)[17:(17+n_cell-1)]<-paste("ligand_",colnames(expression_data)[2:ncol(expression_data)],sep="")

  ind<-match(expressed_net$Receptor.ApprovedSymbol,expression_data$genes)
  expressed_net<-cbind(expressed_net,expression_data[ind,2:ncol(expression_data)])
  colnames(expressed_net)[(17+n_cell):(17+2*n_cell-1)]<-paste("receptor_",colnames(expression_data)[2:ncol(expression_data)],sep="")

  expressed_net$KL<-pairwise.interaction.KL(as.matrix(expressed_net[,17:(17+n_cell-1)]),as.matrix(expressed_net[,(17+n_cell):(17+2*n_cell-1)]),method=KL_method)
  return(expressed_net)
}
