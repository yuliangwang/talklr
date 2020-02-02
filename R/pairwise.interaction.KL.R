#'Calulate KL divergence between observed vs. expected interaction score distributions.
#'
#'Expected interaction score distribution is that all cell types express the same level of ligand and receptor genes, thus the normalized probability is 1/n^2.
#'Observed distribution is the interaction scores based on actual expression values.
#'This is a helper function called by make_expressed_net.R or make_expressed_net_specify_expressed_genes.R
#' @param ligand_mat expression matrix for ligands
#' @param receptor_mat expression matrix for receptors. The same ligand-receptor pair should be in the same order as ligand_mat.
#' @param method  "product" or "min": ligand-receptor interaction score is either defined by multiplying ligand and receptor expression values or take the minimum
#' @return KL divergence for all expressed ligand-receptor pairs.
#' @export
pairwise.interaction.KL<-function(ligand_mat,receptor_mat,method){
  n<-ncol(ligand_mat)
  p<-nrow(ligand_mat)
  ligand_mat<-ligand_mat
  receptor_mat<-receptor_mat
  interaction_mat<-matrix(0,nrow=nrow(ligand_mat),ncol=1)
  if (method=="product"){
    for (i in 1:n){
      temp_mat<-matrix(rep(ligand_mat[,i],n),nrow=p,ncol=n,byrow=F)
      interaction_mat<-cbind(interaction_mat,temp_mat*receptor_mat)
    }
  } else {
    for (i in 1:n){
      for (j in 1:n){
        interaction_mat<-cbind(interaction_mat,apply(cbind(ligand_mat[,i],receptor_mat[,j]),1,min))
      }
    }
  }
  interaction_mat<-interaction_mat[,-1]
  interaction_mat<- interaction_mat/matrix(rep(rowSums(interaction_mat),n^2),nrow=p,ncol=n^2,byrow=F)
  interaction_mat<-interaction_mat*log2(interaction_mat*(n^2))
  interaction_mat[is.nan(interaction_mat)]<-0
  KL<-rowSums(interaction_mat)
  return(KL)
}

