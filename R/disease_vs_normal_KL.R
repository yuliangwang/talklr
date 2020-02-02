#'Identify perturbed ligand-receptor interactions between two conditions
#'
#' This function calculates Kullback-Leibler divergence of ligand-receptor interaction scores between two conditions (e.g., disease vs. normal)
#' @param ligand_mat_disease ligand expression matrix in the disease/perturbed condition
#' @param receptor_mat_disease receptor expression matrix in the disease/perturbed condition. The same ligand-receptor pair should be in the same order as ligand_mat_disease
#' @param ligand_mat_normal ligand expression matrix in the normal/baseline condition. Same order as ligand_mat_disease
#' @param receptor_mat_normal receptor expression matrix in the normal/baseline coniditon. The same ligand-receptor pair should be in the same order as ligand_mat_normal
#' @param pseudo_count a small value added to the expression matrices to avoid dividing by 0.
#' @param method "product" or "min": ligand-receptor interaction score is either defined by multiplying ligand and receptor expression values or take the minimum
#' @return KL divergence for the ligand-receptor interaction scores between disease vs. normal (or perturbed vs. baseline) conditions
#' @export
disease_vs_normal_KL<-function(ligand_mat_disease,receptor_mat_disease,ligand_mat_normal,receptor_mat_normal,pseudo_count,method){
  ligand_mat_diease<-as.matrix(ligand_mat_disease)
  receptor_mat_disease<-as.matrix(receptor_mat_disease)
  ligand_mat_normal<-as.matrix(ligand_mat_normal)
  receptor_mat_normal<-as.matrix(receptor_mat_normal)
  ligand_mat_disease<-ligand_mat_disease+pseudo_count
  ligand_mat_normal<-ligand_mat_normal+pseudo_count
  receptor_mat_disease<-receptor_mat_disease+pseudo_count
  receptor_mat_normal<-receptor_mat_normal+pseudo_count

  n<-ncol(ligand_mat_disease)
  p<-nrow(ligand_mat_disease)
  interaction_mat_disease<-matrix(0,nrow=nrow(ligand_mat_disease),ncol=1)
  interaction_mat_normal<-interaction_mat_disease
  if (method=="product"){
    for (i in 1:n){
      temp_mat<-matrix(rep(ligand_mat_disease[,i],n),nrow=p,ncol=n,byrow=F)
      interaction_mat_disease<-cbind(interaction_mat_disease,temp_mat*receptor_mat_disease)

      temp_mat<-matrix(rep(ligand_mat_normal[,i],n),nrow=p,ncol=n,byrow=F)
      interaction_mat_normal<-cbind(interaction_mat_normal,temp_mat*receptor_mat_normal)
    }
  } else {
    for (i in 1:n){
      for (j in 1:n){
        interaction_mat_disease<-cbind(interaction_mat_disease,apply(cbind(ligand_mat_disease[,i],receptor_mat_disease[,j]),1,min))
        interaction_mat_normal<-cbind(interaction_mat_normal,apply(cbind(ligand_mat_normal[,i],receptor_mat_normal[,j]),1,min))
      }
    }
  }
  interaction_mat_disease<-interaction_mat_disease[,-1]
  interaction_mat_normal<-interaction_mat_normal[,-1]
  interaction_mat_disease<- interaction_mat_disease/matrix(rep(rowSums(interaction_mat_disease),n^2),nrow=p,ncol=n^2,byrow=F)
  interaction_mat_normal<- interaction_mat_normal/matrix(rep(rowSums(interaction_mat_normal),n^2),nrow=p,ncol=n^2,byrow=F)
  interaction_mat<-interaction_mat_disease*log2(interaction_mat_disease/interaction_mat_normal)
  interaction_mat[is.nan(interaction_mat)]<-0
  KL<-rowSums(interaction_mat)
  return(KL)
}

