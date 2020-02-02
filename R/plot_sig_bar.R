#' Make bar plots of ligand/receptor expression and interaction scors
#'
#' This plot makes bar plots for ligand expression, receptor expression, raw interaction scores and normalized interaction scores for a specific ligand-receptor pair
#' @param ligand_exprs a numeric vector of ligand expression values across cell types
#' @param receptor_exprs a numeric vector of receptor expression values across cell types
#' @param cell_labels a character vector of cell type labels
#' @param filenames prefix for the name of the figures
#' @param raw_color color for the bar plot of raw interaction scores
#' @param norm_color color for the bar plot of normalized interaction scores
#' @param ligand_name name of the ligand
#' @param receptor_name name of the receptor
#' @param ylim maximum y axis for the normalized interaction scores.Default is 1
#' @return returns bar plots of ligand and receptor expression, as well as raw and normalized interaction scores
#' @examples plot_sig_bar(c(5,1,200),c(5,1,200),c("A","B","C"),"Figure4A_baseline",'slateblue3','slateblue1',"EGF","EGFR",1)
#' @export
plot_sig_bar<-function(ligand_exprs,receptor_exprs,cell_labels,filenames,raw_color,norm_color,ligand_name,receptor_name,y_lim=1){
  sig_raw<-NULL
  all_labels<-NULL
  for (i in 1:length(ligand_exprs)){
    sig_raw<-c(sig_raw,ligand_exprs[i]*receptor_exprs)
    all_labels<-c(all_labels,paste(cell_labels[i],"->",cell_labels,sep=""))
  }
  total_sig<-sum(ligand_exprs)*sum(receptor_exprs)
  norm_sig<-sig_raw/total_sig
  df<-data.frame(Expression=c(ligand_exprs,receptor_exprs),Cell.Type=rep(cell_labels,2),LR=rep(c(paste("Ligand (",ligand_name,")",sep = ""),paste("Receptor (",receptor_name,")",sep="")),each=length(ligand_exprs)))
  pdf(paste(filenames,"_ligand_receptor_exprs.pdf",sep=""),width=3.5,height=1.5)
  #print(ggplot(df,aes(x=Cell.Type,y=Expression)) + geom_bar(stat='identity',aes(fill=Cell.Type),width=0.5) + xlab("Cell type") + scale_fill_manual(values=c("salmon","yellowgreen","steelblue1","orange")) + facet_wrap(~LR,scales='free_y') + theme(legend.position = 'none') )
  #print(ggplot(df,aes(x=Cell.Type,y=Expression)) + geom_bar(stat='identity',aes(fill=Cell.Type),width=0.5) + xlab("Cell type") + scale_fill_manual(values=c("salmon","yellowgreen","steelblue1","orange")) + facet_wrap(~LR) + theme(legend.position = 'none') )
  print(ggplot2::ggplot(df,aes(x=Cell.Type,y=Expression)) + geom_bar(stat='identity',aes(fill=Cell.Type),width=0.5) + geom_text(aes(label=Expression,color=Cell.Type), vjust=-0.3, size=3.5) + xlab("Cell type") + facet_wrap(~LR) + theme(legend.position = 'none') + ylim(0,1.2*(max(df$Expression))) )
  #print(ggplot(df,aes(x=Cell.Type,y=Expression)) + geom_bar(stat='identity',aes(fill=Cell.Type),width=0.5)  + xlab("Cell type") + scale_fill_manual(values=c("salmon","yellowgreen","steelblue1","orange")) + facet_wrap(~LR,scales="free_y") + theme(legend.position = 'none'))
  dev.off()

  df<-data.frame(Interactions=all_labels,Raw.Scores=sig_raw)
  pdf(paste(filenames,"_raw.pdf",sep=""),width=3.5,height=1.5)
  print(ggplot2::ggplot(df,aes(x=Interactions,y=Raw.Scores)) + geom_bar(stat='identity',fill=raw_color,color=raw_color,width=0.5) + xlab("Cell-cell interaction") + theme_bw() + theme(axis.text.x = element_text(angle=-45),axis.text.y=element_text(color=raw_color,face='bold')))
  dev.off()

  df<-data.frame(Interactions=all_labels,Norm.Scores=norm_sig)
  bar_labels<-round(df$Norm.Scores,digits = 2)
  bar_labels<-ifelse(bar_labels<0.1,"",bar_labels)
  pdf(paste(filenames,"_norm.pdf",sep=""),width=3.5,height=1.5)
  print(ggplot2::ggplot(df,aes(x=Interactions,y=Norm.Scores)) + geom_bar(stat='identity',fill=norm_color,color=norm_color,width=0.5) + geom_text(aes(label=bar_labels),vjust=1.3, color="yellow3", fontface='bold',size=3.5)+ ylim(c(0,y_lim)) + xlab("Cell-cell interaction") + theme_bw() + theme(axis.text.x = element_text(angle=-45),axis.text.y=element_text(color=norm_color,face='bold')))
  dev.off()
}
