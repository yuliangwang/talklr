#' Make scatter plots of GO enrichment results
#'
#' This function makes scatter plots of GO enrichment results. Diagonal line denotes equal enrichment.
#' GO terms with higher enrichment in talklr pairs are in green; GO terms with higher enrichment in DEG approach are in orange
#' Axis are drawn to scale. Longer axis means larger values.
#'
#' @param go_deg enrichment table for ligand-receptor pairs identified by the DEG approach. Two required columns are GO.ID and FDR
#' @param go_product enrichment table for ligand-receptor pairs identified by talklr
#' @param interesting a character vector of GO terms to be visualized
#' @param fdr_thresh FDR threshold above which a GO is considered significant
#' @param diff_thresh if differences of -log10(FDR) between the two appraoches are larger than this threshold, the GO terms will be labeled.
#' @param filename name of the file where the scatter plot will be saved.
#' @return A scatter plot named filename
#' @export
plot_GO_comparison<-function(go_deg,go_product,interesting,fdr_thresh,diff_thresh,filename){
  require(ggplot2)
  require(ggrepel)
  all_gos<-rbind(go_deg,go_product)
  all_gos<-all_gos[!duplicated(all_gos$GO.ID),]
  ind<-match(interesting,all_gos$Term)
  toplot<-data.frame(id=all_gos$GO.ID[ind],name=interesting,fdr_product=1,fdr_deg=1)
  ind<-match(toplot$name,go_product$Term)
  toplot$fdr_product[!is.na(ind)]<-go_product$FDR[ind[!is.na(ind)]]
  ind<-match(toplot$name,go_deg$Term)
  toplot$fdr_deg[!is.na(ind)]<-go_deg$FDR[ind[!is.na(ind)]]

  toplot<-arrange(toplot,fdr_product)
  toplot<-toplot[toplot$fdr_product<fdr_thresh | toplot$fdr_deg<fdr_thresh,]
  p1<-ggplot2::ggplot(toplot,aes(x=-log10(fdr_deg),y=-log10(fdr_product))) + geom_point() + geom_abline(slope=1,lty=2) + geom_text_repel(data=toplot[abs(-log10(toplot$fdr_product) + log10(toplot$fdr_deg))>diff_thresh ,],aes(label=name),size=3,force=4,segment.alpha=0.4) + coord_equal() + ylab("-log10(talkr pairs)") + xlab("-log10(DEG pairs)")

  pdf(filename,width = 5,height=5)
  print(p1)
  dev.off()
  return(p1)
}
