#' Export hpARACNe files
#'
#' This function exports hpARACNe input files from the unified phosphoviper format
#'
#' @param datl Unified phosphoviper format
#' @param output_dir Output directory
#' @param kinases List of kinases (UniProtKB)
#' @param phosphatases List of phosphatases (UniProtKB)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @export
exportARACNe<-function(datl, output_dir, kinases, phosphatases) {
  # format output directory
  dir.create(output_dir)

  # reduce data to top peptide query per phosphosite
  pqp_freq<-data.table(merge(unique(datl[,c("site_id","peptide_id")]), as.data.frame(table(datl$peptide_id),stringsAsFactors = FALSE), by.x="peptide_id", by.y="Var1"))
  pqp_top<-pqp_freq[pqp_freq[, .I[which.max(Freq)], by=site_id]$V1]$peptide_id
  datl<-subset(datl, peptide_id %in% pqp_top)

  # generate matrix
  datm<-dcast(unique(datl[,c("site_id","run_id","peptide_intensity")]), site_id ~ run_id, value.var = "peptide_intensity", fun.aggregate=mean)

  # get annotation
  phosphoAnno<-as.matrix(datm[,1])
  phosphoExp<-as.matrix(datm[,-1])
  phosphoExp[phosphoExp == 0]<-NA

  # write matrix
  write.table(cbind(phosphoAnno,phosphoExp), file=file.path(output_dir,"matrix.txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # write peptides
  write.table(unique(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite")]), file=file.path(output_dir,"peptides.txt"), quote=FALSE, row.names=FALSE, sep="\t")

  # write kinases
  write.table(unique(subset(datl, protein_id %in% kinases)$site_id), file=file.path(output_dir,"kinases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # write kinases + phosphatases
  write.table(unique(c(subset(datl, protein_id %in% kinases)$site_id, subset(datl, protein_id %in% phosphatases)$site_id)), file=file.path(output_dir,"kinases_phosphatases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}
