#' Export hpARACNe files
#'
#' This function exports hpARACNe input files from the unified phosphoviper format
#'
#' @param datl Unified phosphoviper format
#' @param output_dir Output directory
#' @param kinases List of kinases (UniProtKB)
#' @param phosphatases List of phosphatases (UniProtKB)
#' @import data.table
#' @export
export2hparacne<-function(datl, output_dir, kinases, phosphatases) {
  # format output directory
  dir.create(output_dir)

  # reduce data to top peptide query per phosphosite
  pqp_freq<-data.table(merge(unique(datl[,c("site_id","peptide_id")]), as.data.frame(table(datl$peptide_id),stringsAsFactors = FALSE), by.x="peptide_id", by.y="Var1"))
  pqp_top<-pqp_freq[pqp_freq[, .I[which.max(Freq)], by=site_id]$V1]
  datl<-merge(datl, pqp_top[,c("site_id","peptide_id"),which=FALSE], by=c("site_id","peptide_id"))

  # generate matrix
  datm<-dcast(unique(datl[,c("peptide_id","run_id","peptide_intensity")]), peptide_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  phosphoAnno<-as.matrix(datm[,1])
  phosphoExp<-as.matrix(datm[,-1])
  phosphoExp[phosphoExp == 0]<-NA

  # write matrix
  write.table(cbind(phosphoAnno,phosphoExp), file=file.path(output_dir,"matrix.txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # write peptides
  write.table(unique(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite")]), file=file.path(output_dir,"peptides.txt"), quote=FALSE, row.names=FALSE, sep="\t")

  # write kinases
  write.table(unique(subset(datl, protein_id %in% kinases)$peptide_id), file=file.path(output_dir,"kinases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # write kinases + phosphatases
  write.table(unique(c(subset(datl, protein_id %in% kinases)$peptide_id, subset(datl, protein_id %in% phosphatases)$peptide_id)), file=file.path(output_dir,"kinases_phosphatases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}

#' Export to matrix
#'
#' This function exports tables in the unified phosphoviper format to a matrix
#'
#' @param datl Unified phosphoviper format
#' @return matrix
#' @import data.table
#' @export
export2mx<-function(datl) {
  # reduce data to top peptide query per phosphosite
  pqp_freq<-data.table(merge(unique(datl[,c("site_id","peptide_id")]), as.data.frame(table(datl$peptide_id),stringsAsFactors = FALSE), by.x="peptide_id", by.y="Var1"))
  pqp_top<-pqp_freq[pqp_freq[, .I[which.max(Freq)], by=site_id]$V1]
  datl<-merge(datl, pqp_top[,c("site_id","peptide_id"),which=FALSE], by=c("site_id","peptide_id"))

  # generate matrix
  datm<-dcast(unique(datl[,c("site_id","run_id","peptide_intensity")]), site_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  phosphoExp<-as.matrix(datm[,-1])
  phosphoExp[phosphoExp == 0]<-NA

  rownames(phosphoExp)<-datm$site_id

  return(phosphoExp)
}

#' Import hpARACNe network
#'
#' This function imports a hpARACNe network and generates site-specific regulons
#'
#' @param afile hpARACNe network file
#' @param pfile hpARACNe peptides file
#' @import data.table
#' @importFrom plyr dlply .
#' @export
hparacne2regulon<-function(afile, pfile) {
  aracne<-fread(afile)
  peptides<-fread(pfile)

  # Map peptide identifiers to site identififers
  aracne<-merge(aracne,peptides[,c("peptide_id","site_id"), with = FALSE], by.x="Regulator", by.y="peptide_id", allow.cartesian=TRUE)
  aracne<-merge(aracne,peptides[,c("peptide_id","site_id"), with = FALSE], by.x="Target", by.y="peptide_id", allow.cartesian=TRUE)
  aracne<-aracne[,c("site_id.x","site_id.y","MI","Correlation"), with = FALSE]
  names(aracne)<-c("Regulator","Target","MI","Correlation")

  # Compute likelihood from MI
  aracne$likelihood<-(aracne$MI / max(aracne$MI))

  # Generate regulons
  regulons<-dlply(aracne,.(Regulator),function(X){tfmode<-X$Correlation;names(tfmode)<-X$Target;return(list("tfmode"=tfmode, "likelihood"=X$likelihood))})

  return(regulons)
}

#' Convert site-specific regulons to protein regulons
#'
#' This function converts site-specific regulons to a list of protein - protein site-specific regulons
#'
#' @param regulons VIPER regulons
#' @import data.table
#' @importFrom plyr dlply .
#' @export
site2protein<-function(regulons) {
  # Generate regulon map
  regulon_map<-data.table("site_id"=names(regulons),"protein_id"=as.vector(sapply(names(regulons),function(X){strsplit(X,":")[[1]][2]})))

  # For each protein site, generate a unique regulon
  regulon_map[,regulon_id:=c(1:length(unique(site_id))),key="protein_id"]

  # Generate list of regulons
  regulons_list<-dlply(regulon_map,.(regulon_id),function(X){subregulons<-regulons[X$site_id];names(subregulons)<-X$protein_id;return(subregulons)})

  return(regulons_list)
}

#' Convert site-specific regulons to gene regulons
#'
#' This function converts site-specific regulons to a list of gene - gene regulons
#'
#' @param regulons VIPER regulons
#' @import data.table
#' @importFrom plyr dlply .
#' @export
site2gene<-function(regulons) {
  # Generate site_id map
  site_ids<-unique(c(names(regulons),unlist(sapply(regulons,function(X){names(X$tfmode)}))))
  id_map<-data.table("site_id"=site_ids,"gene_id"=as.vector(sapply(ids,function(X){strsplit(X,":")[[1]][1]})))

  # Generate regulon map
  regulon_map<-data.table("site_id"=names(regulons),"gene_id"=as.vector(sapply(names(regulons),function(X){strsplit(X,":")[[1]][1]})))

  # For each protein site, generate a unique regulon
  regulon_map[,regulon_id:=c(1:length(unique(site_id))),key="gene_id"]

  # Generate list of regulons
  regulons_list<-dlply(regulon_map,.(regulon_id),function(X){subregulons<-regulons[X$site_id];names(subregulons)<-X$gene_id;subregulons<-lapply(subregulons,function(Y){ids<-names(Y$tfmode);names(Y$tfmode)<-with(id_map, gene_id[match(ids,site_id)]);return(Y)});return(subregulons)})

  regulons_list_averaged<-lapply(regulons_list,function(X){lapply(X,function(Y){res<-sapply(unique(names(Y$tfmode)),function(Z){idx<-which(names(Y$tfmode)==Z);return(list("tfmode"=mean(Y$tfmode[idx]),"likelihood"=mean(Y$likelihood[idx])))});return(list("tfmode"=unlist(res[1,]),"likelihood"=unname(unlist(res[2,]))))})})

  return(regulons_list_averaged)
}

#' Subset regulons for targets
#'
#' This function subsets regulons to specified targets only
#'
#' @param regulons VIPER regulons
#' @param target Vector of targets
#' @export
subset_regulon<-function(regulons,targets) {
  subregulon<-lapply(regulons,function(X){ids<-which(names(X$tfmode) %in% targets);if(length(ids)>0){return(list("tfmode"=X$tfmode[ids], "likelihood"=X$likelihood[ids]))}})

  return(subregulon[sapply(subregulon,length)>0])
}
