#' Load human reference kinases
#'
#' Based on pkinfam version 2018_03
#'
#' @return List of human kinases in UniProtKB format
#' @references \url{https://www.uniprot.org/docs/pkinfam}
#' @export
loadHumanKinases<-function() {
  # kinases<-readLines("~/Desktop/pv_test/human_kinases.txt")
  # save(kinases,file="data/human_kinases.rda")
  data("human_kinases")
}

#' Load human reference phosphatases
#'
#' Based on Li X, Wilmanns M, Thornton J, Köhn M. Sci Signal. 2013 May 14;6(275):rs10. doi: 10.1126/scisignal.2003203.
#'
#' @return List of human phosphatases in UniProtKB format
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/23674824}
#' @export
loadHumanPhosphatases<-function() {
  # phosphatases<-readLines("~/Desktop/pv_test/human_phosphatases.txt")
  # save(phosphatases,file="data/human_phosphatases.rda")
  data("human_phosphatases")
}

#' Load human reference HSM/P interactions
#'
#' Based on Cunningham JM, Koytiger G, Sorger PK, AlQuraishi M. Nat Methods. 2020 Feb;17(2):175-183. doi: 10.1038/s41592-019-0687-1.
#'
#' @return Table of PPI with their scores
#' @references \url{https://figshare.com/articles/Predictions_-_Domain-Peptide_and_Protein-Protein_Interactions_-_HSM/10084745/2}
#' @export
loadHSMP<-function() {
  # hsmp<-unique(fread("hsm_p.predictions.csv.gz",header=FALSE,col.names=c("regulator","target","confidence")))
  # hsmp2<-hsmp
  # names(hsmp2)<-c("target","regulator","confidence")
  # hsmp<-unique(rbind(hsmp,hsmp2))
  # hsmd<-unique(readRDS("hsmd.rds")[,c("duid","dtype")])
  # names(hsmd)<-c("regulator","regulator_type")
  # hsmd<-ddply(hsmd,.(regulator),function(X){data.frame("regulator_type"=paste(X$regulator_type,collapse=";"))})
  # hsmp<-merge(hsmp, hsmd, by=c("regulator"))
  # hsmp[,confidence:=max(confidence),by=c("regulator","target")]
  # hsmp<-unique(hsmp)
  # save(hsmp,file="../../code/phosphoviper/data/human_hsmp.rda")
  data("human_hsmp")
}

#' Load human reference HSM/P interactions
#'
#' Based on Cunningham JM, Koytiger G, Sorger PK, AlQuraishi M. Nat Methods. 2020 Feb;17(2):175-183. doi: 10.1038/s41592-019-0687-1.
#'
#' @return Table of PPI with their scores
#' @references \url{https://figshare.com/articles/Predictions_-_Domain-Peptide_and_Protein-Protein_Interactions_-_HSM/10084745/2}
#' @export
loadHSMPF<-function() {
  # hsmpf<-unique(fread("hsm_p.predictions.csv.gz",header=FALSE,col.names=c("regulator","target","confidence")))
  # hsmpf2<-hsmpf
  # names(hsmpf2)<-c("target","regulator","confidence")
  # hsmpf<-unique(rbind(hsmpf,hsmpf2))
  # hsmpf[,confidence:=max(confidence),by=c("regulator","target")]
  # hsmpf<-unique(hsmpf)
  # save(hsmpf,file="../../code/phosphoviper/data/human_hsmpf.rda")
  data("human_hsmpf")
}

#' Load human reference HSM/D interactions
#'
#' Based on Cunningham JM, Koytiger G, Sorger PK, AlQuraishi M. Nat Methods. 2020 Feb;17(2):175-183. doi: 10.1038/s41592-019-0687-1.
#'
#' @return Table of DPI with their scores
#' @references \url{https://figshare.com/articles/Predictions_-_Domain-Peptide_and_Protein-Protein_Interactions_-_HSM/10084745/2}
#' @export
loadHSMD<-function() {
  # hsmd<-readRDS("hsmd.rds")[,c("duid","dtype","puid","psid","val")]
  # names(hsmd)<-c("regulator","regulator_type","target","site_id","confidence")
  # hsmd[,confidence:=max(confidence),by=c("regulator","regulator_type","target","site_id")]
  # hsmd<-unique(hsmd)
  # save(hsmd,file="../../code/phosphoviper/data/human_hsmd.rda")
  data("human_hsmd")
}

#' Load human reference PathwayCommons interactions
#'
#' Based on PathwayCommons, version 12.
#'
#' @return Table of PPI with their scores
#' @references \url{https://www.pathwaycommons.org/}
#' @export
loadPCdb<-function() {
  # pcdb<-readRDS("pc.rds")
  # save(pcdb,file="../../code/phosphoviper/data/human_pcdb.rda")
  data("human_pcdb")
}
