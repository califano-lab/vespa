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
