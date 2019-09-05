#' Run VIPER iteratively
#'
#' This function runs VIPER iteratively
#'
#' @param datl Unified phosphoviper format
#' @return matrix
#' @import viper
#' @export
iterviper<-function(eset, regulon, iterations, ...) {

  for (it in 1:iterations) {
    vpres<-viper(eset, regulon, ...)
    eset<-rbind(vpres,eset[which(!(rownames(eset) %in% rownames(vpres))),])
  }

  vpres<-eset

  return(vpres)
}
