#' Benchmark hpARACNe network against PathwayCommons
#'
#' This function conducts a benchmark of a hpARACNe network against PathwayCommons
#'
#' @param network hpARACNe network
#' @param pwfile PathwayCommons SIF file
#' @param negfile Negative interactome file
#' @import data.table
#' @import pROC
#' @export
hparacne2benchmark<-function(network, peptides, pwfile, samplefactor=1) {
  # Parse PC database
  pcd<-fread(pwfile, header=FALSE, col.names=c("regulator","moi","target"))
  pcd<-subset(pcd, moi=="controls-phosphorylation-of")[,c("regulator","target")]
  pcd$GT<-1

  # Prepare network
  net<-fread(network, col.names=c("regulator","target","mi","correlation","pvalue"))
  pep<-fread(peptides)[,c("peptide_id","gene_id")]
  net<-merge(net, pep, by.x="regulator", by.y="peptide_id", allow.cartesian=TRUE)[,c("gene_id","target","mi","correlation","pvalue")]
  names(net)<-c("regulator","target","mi","correlation","pvalue")
  net<-merge(net, pep, by.x="target", by.y="peptide_id", allow.cartesian=TRUE)[,c("regulator","gene_id","mi","correlation","pvalue")]
  names(net)<-c("regulator","target","mi","correlation","pvalue")
  net$regulon<-1

  # Prepare random network
  targets<-unique(fread(peptides)$gene_id)
  regulators<-unique(fread(peptides)[,c("gene_id","protein_id")])
  regulators<-regulators[which(regulators$protein_id %in% c(kinases, phosphatases)),]$gene_id
  combnet<-expand.grid("regulator"=regulators, "target"=targets)

  posnet<-unique(merge(net, pcd, by=c("regulator","target"))[,c("regulator","target")])
  posnet$GT<-1
  negnet<-unique(merge(combnet, pcd, by=c("regulator","target"), all.x=TRUE))
  negnet<-subset(negnet, is.na(GT))
  negnet$GT<-0
  negnet<-negnet[sample(1:dim(negnet)[1], dim(posnet)[1]*samplefactor),]

  # Prepare ground truth
  reference<-rbind(posnet, negnet)
  bnet<-merge(net, reference, by=c("regulator","target"))

  rocl<-dlply(bnet, .(regulon), function(X){roc(X, response="GT", predictor="mi", direction=">")})

  class(rocl)<-"rocl"

  return(rocl)
}

#' Benchmark regulons against PathwayCommons
#'
#' This function conducts a benchmark of gene regulons against PathwayCommons
#'
#' @param regulons Gene-level VIPER regulons
#' @param pwfile PathwayCommons SIF file
#' @import data.table
#' @import pROC
#' @importFrom reshape expand.grid.df
#' @export
regulons2benchmark<-function(regulons, pwfile) {
  # Parse PC database
  pcd<-fread(pwfile, header=FALSE, col.names=c("regulator","moi","target"))
  pcd<-subset(pcd, moi=="controls-phosphorylation-of")[,c("regulator","target")]
  pcd$GT<-1

  # Transform regulons to network
  # If there is only one regulon per regulator
  if ("tfmode" %in% names(regulons)) {
    net<-ldply(regulons, function(X){return(data.frame("target"=names(X$tfmode),"likelihood"=X$likelihood))})
    names(net)<-c("regulator","target","likelihood")
    net$regulon<-1
    net<-net[,c("regulon","regulator","target","likelihood")]
  }
  # If there are multiple regulons (metaVIPER)
  else {
    if (is.null(names(regulons))) {
      names(regulons)<-1:length(regulons)
    }
    net<-ldply(regulons, function(Y){Z<-ldply(Y, function(X){return(data.frame("target"=names(X$tfmode),"likelihood"=X$likelihood))}); names(Z)<-c("regulator","target","likelihood"); return(Z)})
    names(net)<-c("regulon","regulator","target","likelihood")

    # Expand grid
    net<-merge(net, expand.grid.df(as.data.frame(unique(net[,c("regulator","target")])), data.frame("regulon"=unique(net[,c("regulon")]))), by=c("regulon","regulator","target"), all=TRUE)
    net[which(is.na(net$likelihood)),"likelihood"]<-0
  }

  # Prepare random network
  negnet<-data.frame("regulator"=sample(pcd$regulator), "target"=sample(pcd$target))
  negnet<-merge(negnet, pcd, by=c("regulator","target"), all.x=TRUE)
  negnet[which(is.na(negnet$GT)),"GT"]<-0
  negnet<-subset(negnet, GT==0)

  # Prepare ground truth
  reference<-rbind(pcd, negnet)

  bnet<-merge(net, reference, by=c("regulator","target"))

  rocl<-dlply(bnet, .(regulon), function(X){roc(X, response="GT", predictor="likelihood", direction="<")})

  class(rocl)<-"rocl"

  return(rocl)
}

#' Benchmark regulons against PathwayCommons
#'
#' Print summary statistics for benchmark
#'
#' @param rocl Benchmark object
#' @import pROC
#' @export
summary.rocl<-function(rocl) {
  sapply(rocl,print)
}

#' Benchmark regulons against PathwayCommons
#'
#' Plot ROC curves
#'
#' @param rocl Benchmark object
#' @import pROC
#' @export
plot.rocl<-function(rocl) {
  colors<-rainbow(length(rocl))
  for (i in 1:length(rocl)) {
    if (i==1) {
      plot(rocl[[i]], add=FALSE, col=colors[i])
    }
    else {
      plot(rocl[[i]], add=TRUE, col=colors[i])
    }
    legend("bottomright", legend=names(rocl), col=colors, lty=1)
  }
}
