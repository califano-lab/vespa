#' Export hpARACNe files
#'
#' This function exports hpARACNe input files from the unified phosphoviper format
#'
#' @param datl Unified phosphoviper format
#' @param output_dir Output directory
#' @param kinases List of kinases (UniProtKB)
#' @param phosphatases List of phosphatases (UniProtKB)
#' @param interactions HSM/D or HSM/P Interaction table (UniProtKB)
#' @param confidence_threshold Interaction confidence_threshold
#' @import data.table
#' @export
export2hparacne<-function(datl, output_dir, kinases, phosphatases, interactions, confidence_threshold=0.5) {
  # format output directory
  dir.create(output_dir)

  # reduce data to top peptide query per phosphosite
  pqp_freq<-data.table(merge(unique(datl[,c("site_id","peptide_id")]), as.data.frame(table(datl$peptide_id),stringsAsFactors = FALSE), by.x="peptide_id", by.y="Var1"))
  pqp_top<-pqp_freq[pqp_freq[, .I[which.max(Freq)], by=site_id]$V1]
  datl<-merge(datl, pqp_top[,c("site_id","peptide_id"),which=FALSE], by=c("site_id","peptide_id"))

  # generate matrix
  datm<-dcast(data.table(unique(datl[,c("peptide_id","run_id","peptide_intensity")])), peptide_id ~ run_id, value.var = "peptide_intensity")

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

  # write targets (only phosphopeptides or VIPER activity is used for targets, but not protein abundance)
  write.table(unique(subset(datl, phosphosite != "PA")$peptide_id), file=file.path(output_dir,"targets.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # subset to regulators
  # regulator_ppi<-subset(regulator_ppi, regulator %in% c(kinases, phosphatases))

  if ("site_id" %in% names(interactions)) {
    regulator_ppi<-unique(datl[,c("protein_id","peptide_id")])
    names(regulator_ppi)<-c("regulator","regulator_peptide_id")

    target_ppi<-unique(datl[,c("site_id","peptide_id")])
    names(target_ppi)<-c("site_id","target_peptide_id")

    phosphointeractions<-merge(merge(subset(interactions, confidence > confidence_threshold), regulator_ppi, by="regulator", allow.cartesian=TRUE), target_ppi, by="site_id", allow.cartesian=TRUE)[,c("regulator_peptide_id","target_peptide_id","confidence")]
  }
  else {
    regulator_ppi<-unique(subset(datl, phosphosite=="PA")[,c("protein_id","peptide_id")])
    names(regulator_ppi)<-c("regulator","regulator_peptide_id")

    target_ppi<-unique(subset(datl, phosphosite=="PV")[,c("protein_id","peptide_id")])
    names(target_ppi)<-c("target","target_peptide_id")

    phosphointeractions<-unique(merge(merge(subset(interactions, confidence > confidence_threshold), regulator_ppi, by="regulator", allow.cartesian=TRUE), target_ppi, by="target", allow.cartesian=TRUE)[,c("regulator_peptide_id","target_peptide_id","confidence")])
  }

  write.table(phosphointeractions, file=file.path(output_dir,"phosphointeractions.txt"), quote=FALSE, row.names=FALSE, sep="\t")
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

#' Export to unified phosphoviper format from viper matrix
#'
#' This function exports a viper matrix to the unified phosphoviper format
#'
#' @param matrix
#' @param fasta Amino acid FASTA file from UniProt
#' @param tag Site tag
#' @return phosphoviper table
#' @import data.table
#' @export
vmx2pv<-function(vmx, fasta, tag = "PV") {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))

  pv<-reshape2::melt(vmx, value.name="peptide_intensity")
  names(pv)<-c("protein_id","run_id","peptide_intensity")
  pv<-merge(genes_proteins, pv, by="protein_id")
  pv$modified_peptide_sequence<-pv$protein_id
  pv$peptide_sequence<-pv$protein_id
  pv$phosphosite<-tag
  pv$site_id<-paste(pv$gene_id,pv$protein_id,pv$phosphosite,sep=":")
  pv$peptide_id<-paste(pv$gene_id,pv$protein_id,pv$phosphosite,sep=":")

  return(pv)
}

#' metaVIPER regulon selection
#'
#' @param ges Gene Expression Signature (features X samples)
#' @param net.list List object with the networks to be used
#' @return Optimized network with best regulon per gene
#' @export
optimizeRegulon <- function(ges, net.list, min_size=25, pleiotropy=FALSE, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10, penalty = 20, method = "adaptive")) {
  selectRegulon<-function(vip.list, net.list, g) {
    w.mat <- matrix(0L, nrow = num.nets, ncol = 1)
    rownames(w.mat) <- names(net.list)

    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(0)
        }}))
      max.ind <- which.max(abs(nes.vals))
      w.mat[max.ind, 1] <- w.mat[max.ind, 1] + 1
    }

    return(net.list[[row.names(w.mat)[which.max(w.mat)]]][[g]])
  }

  num.nets <- length(net.list)
  num.samps <- ncol(ges)

  ## run VIPER with each network
  print('Generating VIPER matrices...')
  vip.list <- list()
  for (i in 1:num.nets) {
    vip.list[[i]] <- viper(ges, net.list[i], method = 'none', minsize=min_size, pleiotropy = pleiotropy, pleiotropyArgs
 = pleiotropyArgs)
  }
  names(vip.list) <- names(net.list)

  ## select best regulon for each gene
  message("")
  message('Selecting best regulons by max(abs(NES))...')
  uni.genes <- unique(unlist(lapply(vip.list, rownames)))
  regulons<-lapply(uni.genes, function(g){selectRegulon(vip.list, net.list, g)})
  names(regulons)<-uni.genes

  regulons[sapply(regulons, is.null)] <- NULL

  class(regulons)<-"regulon"
  return(regulons)
}

#' Copied from VIPER package
#'
updateRegulon <- function(regul) {
  if (is.null(names(regul[[1]]))) {
    tmp <- lapply(regul, function(x) {
      tmp <- rep(0, length(x))
      names(tmp) <- x
      list(tfmode=tmp, likelihood=rep(1, length(tmp)))
    })
    return(tmp)
  }
  if (names(regul[[1]])[1]=="tfmode") return(regul)
  return(lapply(regul, function(x) list(tfmode=x, likelihood=rep(1, length(x)))))
}

#' Copied from VIPER package
#'
TFmode1 <- function (regulon, expset, method = "spearman") {
  regulon <- updateRegulon(regulon)
  regulon <- regulon[names(regulon) %in% rownames(expset)]
  regulon <- lapply(regulon, function(x, genes) {
    filtro <- names(x$tfmode) %in% genes
    x$tfmode <- x$tfmode[filtro]
    if (length(x$likelihood) == length(filtro))
      x$likelihood <- x$likelihood[filtro]
    return(x)
  }, genes = rownames(expset))
  tf <- unique(names(regulon))
  tg <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  cmat <- cor(t(expset[rownames(expset) %in% tf, ]), t(expset[rownames(expset) %in% tg, ]), method = method)
  reg <- lapply(1:length(regulon), function(i, regulon, cmat) {
    tfscore <- cmat[which(rownames(cmat) == names(regulon)[i]), match(names(regulon[[i]]$tfmode), colnames(cmat))]
    list(tfmode = tfscore, likelihood = regulon[[i]]$likelihood)
  }, regulon = regulon, cmat = cmat)
  names(reg) <- names(regulon)
  return(reg)
}

#' Copied from VIPER package
#'
#' @importFrom mixtools normalmixEM
TFscore <- function (regul, mu = NULL, sigma = NULL, verbose=TRUE) {
  if (length(mu) == 3 & length(sigma) == 3)
    fit <- list(mu = mu, sigma = sigma)
  else {
    tmp <- unlist(lapply(regul, function(x) x$tfmode), use.names = FALSE)
    fit <- list(mu = c(-0.5, 0, 0.5), sigma = c(0.15, 0.25, 0.15), lambda = c(0.2, 0.4, 0.4), all.loglik = rep(0, 1001))
    count <- 0
    while (length(fit$all.loglik) > 1000 & count < 3) {
      fit <- normalmixEM(tmp, mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, mean.constr = c(NA, 0, NA), maxit = 1000, verb = FALSE)
      count <- count + 1
    }
  }
  if (verbose) message("mu: ", paste(fit$mu, collapse = ", "), ". sigma: ", paste(fit$sigma, collapse = ", "))
  regul <- lapply(regul, function(x, fit) {
    g2 <- pnorm(x$tfmode, fit$mu[3], fit$sigma[3], lower.tail = TRUE)
    g1 <- pnorm(x$tfmode, fit$mu[1], fit$sigma[1], lower.tail = FALSE)
    g0 <- pnorm(x$tfmode, fit$mu[2], fit$sigma[2], lower.tail = FALSE)
    g00 <- pnorm(x$tfmode, fit$mu[2], fit$sigma[2], lower.tail = TRUE)
    x$tfmode <- g2/(g1 + g0 + g2) * (x$tfmode >= 0) - g1/(g1 + g00 + g2) * (x$tfmode < 0)
    return(x)
  }, fit = fit)
  return(regul)
}

#' Import hpARACNe network
#'
#' This function imports a hpARACNe network and generates site-specific regulons
#'
#' @param afile hpARACNe network file
#' @param pfile hpARACNe peptides file
#' @param mfile Optional hpARACNe quantitative matrix file for TFmode refinement
#' @param method Correlation method for TFmode refinement
#' @param confidence_threshold Interaction confidence_threshold
#' @param priors Logical, whether prior interaction probabilities should be used.
#' @param verbose Logical, whether progression messages should be printed in the terminal.
#' @import data.table
#' @importFrom plyr dlply .
#' @export
hparacne2regulon<-function(afile, pfile, mfile=NA, method="spearman", confidence_threshold=0, priors=TRUE, verbose=TRUE) {
  aracne<-fread(afile)
  if (!("Prior" %in% names(aracne))) {
    aracne$Prior<-1
  }
  peptides<-fread(pfile)

  # Map peptide identifiers to site identififers
  aracne<-merge(aracne,peptides[,c("peptide_id","site_id"), with = FALSE], by.x="Regulator", by.y="peptide_id", allow.cartesian=TRUE)
  aracne<-merge(aracne,peptides[,c("peptide_id","site_id"), with = FALSE], by.x="Target", by.y="peptide_id", allow.cartesian=TRUE)
  aracne<-aracne[,c("site_id.x","site_id.y","MI","Correlation","Prior"), with = FALSE]
  names(aracne)<-c("Regulator","Target","MI","Correlation","Prior")

  # Apply confidence threshold
  aracne<-subset(aracne, Prior > confidence_threshold)

  # Parse quantitative matrix if specified
  if (class(mfile)=="matrix") {
    mx<-mfile
    aracne<-subset(aracne, Regulator %in% rownames(mx) & Target %in% rownames(mx))
  }
  else if (!is.na(mfile)) {
    mx<-fread(mfile)
    mx<-merge(mx,peptides[,c("peptide_id","site_id"), with = FALSE], by="peptide_id", allow.cartesian=TRUE)
    mx_ids<-mx$site_id
    mx[,peptide_id:=NULL]
    mx[,site_id:=NULL]
    mx<-as.matrix(mx)
    rownames(mx)<-mx_ids
    mx[which(is.na(mx))]<-0
  }

  # Compute likelihood from MI
  aracne$likelihood<-(aracne$MI / max(aracne$MI))

  # Use priors
  if (priors) {
    aracne$likelihood<-aracne$likelihood * aracne$Prior
  }

  # Generate raw regulon
  regulons<-dlply(aracne,.(Regulator),function(X){tfmode<-X$Correlation;names(tfmode)<-X$Target;return(list("tfmode"=tfmode, "likelihood"=X$likelihood))})

  # Compute new TFmode if quantitative matrix is present
  if (!is.na(mfile)) {
    regulons<-TFmode1(regulons, mx, method)
  }

  # Remove missing data
  regulons <- regulons[names(regulons) != "NA"]
  regulons <- lapply(regulons, function(x) {
    filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood))
    x$tfmode <- x$tfmode[filtro]
    x$likelihood <- x$likelihood[filtro]
    return(x)
  })
  regulons <- regulons[sapply(regulons, function(x) length(names(x$tfmode)))>0]

  # Transform TFmode
  regulons <- TFscore(regulons, verbose=verbose)

  class(regulons) <- "regulon"

  return(regulons)
}

#' Convert site-specific regulators to protein regulators
#'
#' This function converts site-specific regulators to site-specific regulators
#'
#' @param regulons VIPER regulons
#' @import data.table
#' @importFrom plyr dlply .
#' @export
regulator2protein<-function(regulons) {
  # Generate regulon map
  regulon_map<-data.table("site_id"=names(regulons),"protein_id"=as.vector(sapply(names(regulons),function(X){strsplit(X,":")[[1]][2]})))

  # For each protein site, generate a unique regulon
  regulon_map[,regulon_id:=c(1:length(unique(site_id))),key="protein_id"]

  # Generate list of regulons
  regulons_list<-dlply(regulon_map,.(regulon_id),function(X){subregulons<-regulons[X$site_id];names(subregulons)<-X$protein_id;return(subregulons)})

  return(regulons_list)
}

#' Convert site-specific regulons to protein regulons
#'
#' This function converts site-specific regulons to a list of protein - protein regulons
#'
#' @param regulons VIPER regulons
#' @import data.table
#' @importFrom plyr dlply .
#' @export
site2protein<-function(regulons, average_regulons = TRUE) {
  # Generate site_id map
  site_ids<-unique(c(names(regulons),unlist(sapply(regulons,function(X){names(X$tfmode)}))))
  id_map<-data.table("site_id"=site_ids,"protein_id"=as.vector(sapply(site_ids,function(X){strsplit(X,":")[[1]][2]})))

  # Generate regulon map
  regulon_map<-data.table("site_id"=names(regulons),"protein_id"=as.vector(sapply(names(regulons),function(X){strsplit(X,":")[[1]][2]})))

  # For each protein site, generate a unique regulon
  if (average_regulons) {
    regulon_map$regulon_id = 1
  }
  else {
    regulon_map[,regulon_id:=c(1:length(unique(site_id))),key="protein_id"]
  }

  # Generate list of regulons
  regulons_list<-dlply(regulon_map,.(regulon_id),function(X){subregulons<-regulons[X$site_id];names(subregulons)<-X$protein_id;subregulons<-lapply(subregulons,function(Y){ids<-names(Y$tfmode);names(Y$tfmode)<-with(id_map, protein_id[match(ids,site_id)]);return(Y)});return(subregulons)})

  regulons_list_averaged<-lapply(regulons_list,function(X){lapply(X,function(Y){res<-sapply(unique(names(Y$tfmode)),function(Z){idx<-which(names(Y$tfmode)==Z);return(list("tfmode"=mean(Y$tfmode[idx]),"likelihood"=mean(Y$likelihood[idx])))});return(list("tfmode"=unlist(res[1,]),"likelihood"=unname(unlist(res[2,]))))})})

  if (average_regulons) {
    regulons_list_averaged = regulons_list_averaged[[1]]
    class(regulons_list_averaged) <- "regulon"
  }

  return(regulons_list_averaged)
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
  regdf<-ldply(names(regulons),function(X){data.frame("site_id"=X, "target_id"=names(regulons[[X]]$tfmode), "tfmode"=regulons[[X]]$tfmode, "likelihood"=regulons[[X]]$likelihood, stringsAsFactors = FALSE)})
  regdf$regulator_id<-as.vector(sapply(regdf$site_id,function(X){strsplit(X,":")[[1]][1]}))
  regdf$target_id<-as.vector(sapply(regdf$target_id,function(X){strsplit(X,":")[[1]][1]}))

  # Average interactions
  regdfa<-ddply(regdf[,c("regulator_id","target_id","tfmode","likelihood")],.(regulator_id, target_id),function(X){data.frame("tfmode"=mean(X$tfmode), "likelihood"=mean(X$likelihood))})

  # Generate regulons
  regulons_genes<-dlply(regdfa,.(regulator_id),function(X){tfmode<-X$tfmode; names(tfmode)<-X$target_id;return(list("tfmode"=tfmode, "likelihood"=X$likelihood))})
  class(regulons_genes) <- "regulon"

  return(regulons_genes)
}

#' Subset regulons for targets
#'
#' This function subsets regulons to specified targets only
#'
#' @param regulons VIPER regulons
#' @param target Vector of targets
#' @param min_size minimum regulon size
#' @export
subset_regulon<-function(regulons,targets,min_size=10) {
  subregulon<-lapply(regulons,function(X){ids<-which(names(X$tfmode) %in% targets);if(length(ids)>min_size){return(list("tfmode"=X$tfmode[ids], "likelihood"=X$likelihood[ids]))}})

  return(subregulon[sapply(subregulon,length)>0])
}
