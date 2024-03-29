#' Export hpARACNe files
#'
#' This function exports hpARACNe input files from the unified vespa format
#'
#' @param datl Unified vespa format
#' @param output_dir Output directory
#' @param kinases List of kinases (UniProtKB)
#' @param phosphatases List of phosphatases (UniProtKB)
#' @param interactions A named list of HSM/D, HSM/P, LPS, LPA or PCdb interaction tables (UniProtKB)
#' @param target_sites (Optional) list of target sites to restrict dataset (PV site_id)
#' @param confidence_threshold Interaction confidence_threshold
#' @param restrict_interactions Restrict interactions to kinase/phosphatase regulator-target subset
#' @param interaction_level Interaction level ("substrate" or "activity")
#' @import data.table
#' @export
export2hparacne<-function(datl, output_dir, kinases, phosphatases, interactions=NULL, target_sites=NULL, confidence_threshold=0, restrict_interactions=TRUE, interaction_level="substrate") {
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

  # select regulators
  if (interaction_level == "activity") {
    regulators <- subset(datl, phosphosite %in% c("PA","PV"))
  } else {
    regulators <- datl
  }

  # select targets (exclude protein abundance)
  targets <- subset(datl, phosphosite!="PA")

  # write kinases
  write.table(unique(subset(regulators, protein_id %in% kinases)$peptide_id), file=file.path(output_dir,"kinases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # write kinases + phosphatases
  write.table(unique(c(subset(regulators, protein_id %in% kinases)$peptide_id, subset(regulators, protein_id %in% phosphatases)$peptide_id)), file=file.path(output_dir,"kinases_phosphatases.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # write targets (only phosphopeptides or VIPER activity are used for targets, but not protein abundance)
  if (!is.null(target_sites)) {
    target_ids<-unique(subset(targets, site_id %in% target_sites)$peptide_id)
  } else {
    target_ids<-unique(targets$peptide_id)
  }
  write.table(target_ids, file=file.path(output_dir,"targets.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")

  # write reference interactions
  if (!is.null(interactions)) {
    sapply(names(interactions), function(X){
      if ("site_id" %in% names(interactions[[X]])) {
        regulator_ppi<-unique(regulators[,c("protein_id","peptide_id")])
        names(regulator_ppi)<-c("regulator","regulator_peptide_id")

        target_ppi<-unique(targets[,c("site_id","peptide_id")])
        names(target_ppi)<-c("site_id","target_peptide_id")

        phosphointeractions<-merge(merge(subset(interactions[[X]], confidence > confidence_threshold), regulator_ppi, by="regulator", allow.cartesian=TRUE), target_ppi, by="site_id", allow.cartesian=TRUE)[,c("regulator_peptide_id","target_peptide_id","confidence")]
      } else {
        regulator_ppi<-unique(regulators[,c("protein_id","peptide_id")])
        names(regulator_ppi)<-c("regulator","regulator_peptide_id")

        target_ppi<-unique(targets[,c("protein_id","peptide_id")])
        names(target_ppi)<-c("target","target_peptide_id")

        phosphointeractions<-merge(merge(subset(interactions[[X]], confidence > confidence_threshold), regulator_ppi, by="regulator", allow.cartesian=TRUE), target_ppi, by="target", allow.cartesian=TRUE)[,c("regulator_peptide_id","target_peptide_id","confidence")]
      }

      if (restrict_interactions) {
        phosphointeractions<-subset(phosphointeractions, regulator_peptide_id %in% unique(c(subset(regulators, protein_id %in% kinases)$peptide_id, subset(regulators, protein_id %in% phosphatases)$peptide_id)))
      }

      write.table(phosphointeractions, file=file.path(output_dir,paste(X,"phosphointeractions.txt",sep="_")), quote=FALSE, row.names=FALSE, sep="\t")
    })
  }
}

#' Export to matrix
#'
#' This function exports tables in the unified vespa format to a matrix
#'
#' @param datl Unified vespa format
#' @param fillvalues Missing value imputation method ("rowmin": Row-wise (peptide) minimum; "colmin": Column-wise (run) minimum; NULL: Skip and use NA)
#' @param jitter Add numerical noise to imputed missing values
#' @return matrix
#' @import data.table
#' @export
export2mx<-function(datl, fillvalues=NA, jitter=TRUE) {
  # reduce data to top peptide query per phosphosite
  pqp_freq<-data.table(merge(unique(datl[,c("site_id","peptide_id")]), as.data.frame(table(datl$peptide_id),stringsAsFactors = FALSE), by.x="peptide_id", by.y="Var1"))
  pqp_top<-pqp_freq[pqp_freq[, .I[which.max(Freq)], by=site_id]$V1]
  datl<-merge(datl, pqp_top[,c("site_id","peptide_id"),which=FALSE], by=c("site_id","peptide_id"))

  # generate matrix
  datm<-data.table::dcast(data.table(unique(datl[,c("site_id","run_id","peptide_intensity")])), site_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  phosphoExp<-as.matrix(datm[,-1])
  phosphoExp[phosphoExp == 0]<-NA

  rownames(phosphoExp)<-datm$site_id

  # fill missing values
  if (is.na(fillvalues) || is.null(fillvalues)) {
    phosphoExp[is.na(phosphoExp)]<-fillvalues
  } else if (fillvalues=="rowmin") {
    # fill missing values row-wise
    if (jitter) {
      phosphoExp<-t(apply(phosphoExp,1,function(X){ja<-diff(sort(unique(X))[1:2]); if (is.na(ja)){ja<-0}; X[is.na(X)]<-jitter(rep(min(X,na.rm=TRUE), sum(is.na(X))), amount=ja); return(X);}))
    } else {
      phosphoExp<-t(apply(phosphoExp,1,function(X){X[is.na(X)]<-min(X,na.rm=TRUE); return(X);}))
    }
  } else if (fillvalues=="colmin") {
    if (jitter) {
      phosphoExp<-apply(phosphoExp,2,function(X){ja<-diff(sort(unique(X))[1:2]); if (is.na(ja)){ja<-0}; X[is.na(X)]<-jitter(rep(min(X,na.rm=TRUE), sum(is.na(X))), amount=ja); return(X);})
    } else {
      phosphoExp<-apply(phosphoExp,2,function(X){X[is.na(X)]<-min(X,na.rm=TRUE); return(X);})
    }
  } else {
    phosphoExp[is.na(phosphoExp)]<-fillvalues
  }

  return(phosphoExp)
}

preprocess_mx<-function(osw, ed) {
  # preprocess peptides
  oswf<-merge(data.table(osw), ed[,c("tag","cl_tag","time_tag","compound_tag")], by.x="run_id", by.y="tag")
  oswf<-unique(oswf[,c("peptide_id","run_id","peptide_intensity","cl_tag","time_tag","compound_tag")])
  oswf<-merge(data.table(oswf), merge(data.table(expand.grid(peptide_id = unique(oswf$peptide_id), run_id = unique(oswf$run_id))), unique(ed[,c("tag","cl_tag","time_tag","compound_tag")]), by.x="run_id", by.y="tag"), by=c("peptide_id","run_id","cl_tag","time_tag","compound_tag"), all=TRUE)

  oswf[, peptide_intensity := nafill(peptide_intensity, type = "const", fill=min(peptide_intensity, na.rm=TRUE)), key="peptide_id"]

  oswm<-merge(data.table(unique(osw[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite")])),data.table(unique(oswf[,c("peptide_id","run_id","peptide_intensity")])), by="peptide_id", allow.cartesian=TRUE)

  return(oswm)
}

#' Export to unified vespa format from viper matrix
#'
#' This function exports a viper matrix to the unified vespa format
#'
#' @param matrix phosphoVIPER quantitative matrix
#' @param fasta Amino acid FASTA file from UniProt
#' @param tag Site tag
#' @return vespa table
#' @import data.table
#' @export
vmx2pv<-function(vmx, fasta=NULL, tag = "PV") {
  # load reference fasta library
  if (!is.null(fasta)) {
    message("Loading FASTA DB")
    fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = TRUE)
    genes_proteins<-data.table("gene_id"=as.vector(sapply(fasta,function(X){strsplit(strsplit(attributes(X)$Annot,"GN=")[[1]][2]," ")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))

    pv<-reshape2::melt(vmx, value.name="peptide_intensity")
    names(pv)<-c("protein_id","run_id","peptide_intensity")
    pv<-merge(genes_proteins, pv, by="protein_id")
    pv$modified_peptide_sequence<-pv$protein_id
    pv$peptide_sequence<-pv$protein_id
    pv$phosphosite<-tag
    pv$site_id<-paste(pv$gene_id,pv$protein_id,pv$phosphosite,sep=":")
    pv$peptide_id<-paste(pv$gene_id,pv$protein_id,pv$phosphosite,sep=":")
  } else {
    if (is.null(tag)) {
      pv<-reshape2::melt(vmx, value.name="peptide_intensity")
      names(pv)<-c("site_id","run_id","peptide_intensity")

      regulon_map<-data.table("site_id"=as.character(unique(pv$site_id)),"gene_id"=as.vector(sapply(as.character(unique(pv$site_id)),function(X){strsplit(X,":")[[1]][1]})),"protein_id"=as.vector(sapply(as.character(unique(pv$site_id)),function(X){strsplit(X,":")[[1]][2]})),"phosphosite"=as.vector(sapply(as.character(unique(pv$site_id)),function(X){strsplit(X,":")[[1]][3]})))

      pv<-merge(pv, regulon_map, by="site_id")
      pv$modified_peptide_sequence<-pv$site_id
      pv$peptide_sequence<-pv$site_id
      pv$peptide_id<-pv$site_id
    } else {
      pv<-reshape2::melt(vmx, value.name="peptide_intensity")
      names(pv)<-c("site_id","run_id","peptide_intensity")

      regulon_map<-data.table("site_id"=as.character(unique(pv$site_id)),"gene_id"=as.vector(sapply(as.character(unique(pv$site_id)),function(X){strsplit(X,":")[[1]][1]})),"protein_id"=as.vector(sapply(as.character(unique(pv$site_id)),function(X){strsplit(X,":")[[1]][2]})))

      pv<-merge(pv, regulon_map, by="site_id")
      pv$site_id<-paste(pv$site_id,tag,sep=":")
      pv$modified_peptide_sequence<-pv$site_id
      pv$peptide_sequence<-pv$site_id
      pv$phosphosite<-tag
      pv$peptide_id<-pv$site_id
    }
  }

  return(pv)
}

#' metaVIPER regulon selection
#'
#' @param ges Gene Expression Signature (features X samples)
#' @param net.list List object with the networks to be used
#' @return Optimized network with best regulon per gene
#' @export
optimizeRegulon <- function(ges, net.list, min_size=10) {
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

    selected_regulon<-net.list[[row.names(w.mat)[which.max(w.mat)]]][[g]]
    if (is.na(as.numeric(names(net.list)[1]))) {
      selected_regulon$meta$origin<-row.names(w.mat)[which.max(w.mat)]
    }

    return(selected_regulon)
  }

  num.nets <- length(net.list)

  ## run VIPER with each network
  print('Generating VIPER matrices...')
  vip.list <- list()
  for (i in 1:num.nets) {
    vip.list[[i]] <- viper(ges, net.list[i], method = 'none')
  }
  names(vip.list) <- names(net.list)

  num.samps <- ncol(vip.list[[1]])

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

#' Orthogonal site-specific regulon selection
#'
#' @param ges Gene Expression Signature (features X samples)
#' @param regulon Site-specific regulons
#' @param min_size Minimum regulon size
#' @param cutoff Pairwise correlation cutoff to remove regulons
#' @return Optimized network containing only orthogonal site-specific regulons for each protein
#' @import plyr
#' @importFrom caret findCorrelation
#' @export
orthogonalRegulon <- function(ges, regulon, min_size=10, cutoff=0.5) {
  # generate VIPER matrix
  subregulon <- subsetRegulon(regulon, min_size = min_size, target = row.names(ges))
  vpmx <- viper(ges, subregulon, method = 'none', minsize=min_size)

  # generate regulon map
  mapregulon<-data.frame("site_id"=names(subregulon),"protein_id"=as.vector(sapply(names(subregulon),function(X){strsplit(X,":")[[1]][2]})))

  # identify correlated regulons
  correlated_regulons<-ddply(mapregulon,.(protein_id),function(X){if (dim(X)[1]>1){data.frame("site_id"=findCorrelation(cor(t(vpmx[X$site_id,])), cutoff = cutoff, names=TRUE))}})

  # report number of regulons to be removed
  message(paste(length(correlated_regulons$site_id),"regulons out of",length(regulon),"are correlated and will be removed."))

  # remove regulons
  regulon[correlated_regulons$site_id]<-NULL

  return(regulon)
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
#' @param proteinlevel Whether protein-level regulator identifiers are used instead of peptide-level.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @import data.table
#' @importFrom plyr dlply .
#' @export
hparacne2regulon<-function(afile, pfile, proteinlevel=FALSE, verbose=TRUE) {
  aracne<-fread(afile)
  if (!("Prior" %in% names(aracne))) {
    aracne$Prior<-NA
  }
  if (!("Correlation" %in% names(aracne))) {
    aracne$Correlation<-NA
  }

  # Load peptides
  peptides<-fread(pfile)

  # Map peptide identifiers to site identifiers
  if (!proteinlevel) {
    regulatorpeptides<-peptides[,c("protein_id","peptide_id","site_id"), with = FALSE]
    names(regulatorpeptides)<-c("regulator_protein_id","regulator_peptide_id","regulator_site_id")
    aracne<-merge(aracne, regulatorpeptides, by.x="Regulator", by.y="regulator_peptide_id", allow.cartesian=TRUE)
  } else {
    aracne$regulator_protein_id<-aracne$regulator_site_id<-aracne$Regulator
  }
  targetpeptides<-peptides[,c("protein_id","peptide_id","site_id"), with = FALSE]
  names(targetpeptides)<-c("target_protein_id","target_peptide_id","target_site_id")
  aracne<-merge(aracne, targetpeptides, by.x="Target", by.y="target_peptide_id", allow.cartesian=TRUE)

  # Remove interactions between sites of the same protein
  aracne<-subset(aracne, regulator_protein_id != target_protein_id)

  aracne<-aracne[,c("regulator_site_id","target_site_id","MI","Correlation","Prior"), with = FALSE]
  names(aracne)<-c("Regulator","Target","MI","Correlation","Prior")

  # Compute likelihood from MI
  aracne$likelihood<-(aracne$MI / max(aracne$MI))

  # Normalize priors by measured maximum PBD-peptide or PPI interaction
  aracne<-data.table(ddply(aracne,.(Regulator),function(X){data.frame("Target"=X$Target, "MI"=X$MI, "Correlation"=X$Correlation, "likelihood"=X$likelihood, "Prior"=X$Prior/max(X$Prior))}))

  # Generate raw regulon
  regulons<-dlply(aracne,.(Regulator),function(X){tfmode<-X$Correlation;names(tfmode)<-X$Target;return(list("tfmode"=tfmode, "likelihood"=X$likelihood, "prior"=X$Prior, "meta"=list("regulator"=unique(X$Regulator), "dataset"=afile)))})

  # Remove missing data
  regulons <- regulons[names(regulons) != "NA"]
  regulons <- lapply(regulons, function(x) {
    filtro <- !(names(x$tfmode)=="NA" | is.na(x$tfmode) | is.na(x$likelihood) | is.na(x$prior))
    x$tfmode <- x$tfmode[filtro]
    x$likelihood <- x$likelihood[filtro]
    x$prior <- x$prior[filtro]
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

  regulons_list_averaged<-lapply(regulons_list,function(X){lapply(X,function(Y){res<-sapply(unique(names(Y$tfmode)),function(Z){idx<-which(names(Y$tfmode)==Z);return(list("tfmode"=mean(Y$tfmode[idx]),"likelihood"=mean(Y$likelihood[idx]),"prior"=mean(Y$prior[idx])))});return(list("tfmode"=unlist(res[1,]),"likelihood"=unname(unlist(res[2,])),"prior"=unname(unlist(res[3,]))))})})

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
  regdf<-ldply(names(regulons),function(X){data.frame("site_id"=X, "target_id"=names(regulons[[X]]$tfmode), "tfmode"=regulons[[X]]$tfmode, "likelihood"=regulons[[X]]$likelihood, "prior"=regulons[[X]]$prior, stringsAsFactors = FALSE)})
  regdf$regulator_id<-as.vector(sapply(regdf$site_id,function(X){strsplit(X,":")[[1]][1]}))
  regdf$target_id<-as.vector(sapply(regdf$target_id,function(X){strsplit(X,":")[[1]][1]}))

  # Average interactions
  regdfa<-ddply(regdf[,c("regulator_id","target_id","tfmode","likelihood","prior")],.(regulator_id, target_id),function(X){data.frame("tfmode"=mean(X$tfmode), "likelihood"=mean(X$likelihood), "prior"=mean(X$prior))})

  # Generate regulons
  regulons_genes<-dlply(regdfa,.(regulator_id),function(X){tfmode<-X$tfmode; names(tfmode)<-X$target_id;return(list("tfmode"=tfmode, "likelihood"=X$likelihood, "prior"=X$prior))})
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
subsetRegulon<-function(regulons,targets,min_size=10) {
  subregulon<-lapply(regulons,function(X){ids<-which(names(X$tfmode) %in% targets);if(length(ids)>min_size){return(list("tfmode"=X$tfmode[ids], "likelihood"=X$likelihood[ids], "prior"=X$prior[ids], "meta"=X$meta))}})

  return(subregulon[sapply(subregulon,length)>0])
}

##########
#' Prune Regulons
#'
#' This function limits the maximum size of the regulons
#'
#' @param regulon Object of class regulon
#' @param cutoff Number indicating the maximum size for the regulons (maximum number of target genes)
#' @param adaptive Logical, whether adaptive size should be used (i.e. sum(likelihood^2))
#' @param eliminate Logical whether regulons smalles than \code{cutoff} should be eliminated
#' @param wm Optional numeric vector of weights (0; 1) for the genes
#' @return Pruned regulon
#' @seealso \code{\link{viper}}, \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' hist(sapply(regulon, function(x) sum(x$likelihood)/max(x$likelihood)), nclass=20)
#' preg <- pruneRegulon(regulon, 400)
#' hist(sapply(preg, function(x) sum(x$likelihood)/max(x$likelihood)), nclass=20)
#' @export
pruneRegulon <- function(regulon, cutoff=50, adaptive=TRUE, eliminate=FALSE) {
  if (adaptive) {
    regulon <- lapply(regulon, function(x, cutoff, wm) {
      likelihood <- x$likelihood
      likelihood <- likelihood * x$prior
      pos <- order(likelihood, decreasing=TRUE)
      ws <- (likelihood/max(likelihood))^2
      pos <- pos[cumsum(ws[pos])<=cutoff]
      return(list(tfmode=x$tfmode[pos], likelihood=x$likelihood[pos], prior=x$prior[pos], meta=x$meta))
    }, cutoff=cutoff, wm=wm)
  }
  else {
    regulon <- lapply(regulon, function(x, cutoff) {
      pos <- order(x$likelihood, decreasing=TRUE)
      pos <- pos[1:min(length(pos), cutoff)]
      return(list(tfmode=x$tfmode[pos], likelihood=x$likelihood[pos], prior=x$prior[pos], meta=x$meta))
    }, cutoff=cutoff)
    if (eliminate) regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=cutoff]
  }
  class(regulon) <- "regulon"
  return(regulon)
}
