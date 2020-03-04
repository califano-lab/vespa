replace_sanger_ptms<-function(X){
  X<-str_replace_all(X,"pS","s")
  X<-str_replace_all(X,"pT","t")
  X<-str_replace_all(X,"pY","y")
  return(X)
}

replace_spectronaut_ptms<-function(X){
  X<-str_replace_all(X,"S\\[Phospho \\(STY\\)\\]","s")
  X<-str_replace_all(X,"T\\[Phospho \\(STY\\)\\]","t")
  X<-str_replace_all(X,"Y\\[Phospho \\(STY\\)\\]","y")
  X<-str_replace_all(X,"\\_","")
  X<-gsub("\\s*\\([^\\)]+\\)","",as.character(X))
  return(X)
}

replace_ptms<-function(X){
  X<-str_replace_all(X,"S\\(Phospho\\)","s")
  X<-str_replace_all(X,"T\\(Phospho\\)","t")
  X<-str_replace_all(X,"Y\\(Phospho\\)","y")
  X<-str_replace_all(X,"\\.","")
  X<-gsub("\\s*\\([^\\)]+\\)","",as.character(X))
  return(X)
}

strip_ptms<-function(X){
  X<-str_replace_all(X,"\\.","")
  X<-gsub("\\s*\\([^\\)]+\\)","",as.character(X))
  return(X)
}

convert_sites<-function(fasta, gene_id, protein_id, modified_peptide_sequence){
  peptide_start<-str_locate(fasta[protein_id], str_to_upper(modified_peptide_sequence))[1]
  peptide_sites<-str_locate_all(modified_peptide_sequence, c("s","t","y"))

  s_sites<-peptide_sites[[1]][,1] + (peptide_start - 1)
  t_sites<-peptide_sites[[2]][,1] + (peptide_start - 1)
  y_sites<-peptide_sites[[3]][,1] + (peptide_start - 1)

  protein_sites<-rbind(
    data.frame("aa"=rep("S", length(s_sites)), "site"=s_sites, stringsAsFactors = FALSE),
    data.frame("aa"=rep("T", length(t_sites)), "site"=t_sites, stringsAsFactors = FALSE),
    data.frame("aa"=rep("Y", length(y_sites)), "site"=y_sites, stringsAsFactors = FALSE)
    , make.row.names = FALSE)

  protein_sites$phosphosite<-paste(protein_sites$aa,protein_sites$site,sep="")
  protein_sites$site_id<-paste(gene_id,protein_id,protein_sites$phosphosite,sep=":")

  return(protein_sites[,c("site_id","phosphosite")])
}

subsample_pvt<-function(pvt, level=c("protein_id","peptide_id"), batches) {
  subsample_peptides<-function(pvt, level=c("protein_id","peptide_id")) {
    # Filter to overlapping peptides or proteins
    pvt[,ds:=length(unique(ds_id)), by=level]
    pvt<-subset(pvt, ds==length(unique(pvt$ds_id)))

    # Add group variable
    pvt$group_id<-paste(pvt$run_id,pvt$protein_id,pvt$peptide_id,sep="_")

    # Compute peptide or protein observability
    pvt[,obs:=length(unique(group_id)), by=c("ds_id", level)]
    pvt[,total:=length(unique(group_id)), by="ds_id"]
    pvt$freq<-pvt$obs/pvt$total
    pvt[,min_freq:=min(freq), by=level]
    pvt[,samples:=floor(min_freq*total)]
    pvt[,subsample_id:=sample(length(group_id)), by=c("ds_id", level)]

    pvt<-subset(pvt, subsample_id <= samples)
    pvt[,c("ds","group_id","obs","total","freq","min_freq","samples","subsample_id") := NULL]

    return(pvt)
  }

  # merge pvt with batches
  pvt_batches<-merge(pvt, batches[,c("run_id","aggregator_id","ds_id")], by="run_id")

  # subsample pvt batches over aggregator, e.g. separate sample
  pvt_sub<-ddply(pvt_batches,.(aggregator_id),function(X){subsample_peptides(data.table(X), level)})

  return(pvt_sub)
}

subset_pvt<-function(pvt, level=c("protein_id","peptide_id"), batches) {
  subsample_peptides<-function(pvt, level=c("protein_id","peptide_id")) {
    # Filter to overlapping peptides or proteins
    pvt[,ds:=length(unique(ds_id)), by=level]
    pvt<-subset(pvt, ds==length(unique(pvt$ds_id)))

    pvt[,c("ds") := NULL]

    return(pvt)
  }

  # merge pvt with batches
  pvt_batches<-merge(pvt, batches[,c("run_id","aggregator_id","ds_id")], by="run_id")

  # subsample pvt batches over aggregator, e.g. separate sample
  pvt_sub<-ddply(pvt_batches,.(aggregator_id),function(X){subsample_peptides(data.table(X), level)})

  return(pvt_sub)
}

#' Import OpenSWATH TSV file
#'
#' This function imports a OpenSWATH TSV file and converts the data to the unified phosphoviper format
#'
#' @param file OpenSWATH TSV file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @param normalization Either "FALSE" (skip normalization), "quantile" or "cyclicLoess" normalization
#' @param batchcorrection Whether batchcorrection should be skipped ("FALSE") or applied by subsetting ("subset") or subsampling ("subsample").
#' @param batchlevel On which level batch correction should be applied
#' @param batchfile A data.table with batch correction annotation (columns: run_id, aggregator_id, ds_id)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @importFrom limma normalizeCyclicLoess
#' @export
importOpenSWATH<-function(file, fasta, cores = 1, normalization = FALSE, batchcorrection = FALSE, batchlevel = c("protein_id", "peptide_id"), batchfile) {

  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load OpenSWATH data
  message("Loading OpenSWATH data")
  dat<-fread(file)
  dat<-dat[,c("filename","Sequence","FullPeptideName","transition_group_id","Intensity"),with = FALSE]
  names(dat)<-c("run_id","peptide_sequence","modified_peptide_sequence","peptide_id","peptide_intensity")

  # map peptides to proteins
  dat$modified_peptide_sequence<-strip_ptms(replace_ptms(dat$modified_peptide_sequence))
  peptides_proteins<-data.table("peptide_sequence"=unique(dat$peptide_sequence))
  message("Mapping OpenSWATH peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datl<-merge(merge(peptides_proteins, dat, by="peptide_sequence"), genes_proteins, by="protein_id")

  # modify run identifier
  datl$run_id<-as.vector(sapply(datl$run_id,function(X){strsplit(X,"\\/\\/")[[1]][2]}))
  datl$run_id<-as.vector(sapply(datl$run_id,function(X){strsplit(X,"\\.")[[1]][1]}))

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","modified_peptide_sequence")]),.(modified_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$modified_peptide_sequence)}))
  site_mapping<-merge(site_mapping, unique(datl[,c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence")]),by="modified_peptide_sequence")

  # conduct batch correction
  if (batchcorrection=="subset") {
    datl<-subset_pvt(datl, batchlevel, batchfile)
  } else if (batchcorrection=="subsample") {
    datl<-subsample_pvt(datl, batchlevel, batchfile)
  }

  # generate matrix
  datmx<-dcast(unique(datl[,c("peptide_id","run_id","peptide_intensity")]), peptide_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  rowids<-datmx$peptide_id
  colids<-colnames(datmx[,-1])
  datmx<-as.matrix(datmx[,-1])
  datmx[datmx == 0]<-NA

  # log transform matrix
  datmx<-log10(datmx)
  datmx[datmx<0]<-NA

  # normalize matrix
  if (normalization==FALSE) {
    datmxn<-datmx
  }
  else if (normalization=="quantile") {
    datmxn<-normalize.quantiles(datmx)
    colnames(datmxn)<-colids
    datmxn<-data.table(datmxn)
    datmxn$peptide_id<-rowids
  }
  else if (normalization=="cyclicLoess") {
    datmxn<-normalizeCyclicLoess(datmx)
    colnames(datmxn)<-colids
    datmxn<-data.table(datmxn)
    datmxn$peptide_id<-rowids
  }
  else {
    stop("Error: Unknown normalization method")
  }

  datln<-melt(datmxn, id.vars="peptide_id", variable.name="run_id", value.name="peptide_intensity")

  datl<-merge(datln, site_mapping,by="peptide_id", allow.cartesian=TRUE)

  return(subset(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")],!is.na(peptide_intensity)))
}

#' Import Spectronaut TXT file
#'
#' This function imports a Spectronaut TXT file and converts the data to the unified phosphoviper format
#'
#' @param file Spectronaut TXT file
#' @param fasta Amino acid FASTA file from UniProt
#' @param run_ids Column names of runs to extract
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @export
importSpectronaut<-function(file, fasta, run_ids, cores = 1) {

  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load Spectronaut data
  message("Loading Spectronaut data")
  dat<-fread(file)
  dat<-dat[,c("EG.PrecursorId_Phos", "ID_Phos", run_ids),with = FALSE]
  dat<-melt(dat, id.vars=c("EG.PrecursorId_Phos","ID_Phos"), measure.vars=run_ids, variable.name="run_id", value.name="peptide_intensity")

  names(dat)<-c("peptide_id","modified_peptide_sequence","run_id","peptide_intensity")

  # map peptides to proteins
  dat$modified_peptide_sequence<-replace_spectronaut_ptms(dat$modified_peptide_sequence)
  dat$peptide_sequence<-str_to_upper(strip_ptms(dat$modified_peptide_sequence))
  peptides_proteins<-data.table("peptide_sequence"=unique(dat$peptide_sequence))
  message("Mapping Spectronaut peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datl<-merge(merge(peptides_proteins, dat, by="peptide_sequence"), genes_proteins, by="protein_id")

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","modified_peptide_sequence")]),.(modified_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$modified_peptide_sequence)}))
  site_mapping<-merge(site_mapping, unique(datl[,c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence")]),by="modified_peptide_sequence")

  # generate matrix
  datmx<-dcast(unique(datl[,c("peptide_id","run_id","peptide_intensity")]), peptide_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  rowids<-datmx$peptide_id
  colids<-colnames(datmx[,-1])
  datmx<-as.matrix(datmx[,-1])
  datmx[datmx == 0]<-NA

  # normalize matrix
  datmxn<-log10(normalize.quantiles(datmx))
  colnames(datmxn)<-colids
  datmxn<-data.table(datmxn)
  datmxn$peptide_id<-rowids

  datln<-melt(datmxn, id.vars="peptide_id", variable.name="run_id", value.name="peptide_intensity")

  datl<-merge(datln, site_mapping,by="peptide_id", allow.cartesian=TRUE)

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import CPTAC file
#'
#' This function imports a CPTAC file and converts the data to the unified phosphoviper format
#'
#' @param file CPTAC file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @export
importCPTAC<-function(file, fasta, cores = 1) {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load CPTAC data
  message("Loading CPTAC data")
  dat<-fread(file)
  dat<-dat[,c("Phosphopeptide", str_subset(names(dat),"Log Ratio")), with=FALSE]

  # map peptides to proteins
  peptides_phosphopeptides<-data.table("peptide_sequence"=str_to_upper(dat$Phosphopeptide), "modified_peptide_sequence"=dat$Phosphopeptide)
  peptides_phosphopeptides$peptide_id<-row.names(peptides_phosphopeptides)
  peptides_proteins<-data.table("peptide_sequence"=unique(peptides_phosphopeptides$peptide_sequence))
  message("Mapping CPTAC peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datan<-merge(merge(merge(peptides_proteins, peptides_phosphopeptides, by="peptide_sequence"),dat, by.x="modified_peptide_sequence", by.y="Phosphopeptide"), genes_proteins, by="protein_id")

  # transform data to list
  datl<-melt(datan, id.vars=c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl),]

  # modify run identifier
  datl$run_id<-sapply(as.character(datl$run_id),function(X){strsplit(X," ")[[1]][1]})

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","modified_peptide_sequence")]),.(modified_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$modified_peptide_sequence)}))

  datl<-merge(datl,site_mapping[,c("modified_peptide_sequence","phosphosite","site_id")],by="modified_peptide_sequence", allow.cartesian=TRUE)

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import CCT file
#'
#' This function imports a CCT file and converts the data to the unified phosphoviper format
#'
#' @param file CPTAC file
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @export
importCCT<-function(file) {
  # load CCT data
  message("Loading CCT data")
  dat<-fread(file)
  dat$gene_id<-sapply(dat$V1,function(X){strsplit(X,"_")[[1]][1]})
  dat$protein_id<-sapply(sapply(dat$V1,function(X){strsplit(X,"__")[[1]][2]}),function(X){strsplit(X,":")[[1]][1]})
  dat$phosphosite<-sapply(sapply(dat$V1,function(X){strsplit(X,"__")[[1]][2]}),function(X){strsplit(X,":")[[1]][2]})
  dat$peptide_id<-dat$V1
  dat$peptide_sequence<-dat$V1
  dat$modified_peptide_sequence<-dat$V1
  dat[,V1:=NULL]
  dat$site_id<-paste(dat$gene_id,dat$protein_id,dat$phosphosite,sep=":")

  # transform data to list
  datl<-melt(dat, id.vars=c("gene_id","protein_id","site_id","peptide_id","peptide_sequence","modified_peptide_sequence","phosphosite"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl),]

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import SANGER file
#'
#' This function imports a SANGER file and converts the data to the unified phosphoviper format
#'
#' @param file SANGER file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @export
importSANGER<-function(file, fasta, cores = 1) {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load SANGER data
  message("Loading SANGER data")
  dat<-fread(file)
  peptide_ids<-dat[,"Annotated sequence"]
  names(peptide_ids)<-"V1"
  dat$peptide_id<-replace_sanger_ptms(peptide_ids$V1)
  dat[,c("Protein accession","Gene name","Annotated sequence","Protein site","Regulatroy kinases","KEGG name"):=NULL]

  # map peptides to proteins
  peptides_phosphopeptides<-data.table("peptide_id"=dat$peptide_id, "peptide_sequence"=str_to_upper(dat$peptide_id), "modified_peptide_sequence"=dat$peptide_id)
  peptides_proteins<-data.table("peptide_sequence"=unique(peptides_phosphopeptides$peptide_sequence))
  message("Mapping SANGER peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datan<-merge(merge(merge(peptides_proteins, peptides_phosphopeptides, by="peptide_sequence"),dat, by="peptide_id"), genes_proteins, by="protein_id")

  # transform data to list
  datl<-melt(datan, id.vars=c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl),]

  # modify run identifier
  datl$run_id<-sapply(as.character(datl$run_id),function(X){strsplit(X," ")[[1]][1]})

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","modified_peptide_sequence")]),.(modified_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$modified_peptide_sequence)}))

  datl<-merge(datl,site_mapping[,c("modified_peptide_sequence","phosphosite","site_id")],by="modified_peptide_sequence", allow.cartesian=TRUE)

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}
