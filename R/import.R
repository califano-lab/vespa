replace_sanger_ptms<-function(X){
  X<-str_replace_all(X,"pS","s")
  X<-str_replace_all(X,"pT","t")
  X<-str_replace_all(X,"pY","y")
  return(X)
}

replace_tpp_ptms<-function(X){
  X<-str_replace_all(X,"S\\[79.9663\\]","s")
  X<-str_replace_all(X,"T\\[79.9663\\]","t")
  X<-str_replace_all(X,"Y\\[79.9663\\]","y")
  X<-str_replace_all(X,"n\\[42.0106\\]","")
  X<-str_replace_all(X,"\\.","")
  X<-gsub("\\s*\\[[^\\)]+\\]","",as.character(X))
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

replace_oldmq_ptms<-function(X){
  X<-str_replace_all(X,"S\\(Phospho \\(STY\\)\\)","s")
  X<-str_replace_all(X,"T\\(Phospho \\(STY\\)\\)","t")
  X<-str_replace_all(X,"Y\\(Phospho \\(STY\\)\\)","y")
  X<-str_replace_all(X,"\\(Oxidation \\(M\\)\\)","")
  X<-str_replace_all(X,"\\(Acetyl \\(Protein N-term\\)\\)","")
  X<-str_replace_all(X,"\\_","")
  X<-gsub("\\s*\\([^\\)]+\\)","",as.character(X))
  return(X)
}

replace_mq_ptms<-function(X){
  X<-str_replace_all(X,"pS","s")
  X<-str_replace_all(X,"pT","t")
  X<-str_replace_all(X,"pY","y")
  X<-str_replace_all(X,"S\\(ph\\)","s")
  X<-str_replace_all(X,"T\\(ph\\)","t")
  X<-str_replace_all(X,"Y\\(ph\\)","y")
  X<-str_replace_all(X,"\\.","")
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

  if (dim(protein_sites)[1]>0) {
    protein_sites$phosphosite<-paste(protein_sites$aa,protein_sites$site,sep="")
  }
  else {
    protein_sites<-data.frame("phosphosite"=NA, stringsAsFactors = FALSE)
  }
  protein_sites$site_id<-paste(gene_id,protein_id,protein_sites$phosphosite,sep=":")

  return(protein_sites[,c("site_id","phosphosite")])
}

normalize_pvt<-function(pvt, normalization_method="quantile", normalization_center=TRUE) {
  # generate matrix
  datmx<-dcast(unique(pvt[,c("peptide_id","tag","peptide_intensity")]), peptide_id ~ tag, value.var = "peptide_intensity")

  # get annotation
  rowids<-datmx$peptide_id
  colids<-colnames(datmx[,-1])
  datmx<-as.matrix(datmx[,-1])
  datmx[datmx == 0]<-NA

  # log transform matrix
  datmx<-log2(datmx)
  datmx[datmx<0]<-NA

  # center run-wise by median
  if (normalization_center) {
    datmx<-apply(datmx,2,function(X){return(X-median(X, na.rm=TRUE))})
  }

  # normalize matrix
  if (normalization_method==FALSE) {
    datmxn<-data.table(datmx)
    datmxn$peptide_id<-rowids
  }
  else if (normalization_method=="quantile") {
    datmxn<-normalize.quantiles(datmx)
    colnames(datmxn)<-colids
    datmxn<-data.table(datmxn)
    datmxn$peptide_id<-rowids
  }
  else if (normalization_method=="cyclicLoess") {
    datmxn<-normalizeCyclicLoess(datmx)
    colnames(datmxn)<-colids
    datmxn<-data.table(datmxn)
    datmxn$peptide_id<-rowids
  }
  else {
    stop("Error: Unknown normalization method")
  }

  pvtn<-melt(datmxn, id.vars="peptide_id", variable.name="tag", value.name="peptide_intensity")

  return(pvtn)
}

#' Import OpenSWATH TSV file
#'
#' This function imports a OpenSWATH TSV file and converts the data to the unified phosphoviper format
#'
#' @param file OpenSWATH TSV file
#' @param fasta Amino acid FASTA file from UniProt
#' @param normalization_method Either "FALSE" (skip normalization), "quantile" or "cyclicLoess" normalization
#' @param normalization_center Either FALSE (skip center normalization) or TRUE
#' @param batchfile A data.table with tags (replacement for run_id) and batch annotation for separate normalization (columns: run_id, tag, aggregator_id, ds_id)
#' @param cores Integer indicating the number of cores to use
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @importFrom limma normalizeCyclicLoess
#' @export
importOpenSWATH<-function(file, fasta, normalization_method = FALSE, normalization_center=TRUE, batchfile, cores = 1) {

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

  # modify run identifier
  dat$run_id<-as.vector(sapply(dat$run_id,function(X){strsplit(X,"\\/\\/")[[1]][2]}))
  dat$run_id<-as.vector(sapply(dat$run_id,function(X){strsplit(X,"\\.")[[1]][1]}))

  # append tag identifier
  dat<-merge(dat, batchfile[,c("run_id","tag")], by="run_id", allow.cartesian=TRUE)
  dat$run_id<-dat$tag

  # map peptides to proteins
  dat$phospho_peptide_sequence<-strip_ptms(replace_ptms(dat$modified_peptide_sequence))
  peptides_proteins<-data.table("peptide_sequence"=unique(dat$peptide_sequence))
  message("Mapping OpenSWATH peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datl<-merge(merge(peptides_proteins, dat, by="peptide_sequence"), genes_proteins, by="protein_id")

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","phospho_peptide_sequence")]),.(phospho_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$phospho_peptide_sequence)}))
  site_mapping<-merge(site_mapping, unique(datl[,c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence","phospho_peptide_sequence")]),by="phospho_peptide_sequence")

  # normalize data
  if ("aggregator_id" %in% colnames(batchfile)) {
    datln<-ddply(batchfile,.(aggregator_id),function(X){return(normalize_pvt(subset(datl, tag %in% X$tag), normalization_method, normalization_center))})
  }
  else {
    datln<-normalize_pvt(datl, normalization_method)
  }

  datl<-merge(datln, site_mapping,by="peptide_id", allow.cartesian=TRUE)

  datl<-subset(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","tag","peptide_intensity")],!is.na(peptide_intensity))
  names(datl)<-c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")

  return(datl)
}

#' Import IonQuant file
#'
#' This function imports an IonQuant mbr_ion.tsv or MSstats.csv file and converts the data to the unified phosphoviper format
#'
#' @param file IonQuant mbr_ion.tsv or MSstats.csv file
#' @param fasta Amino acid FASTA file from UniProt
#' @param normalization_method Either "FALSE" (skip normalization), "quantile" or "cyclicLoess" normalization
#' @param normalization_center Either FALSE (skip center normalization) or TRUE
#' @param batchfile A data.table with tags (replacement for run_id) and batch annotation for separate normalization (columns: run_id, tag, aggregator_id, ds_id)
#' @param cores Integer indicating the number of cores to use
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @importFrom limma normalizeCyclicLoess
#' @export
importIonQuant<-function(file, fasta, normalization_method = FALSE, normalization_center=TRUE, batchfile, cores = 1) {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load IonQuant data
  dat<-fread(file)
  if ("PeptideSequence" %in% colnames(dat)) {
    message("Loading IonQuant MSstats.csv data")
    dat<-dat[,c("PeptideSequence","PrecursorCharge","Run","Intensity")]
    names(dat)<-c("modified_peptide_sequence","precursor_charge","run_id","peptide_intensity")
    dat$pep<-0
  } else if ("modified_peptide" %in% colnames(dat)) {
    message("Loading IonQuant mbr_ion.tsv data")
    dat<-fread(file)[,c("modified_peptide","charge","acceptor_run_name","probability","intensity")]
    names(dat)<-c("modified_peptide_sequence","precursor_charge","run_id","pep","peptide_intensity")
    dat$pep<-(1-dat$pep)
  } else {
    stop("Error: IonQuant file format unknown")
  }

  dat<-subset(dat, !is.na(peptide_intensity))
  dat$peptide_id<-paste(dat$modified_peptide_sequence,dat$precursor_charge,sep="_")

  # select best scoring peptide per run if multiple are available
  message("Selecting peptides")
  dat<-ddply(dat,.(peptide_id, run_id),function(X){return(X[which(X$pep==min(X$pep))[1],])})

  # append tag identifier
  dat<-merge(dat, batchfile[,c("run_id","tag")], by="run_id", allow.cartesian=TRUE)
  dat$run_id<-dat$tag

  # map peptides to proteins
  dat$phospho_peptide_sequence<-replace_tpp_ptms(dat$modified_peptide_sequence)
  dat$peptide_sequence<-str_to_upper(dat$phospho_peptide_sequence)
  # dat<-subset(dat, peptide_sequence != modified_peptide_sequence)
  peptides_proteins<-data.table("peptide_sequence"=unique(dat$peptide_sequence))
  message("Mapping IonQuant peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datl<-merge(merge(peptides_proteins, dat, by="peptide_sequence"), genes_proteins, by="protein_id")

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","phospho_peptide_sequence")]),.(phospho_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$phospho_peptide_sequence)}))
  site_mapping<-merge(site_mapping, unique(datl[,c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence","phospho_peptide_sequence")]),by="phospho_peptide_sequence")

  # normalize data
  if ("aggregator_id" %in% colnames(batchfile)) {
    datln<-ddply(batchfile,.(aggregator_id),function(X){return(normalize_pvt(subset(datl, tag %in% X$tag), normalization_method, normalization_center))})
  }
  else {
    datln<-normalize_pvt(datl, normalization_method, normalization_center)
  }

  datl<-merge(merge(datln, unique(datl[,c("peptide_id","tag","pep")]), by=c("peptide_id","tag")), site_mapping,by="peptide_id", allow.cartesian=TRUE)

  datl<-subset(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","tag","peptide_intensity","pep")],!is.na(peptide_intensity))
  names(datl)<-c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity","pep")

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity","pep")])
}

#' Import protein MaxQuant proteinGroups TXT file
#'
#' This function imports a protein MaxQuant proteinGroups TXT file and converts the data to the unified phosphoviper format
#'
#' @param file protein MaxQuant proteinGroups TXT file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @export
importProteoMaxQuant<-function(file, fasta, cores = 1) {

  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))

  # load MaxQuant data
  message("Loading MaxQuant data")
  dat<-fread(file)

  # parse run identifiers
  run_ids<-colnames(dat)[str_detect(colnames(dat),"Intensity ")]

  dat<-dat[,c("Protein IDs","Number of proteins", run_ids),with = FALSE]
  colnames(dat)<-c("protein_id","n_proteins",str_replace(run_ids,"Intensity ",""))

  # proteotypic peptides only
  dat<-subset(dat, n_proteins==1)[,c("protein_id", str_replace(run_ids,"Intensity ","")),with = FALSE]

  datl<-melt(dat, id.vars="protein_id", variable.name="run_id", value.name="peptide_intensity")
  datl<-subset(datl, !is.na(peptide_intensity) & peptide_intensity > 0)

  datl<-merge(datl, genes_proteins,by="protein_id", allow.cartesian=TRUE)

  datl$peptide_id<-datl$protein_id
  datl$modified_peptide_sequence<-datl$protein_id
  datl$peptide_sequence<-datl$protein_id
  datl$phosphosite<-"PA"
  datl$site_id<-paste(datl$gene_id,datl$protein_id,datl$phosphosite,sep=":")
  datl$peptide_intensity<-log10(datl$peptide_intensity)

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import phospho MaxQuant Evidence TXT file
#'
#' This function imports a phospho MaxQuant Evidence TXT file and converts the data to the unified phosphoviper format
#'
#' @param file phospho MaxQuant Evidence TXT file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @import preprocessCore
#' @export
importPhosphoMaxQuant<-function(file, fasta, cores = 1) {

  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))
  names(fasta)<-as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]}))

  # load MaxQuant data
  message("Loading MaxQuant data")
  dat<-fread(file)

  legacy_mode<-FALSE
  if (!("Experiment" %in% colnames(dat))) {
    colids<-c("Modified sequence","Charge","Raw file","Intensity","PEP")
    legacy_mode<-TRUE
  } else if ("Modified sequence" %in% colnames(dat)) {
    colids<-c("Modified sequence","Charge","Experiment","Intensity","PEP")
  } else {
    colids<-c("Modified Sequence","Charge","Experiment","Intensity","PEP")
  }

  dat<-dat[,colids,with = FALSE]

  names(dat)<-c("peptide_id","precursor_charge","run_id","peptide_intensity","pep")

  dat$modified_peptide_sequence<-sapply(dat$peptide_id,function(X){str_replace_all(X,"\\_","")})
  dat$peptide_id<-paste(dat$peptide_id,dat$precursor_charge,sep=".")

  # map peptides to proteins
  if (legacy_mode) {
    dat$modified_peptide_sequence<-replace_oldmq_ptms(dat$modified_peptide_sequence)
  } else {
    dat$modified_peptide_sequence<-replace_mq_ptms(dat$modified_peptide_sequence)
  }
  dat$peptide_sequence<-str_to_upper(strip_ptms(dat$modified_peptide_sequence))
  dat<-subset(dat, peptide_sequence != modified_peptide_sequence)
  peptides_proteins<-data.table("peptide_sequence"=unique(dat$peptide_sequence))
  message("Mapping MaxQuant peptides to FASTA")
  peptides_proteins$protein_id<-names(fasta)[pbsapply(peptides_proteins$peptide_sequence,function(X){index<-str_which(fasta,X);if(length(index)==1){return(index)}else{return(length(fasta)+1)}}, cl=cores)]
  datl<-merge(merge(peptides_proteins, dat, by="peptide_sequence"), genes_proteins, by="protein_id")

  # map phosphosites
  message("Mapping phosphopeptide sites")
  site_mapping<-as.data.table(ddply(unique(datl[,c("gene_id","protein_id","modified_peptide_sequence")]),.(modified_peptide_sequence),function(X){convert_sites(fasta,X$gene_id,X$protein_id,X$modified_peptide_sequence)}))
  site_mapping<-merge(site_mapping, unique(datl[,c("gene_id","protein_id","peptide_id","peptide_sequence","modified_peptide_sequence")]),by="modified_peptide_sequence")

  # reduce list to best results
  datl<-datl[which(is.finite(datl$peptide_intensity)),]
  datl[which(!is.finite(datl$pep)),"pep"]<-1
  datl<-ddply(datl[which(is.finite(datl$peptide_intensity)),],.(peptide_id, run_id), function(X){return(data.frame("peptide_intensity"=X[which(X$pep==min(X$pep))[1],"peptide_intensity"]))})

  datl$peptide_intensity<-log10(datl$peptide_intensity)

  # generate matrix
  datmx<-dcast(data.table(datl[,c("peptide_id","run_id","peptide_intensity")]), peptide_id ~ run_id, value.var = "peptide_intensity")

  # get annotation
  rowids<-datmx$peptide_id
  colids<-colnames(datmx[,-1])
  datmx<-as.matrix(datmx[,-1])
  datmx[datmx == 0]<-NA

  # normalize matrix
  datmxn<-normalize.quantiles(datmx)
  colnames(datmxn)<-colids
  datmxn<-data.table(datmxn)
  datmxn$peptide_id<-rowids

  datln<-melt(datmxn, id.vars="peptide_id", variable.name="run_id", value.name="peptide_intensity")
  datln<-subset(datln, !is.na(peptide_intensity))

  datl<-merge(datln, site_mapping,by="peptide_id", allow.cartesian=TRUE)

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import protein CPTAC file
#'
#' This function imports a protein CPTAC file and converts the data to the unified phosphoviper format
#'
#' @param file protein-level CPTAC file
#' @param fasta Amino acid FASTA file from UniProt
#' @param hgnc HGNC mapping file
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @import stringr
#' @importFrom plyr ddply .
#' @importFrom tidyr separate_rows
#' @export
importProteoCPTAC<-function(file, fasta, hgnc) {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))

  # load HGNC data
  message("Loading HGNC data")
  hgnc_proteins<-fread(hgnc)[,c("hgnc_id","uniprot_id")]
  hgnc_proteins<-separate_rows(hgnc_proteins, uniprot_id)
  names(hgnc_proteins)<-c("hgnc_id","protein_id")
  hgnc_genes_proteins<-merge(genes_proteins, hgnc_proteins, by="protein_id")

  # load CPTAC data
  message("Loading CPTAC data")
  dat<-fread(file)
  dat<-dat[,c("Authority",names(dat)[str_detect(names(dat), "Unshared Log Ratio")]),with=FALSE]

  dat$hgnc_id<-dat$Authority
  dat$gene_id<-hgnc_genes_proteins[match(dat$hgnc_id,hgnc_genes_proteins$hgnc_id),]$gene_id
  dat$protein_id<-hgnc_genes_proteins[match(dat$hgnc_id,hgnc_genes_proteins$hgnc_id),]$protein_id
  dat$phosphosite<-"PA"
  dat$peptide_id<-dat$protein_id
  dat$peptide_sequence<-dat$protein_id
  dat$modified_peptide_sequence<-dat$protein_id
  dat[,hgnc_id:=NULL]
  dat[,Authority:=NULL]
  dat$site_id<-paste(dat$gene_id,dat$protein_id,dat$phosphosite,sep=":")

  # transform data to list
  datl<-melt(dat, id.vars=c("gene_id","protein_id","site_id","peptide_id","peptide_sequence","modified_peptide_sequence","phosphosite"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl$protein_id),]
  datl<-subset(datl, protein_id != "")

  # modify run identifier
  datl$run_id<-sapply(as.character(datl$run_id),function(X){strsplit(X," ")[[1]][1]})

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import phospho CPTAC file
#'
#' This function imports a phospho CPTAC file and converts the data to the unified phosphoviper format
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
importPhosphoCPTAC<-function(file, fasta, cores = 1) {
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
  peptides_phosphopeptides$peptide_id<-paste("PEP",row.names(peptides_phosphopeptides),sep="_")
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

#' Import protein CCT file
#'
#' This function imports a protein CCT file and converts the data to the unified phosphoviper format
#'
#' @param file protein-level CCT file
#' @param fasta Amino acid FASTA file from UniProt
#' @param hgnc HGNC mapping file
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @importFrom tidyr separate_rows
#' @export
importProteoCCT<-function(file, fasta, hgnc) {
  # load reference fasta library
  message("Loading FASTA DB")
  fasta<-read.fasta(fasta, seqtype="AA", as.string = TRUE, set.attributes = FALSE)
  genes_proteins<-data.table("gene_id"=as.vector(sapply(names(fasta),function(X){strsplit(strsplit(X,"\\|")[[1]][3],"_")[[1]][1]})), "protein_id"=as.vector(sapply(names(fasta),function(X){strsplit(X,"\\|")[[1]][2]})))

  # load HGNC data
  message("Loading HGNC data")
  hgnc_proteins<-fread(hgnc)[,c("symbol_id","uniprot_id")]
  hgnc_proteins<-separate_rows(hgnc_proteins, uniprot_id)
  names(hgnc_proteins)<-c("symbol_id","protein_id")
  hgnc_genes_proteins<-merge(genes_proteins, hgnc_proteins, by="protein_id")

  # load CCT data
  message("Loading CCT data")
  dat<-fread(file)
  dat$symbol_id<-dat$attrib_name
  dat$gene_id<-hgnc_genes_proteins[match(dat$symbol_id,hgnc_genes_proteins$symbol_id),]$gene_id
  dat$protein_id<-hgnc_genes_proteins[match(dat$symbol_id,hgnc_genes_proteins$symbol_id),]$protein_id
  dat$phosphosite<-"PA"
  dat$peptide_id<-dat$protein_id
  dat$peptide_sequence<-dat$protein_id
  dat$modified_peptide_sequence<-dat$protein_id
  dat[,attrib_name:=NULL]
  dat[,symbol_id:=NULL]
  dat$site_id<-paste(dat$gene_id,dat$protein_id,dat$phosphosite,sep=":")

  # transform data to list
  datl<-melt(dat, id.vars=c("gene_id","protein_id","site_id","peptide_id","peptide_sequence","modified_peptide_sequence","phosphosite"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl$protein_id),]
  datl<-subset(datl, protein_id != "")

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import phospho CCT file
#'
#' This function imports a phospho CCT file and converts the data to the unified phosphoviper format
#'
#' @param file CCT file
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @export
importPhosphoCCT<-function(file) {
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

#' Import protein SANGER file
#'
#' This function imports a protein SANGER file and converts the data to the unified phosphoviper format
#'
#' @param file protein SANGER file
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @export
importProteoSANGER<-function(file) {
  # load SANGER data
  message("Loading SANGER data")
  dat<-fread(file)
  dat$gene_id<-dat$`Gene name`
  dat$protein_id<-dat$Accession
  dat$phosphosite<-"PA"
  dat$peptide_id<-dat$protein_id
  dat$peptide_sequence<-dat$protein_id
  dat$modified_peptide_sequence<-dat$protein_id
  dat[,Accession:=NULL]
  dat[,`Gene name`:=NULL]
  dat$site_id<-paste(dat$gene_id,dat$protein_id,dat$phosphosite,sep=":")

  # transform data to list
  datl<-melt(dat, id.vars=c("gene_id","protein_id","site_id","peptide_id","peptide_sequence","modified_peptide_sequence","phosphosite"), variable.name="run_id", value.name="peptide_intensity")
  datl<-datl[complete.cases(datl$protein_id),]
  datl<-subset(datl, protein_id != "")

  return(datl[,c("gene_id","protein_id","peptide_id","site_id","modified_peptide_sequence","peptide_sequence","phosphosite","run_id","peptide_intensity")])
}

#' Import phospho SANGER file
#'
#' This function imports a phospho SANGER file and converts the data to the unified phosphoviper format
#'
#' @param file phospho SANGER file
#' @param fasta Amino acid FASTA file from UniProt
#' @param cores Integer indicating the number of cores to use (only 1 in Windows-based systems)
#' @return peptide-level phosphoviper data.table
#' @import data.table
#' @importFrom plyr ddply .
#' @import seqinr
#' @import stringr
#' @import pbapply
#' @export
importPhosphoSANGER<-function(file, fasta, cores = 1) {
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
