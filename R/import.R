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
