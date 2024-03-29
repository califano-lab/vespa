% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{importIonQuant}
\alias{importIonQuant}
\title{Import IonQuant file}
\usage{
importIonQuant(
  file,
  fasta,
  normalization_method = FALSE,
  normalization_center = FALSE,
  batchfile = NULL,
  cores = 1
)
}
\arguments{
\item{file}{IonQuant mbr_ion.tsv or MSstats.csv file}

\item{fasta}{Amino acid FASTA file from UniProt}

\item{normalization_method}{Either "FALSE" (skip normalization), "quantile" or "cyclicLoess" normalization}

\item{normalization_center}{Either FALSE (skip center normalization) or TRUE}

\item{batchfile}{An optional data.table with tags ("tag"; replacement for "run_id") and batch annotation for separate normalization (columns: "run_id", "tag", "aggregator_id", "ds_id")}

\item{cores}{Integer indicating the number of cores to use}
}
\value{
peptide-level vespa data.table
}
\description{
This function imports an IonQuant mbr_ion.tsv or MSstats.csv file and converts the data to the unified vespa format
}
