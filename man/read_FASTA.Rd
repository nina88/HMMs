\name{read_FASTA}
\alias{combine_FASTA}
\alias{convert}
\alias{read_FASTA}
\alias{read_FASTA_string}
\alias{reduce_FASTA}
\title{Makes hmm_fasta object}
\usage{
  read_FASTA(filename, level)

  read_FASTA_string(fasta_seq, level)

  combine_FASTA(protein1, protein2)
}
\arguments{
  \item{filename}{location of file}

  \item{level}{value of f}

  \item{fasta_seq}{a vector of characters containing the
  FASTA sequence}

  \item{protein1,}{protein2 an hmm_fasta object generated
  from the read_FASTA function}
}
\value{
  \item{object}{hmm_fasta object}
}
\description{
  Makes hmm_fasta object
}
\author{
  Nina Wilkinson
}
\keyword{character}

