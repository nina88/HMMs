#' @name p53
#' @aliases mdm2 CBP TAF15 FUS EWS Bnip3 NBR1 SQSTM1
#' @title Example FASTA sequences
#' @description The amino acid sequence of p53, mdm2, CBP, TAF15, FUS, EWS, 
#' Bnip3, NBR1, SQSTM1. To load a particular sequence, use the data function. 
#' For example, data(p53). For use in an analysis the data must be of class 
#' "hmm_fasta" which is done using the read_FASTA or read_FASTA_string functions. 
#' These functions also convert the amino acid sequences another code which depends 
#' on hydrophobicity (f=2), charge (f=3) or both (f=4).  
#' In the analysis you can combine proteins together to increase the information for the model. 
#' This is done using the combine_FASTA function. 
#' @docType data
#' @format Vector of strings
#' @source Obtained from an online database NCBI http://www.ncbi.nlm.nih.gov/protein
#' @examples 
#' data(TAF15)
#' data(FUS)
#' f=4
#' x1 = read_FASTA_string(TAF15,f)
#' x2 = read_FASTA_string(FUS,f)
#' y = combine_FASTA(x1,x2)
NULL
