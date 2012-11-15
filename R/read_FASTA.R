####################################################
#Private functions
####################################################
reduce_FASTA = function(x, from, to) {
    to[match(x,from)]
}

convert = function(fasta_seq, level) {
    if(level==2)
        fasta_seq = reduce_FASTA(fasta_seq,
                                 LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],
                                 to=c(1,1,2,2,1,1,2,1,2,1,1,2,1,2,2,2,2,1,1,1))
    else if(level == 3)
        fasta_seq = reduce_FASTA(fasta_seq,
                                 LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],
                                 to=c(1,1,3,3,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1))
    else if (level ==4)
        fasta_seq = reduce_FASTA(fasta_seq,
                                 LETTERS[c(1,3:9,11:14, 16:20,22,23,25)],
                                 to=c(1,1,4,4,1,1,2,1,3,1,1,2,1,2,3,2,2,1,1,1))
    
    else
        stop("Error: level out of bounds")
    return(fasta_seq)
}

####################################################
#Public functions
####################################################
#' Makes hmm_fasta object
#'
#' @aliases convert reduce_FASTA
#' @param filename location of file
#' @param level value of f
#' @author Nina Wilkinson
#' @return \item{object}{hmm_fasta object}
#' @keywords character
#' @export
read_FASTA <- function(filename, level)
{
    fasta_seq=strsplit(paste(scan(filename,skip=1,what="character",comment.char=";"),collapse=""),"")[[1]]
    fasta_seq=convert(fasta_seq, level)
    join=c(1,length(fasta_seq))
    object=list(fasta_seq=fasta_seq,level=level,join=join)
    class(object) = "hmm_fasta"
    return(object)
}
#' Makes hmm_fasta object
#'
#' @aliases convert reduce_FASTA
#' @param name of vector of strings
#' @param level value of f
#' @author Nina Wilkinson
#' @return \item{object}{hmm_fasta object}
#' @keywords character
#' @export
read_FASTA_string <- function(name, level)
{
  fasta_seq = name
  join=c(1,length(fasta_seq))
  object=list(fasta_seq=fasta_seq,level=level,join=join)
  class(object) = "hmm_fasta"
  return(object)
}

#' Combines hmm_fasta objects
#'
#' @param protein1 hmm_fasta object
#' @param protein2 hmm_fasta object
#' @author Nina Wilkinson
#' @return \item{object}{hmm_fasta object}
#' @keywords character
#' @export
combine_FASTA = function(protein1,protein2){
    if (identical(class(protein1),class(protein2))!=TRUE) {
        stop("Objects are from different classes")
    }
    if (identical(protein1$level,protein2$level)!=TRUE){
        stop("Objects not have all same level")
    }
    vec=c(protein1$fasta_seq, protein2$fasta_seq)
    joins=c(protein1$join,protein2$join[2:length(protein2$join)]+protein1$join[length(protein1$join)])
    combine=list(fasta_seq=vec, level=protein1$level,join=joins)
    class(combine) <- "hmm_fasta"
    return(combine)
}
