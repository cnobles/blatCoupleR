#' read psl.gz files, assuming psl gz files don't have column header
#' @param pslFile character vector of file name(s).
#' @param toNull  character vector of column names to remove.
#' @return data.frame, data.table of the psl table
#' @author Christopher Nobles, Ph.D. and Yinghua Wu, Ph.D.

readPSL <- function(pslFile, toNull = NULL) {
  stopifnot(require("data.table"))
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  psl <- lapply(pslFile, function(f){
    message("Reading ",f)
    suppressMessages(
      try(data.table::fread(paste("zcat", f), sep = "\t"), silent = TRUE))
  })
  psl <- lapply(psl, function(p){
    if(class(p) == "try-error" & grepl("File is empty:", p)){
      p <- read.table(text = "", col.names = cols, colClasses = cols.class)
      return(as.data.table(p))
    }else if(class(p) == "try-error"){
      stop("Error in loading psl files. Check output from alignment (blat).")
    }else{
      return(p)
    }
  })
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- cols
  
  if(length(toNull)>0) psl[, toNull] <- NULL
  
  return(as.data.frame(psl))
}