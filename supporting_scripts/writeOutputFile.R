#' Write output files based on file extensions.
#' @param object
#' @param file character File name with extension. Supported formats include:
#' csv, tsv, rds, and RData.
#' @author Christopher Nobles, Ph.D.

writeOutputFile <- function(object, file){
  fileType <- stringr::str_extract(file, "[\\w]+$")
  if(fileType == "csv"){
    write.csv(
      as.data.frame(object), 
      file = file, 
      quote = FALSE, 
      row.names = FALSE)
  }else if(fileType == "tsv"){
    write.table(
      as.data.frame(object), 
      file = file, 
      sep = "\t",
      quote = FALSE, 
      row.names = FALSE)
    )
  }else if(fileType == "rds"){
    saveRDS(object, file)
  }else if(fileType == "RData"){
    save(object, file)
  }else{
    stop("Unsupported file type. Supported: csv, tsv, rds, and RData.")
  }
}