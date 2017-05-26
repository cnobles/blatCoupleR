#' Write output files based on file extensions.
#' @param object any object that can be coerced to a data.frame with 
#' as.data.frame() methods. Can also be a list of objects with the same 
#' requirements. In the latter case, the file name will be parsed into a prefix
#' and file extension, and the name of the object will be inserted between them
#' when writen to either a csv or tsv file.
#' @param file character File name with extension. Supported formats include:
#' csv, tsv, rds, and RData.
#' @author Christopher Nobles, Ph.D.

writeOutputFile <- function(object, file){
  fileType <- stringr::str_extract(file, "[\\w]+$")
  objectClass <- class(object)
  if(fileType == "csv" & objectClass != "list"){
    write.csv(
      as.data.frame(object, row.names = NULL), 
      file = file, 
      quote = FALSE, 
      row.names = FALSE)
  }else if(fileType == "tsv" & objectClass != "list"){
    write.table(
      as.data.frame(object, row.names = NULL), 
      file = file, 
      sep = "\t",
      quote = FALSE, 
      row.names = FALSE)
  }else if(fileType == "csv" & objectClass == "list"){
    prefix <- gsub("csv$", "", file)
    if(is.null(names(object))) names(object) <- 1:length(object)
    objectNames <- names(object)
    fileNames <- paste0(prefix, objectNames, ".csv")
    null <- mapply(function(ob, fileName){
        ob <- as.data.frame(ob, row.names = NULL)
        write.csv(ob, fileName, quote = FALSE, row.names = FALSE)
      }, ob = object, fileName = fileNames)
  }else if(fileType == "tsv" & objectClass == "list"){
    prefix <- gsub("tsv$", "", file)
    objectNames <- names(object)
    fileNames <- paste0(prefix, objectNames, ".tsv")
    null <- mapply(function(ob, fileName){
        ob <- as.data.frame(ob, row.names = NULL)
        write.table(ob, fileName, sep = "\t", quote = FALSE, row.names = FALSE)
      }, ob = object, fileName = fileNames)
  }else if(fileType == "rds"){
    saveRDS(object, file)
  }else if(fileType == "RData"){
    save(object, file)
  }else{
    stop("Unsupported file type. Supported: csv, tsv, rds, and RData.")
  }
}