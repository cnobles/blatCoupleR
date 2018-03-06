#' Read key files from multiple formats.
#' @param keyFile character File name (path/to/name.format).
#' @param format character Acceptable formats include: "csv", "tsv", "rds", and 
#' "RData".
#' @author Christopher Nobles, Ph.D.
#'

readKeyFile <- function(keyFile, format){
  cols <- c("readNames", "seqID")
  cols.class <- c("character", "character")
  if(format == "csv"){
    file <- try(data.table::fread(keyFile, sep = ",", data.table = FALSE))
    if(class(file) == "try-error" & grepl("File is empty:", file[1])){
      file <- read.table(text = "", col.names = cols, colClasses = cols.class)
    }
    return(file)
  }else if(format == "tsv"){
    file <- try(data.table::fread(keyFile, sep = "\t", data.table = FALSE))
    if(class(file) == "try-error" & grepl("no lines available", file[1])){
      file <- read.table(text = "", col.names = cols, colClasses = cols.class)
    }
    return(file)
  }else if(format == "rds"){
    file <- readRDS(keyFile)
    return(as.data.frame(file))
  }else if(format == "RData"){
    env <- new.env()
    load(keyFile, envir = env)
    if(length(env) > 1) stop("More than one object present in keyFile (RData).")
    return(env[[ls(env)]])
  }else{
    stop("Key file not in a supported file type, convert to csv, tsv, rds, or RData.")
  }
}
