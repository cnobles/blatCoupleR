#' Pander out head portions of objects for checking.
#' @param object object which will be coerced into a data.frame and will be 
#' passed to head().
#' @param title character Title to be included above table
#' @param caption character Caption to be included beneath table.
#' @param row.names logical TRUE includes while FALSE (default) drops row names.
#' @author Christopher Nobles, Ph.D.
#' 

panderHead <- function(object, title = NULL, caption = NULL, row.names = FALSE){
  stopifnot(!class(object) == "list")
  if(!is.null(title)) pandoc.title(title)
  
  if(!row.names){
    df <- as.data.frame(head(object), row.names = NULL)
  }else{
    df <- as.data.frame(head(object))
  }
  pandoc.table(df, style = "simple", split.tables = Inf, caption = caption)
}