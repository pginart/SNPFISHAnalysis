#' Make catagory names friendly to ggplot2
#'
#' Changes a string so it is enclosed in back-ticks.
#' This can be used to make column names that have spaces (blanks)
#' or non-letter characters acceptable to ggplot2.
#' This version of the function is vectorized with sapply.
#'
#' @param x string to make acceptable to ggplot2
#' @return x string made acceptable to ggplot2
#'
#' @examples
#'  ggname("B6 Allele")
ggname <- function(x) {
  if (class(x) != "character") {
    return(x)
  }
  y <- sapply(x, function(s) {
    if (!grepl("^`", s)) {
      s <- paste("`", s, sep="", collapse="")
    }
    if (!grepl("`$", s)) {
      s <- paste(s, "`", sep="", collapse="")
    }
  }
  )
  y
}
