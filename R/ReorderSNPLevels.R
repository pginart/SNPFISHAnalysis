#' Reorders the label's levels in table such that undetected spots and
#' three color are stacked on top of the SNP labels
#'
#' @param snp.table the standard SNP table whose levels will be reordered
#' @param undetec.label label for undetected spots. Default is "undetec".
#' @param threecolor.label label used for three-color spots. Default is "3-color".
#
#' @return snp.table with levels reordered
ReorderSNPLevels <- function(snp.table, undetec.label = "undetec",
                             threecolor.label = "3-color") {

  initial.levels <- levels(snp.table$labels)

  undetec.index <-  which(initial.levels == undetec.label)
  three.color.index <-  which(initial.levels == threecolor.label)
  initial.levels <- initial.levels[-undetec.index]
  initial.levels <- initial.levels[-three.color.index]
  final.levels <- append(initial.levels, c(undetec.label, threecolor.label),
                         after = length(initial.levels))

  snp.table$labels <- factor(snp.table$labels, levels = final.levels)
  return(snp.table)

}
