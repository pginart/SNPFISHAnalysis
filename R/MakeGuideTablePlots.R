#' Makes a standard set of graphs for guide table
#'
#' @param exp.Data experiment data extracted from YAML file
#' @param guide.table SNP Table for guide probe
#' @param outdir directory in which to output SNP plots
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom scales percent
MakeGuideTablePlots <- function(exp.data, guide.table, outdir){
  MakeStackedCellCountPlot(exp.data, guide.table, outdir)
  MakeTotalCellCountPie(exp.data, guide.table, outdir)
  MakeOverallImbalanceBar(exp.data, guide.table, outdir)
  MakeAlleleScatterPlot(exp.data, guide.table, outdir)
}


