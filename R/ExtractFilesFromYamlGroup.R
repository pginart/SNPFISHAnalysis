#' Loads Yaml file and SNPTables from a YamlGroup
#'
#' SNPTables are extracted from a YamlGroup according to the headers used in
#' the file names. Modify defaults as required. Objects not present will have
#' their places held by NA.
#'
#' @param yaml.group the list of directories of a yaml file and its associated
#'  SNP tables
#' @param masterDir directory where yamlGroup was called
#' @param guide.header tag to search for SNP guide table. The default is
#' "_guide_".
#' @param cy.header tag to search for tmr or SNP 1 table. The default is
#' "_tmr_"
#' @param tmr.header tag to search for cy or SNP 2 table. The default is
#' "_cy_"
#'
#' @return A list ordered in the following way:
#'  The first object will the the experimental data loaded from the Yaml file
#'  The second object will the guide SNP data table
#'  The third object will be the tmr SNP data table.
#'  The fourth object will be the cy SNP data table.

ExtractFilesFromYamlGroup <- function(yaml.group, masterDir,  guide.header = "_guide_",
                                      cy.header = "_cy_", tmr.header = "_tmr_"){
  guide.index <- grepl(guide.header, yaml.group)
  cy.index <- grepl(cy.header, yaml.group)
  tmr.index <- grepl(tmr.header, yaml.group)


  # Get table directories if they exist
  yaml.dir <- paste(masterDir, yaml.group[1], sep = "")

  if (!any(cy.index)){
    error.message <- paste("No SNP data table found matching ",
                           cy.header, sep = "" )
    stop(error.message)
    cy.table.dir <- NA

  } else {
    cy.table.dir <- paste(masterDir, yaml.group[cy.index], sep = "")
  }
  if (!any(tmr.index)){
    error.message <- paste("No SNP data table found matching ",
                           tmr.header, sep = "" )
    stop(error.message)
    tmr.table.dir <- NA

  } else {
    tmr.table.dir <- paste(masterDir, yaml.group[tmr.index], sep = "")
  }

  if (!any(guide.index)){
    error.message <- paste("No SNP data table found matching ",
                           guide.header, sep = "" )
    stop(error.message)
  } else {
    guide.table.dir <- paste(masterDir, yaml.group[guide.index], sep = "")
  }


  #Load data files
  exp.data <- yaml::yaml.load_file(yaml.dir)
  guide.table <- dplyr::tbl_df(read.csv(guide.table.dir , stringsAsFactors = T))
  guide.table <- ReorderSNPLevels(guide.table)

  if (is.na(cy.table.dir)){
    cy.table <- NA
  } else {
    cy.table <- dplyr::tbl_df(read.csv(cy.table.dir , stringsAsFactors = T))
    cy.table <- ReorderSNPLevels(cy.table)
  }

  if (is.na(tmr.table.dir)){
    tmr.table <- NA
  } else {
    tmr.table <- dplyr::tbl_df(read.csv(tmr.table.dir , stringsAsFactors = T))
    tmr.table <- ReorderSNPLevels(tmr.table)
  }

  List.to.return <- list()
  List.to.return[[1]] <- exp.data
  List.to.return[[2]] <- guide.table
  List.to.return[[3]] <- tmr.table
  List.to.return[[4]] <- cy.table

  return(List.to.return)
}


