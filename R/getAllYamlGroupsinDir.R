#' Searches the given directory for all readMe.yaml files
#'
#' Searches the given directory for all readMe.yaml files within all
#' subdirectories as well as their associated files their associated
#' "SNPTable_" files in that directory.
#'   NOTE: This function assumes that your YAML readme is in the same
#'          directory as its paired SNP CSV data. If thatis not the case,
#'          it will not work correctly.
#'
#' @param dir the directory that you are searching within
#' @param snptable.header the key text to grep when searching for the SNP csv.
#'     Default is "SNPTable_".
#' @param  yamlfile.header the key text to grep when searching for the readMe.yaml.
#'     Default is "readMe.yaml".
#'
#' @return A list with each element corresponding to a unique readMe.yaml file. For
#'   each of those elements, the first entry will be that yaml file and
#'   subsequent entries will be the associated SNP csv files. Will return NA
#'   if there are no paired SNP csv files found with a given readMe.yaml file.


GetAllYamlGroupsinDir <- function(dir, snptable.header = "SNPTable_",
                                  yamlfile.header = "readMe.yaml"){


  all.files <- list.files(dir, recursive = T);

  SNPtable.indexes <- grep(snptable.header, all.files) #Match files in header
  yaml.indexes <- grep(yamlfile.header, all.files)

  directory.SNPtable.matrix <- strsplit(all.files[SNPtable.indexes], "/")
  directory.yaml.matrix <- strsplit(all.files[yaml.indexes], "/")

  yaml.groups  <- list()  #Initialize Empty List

  for (i in 1:length(directory.yaml.matrix)) {
    # For every readMe.yaml file that was found...
    current.row <- directory.yaml.matrix[[i]]
    current.row.Length <- length(current.row)
    comparison.matrix <- t(sapply(directory.SNPtable.matrix,
                                  function(x) current.row == x))
    next.comparison.matrix <-  as.data.frame(
      comparison.matrix[,-ncol(comparison.matrix)])
    #Use all function to see if all comparisions true, if so, then part of group
    in.yaml.group <- apply(next.comparison.matrix, 1, function(x) all(x))

    #Paste Back into cumulative list
    SNPdata.tables.in.yamlgroup <- sapply(
      directory.SNPtable.matrix[in.yaml.group],
      function(x) paste(x, collapse = "/"))
    len.yamlgroup <- length(SNPdata.tables.in.yamlgroup )

    yaml.groups[[i]] <- paste(directory.yaml.matrix[[i]], collapse = "/")
    if (len.yamlgroup >0){
      yaml.groups[[i]][1:len.yamlgroup+1] <- SNPdata.tables.in.yamlgroup
    } else {
      yaml.groups[[i]][2] <- NA
    }
  }

  return(yaml.groups)

}
