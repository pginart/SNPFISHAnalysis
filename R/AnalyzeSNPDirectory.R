#' Does a rough, first pass SNP FISH analysis for a given directory
#'
#' @param masterdir directory to search for all sub yaml-files for standard
#'        analysis
#' @param outdir directory in which to output SNP plots
#' @param SNPTableFlag Will produce plots on SNP tables if they exist. Default
#'  is false.
AnalyzeSNPDirectory <-
  function(masterdir, outdir, SNPTableFlag = FALSE) {
    yaml.groups <- GetAllYamlGroupsinDir(masterdir)
    for (i in 1:length(yaml.groups)) {
      data.Files <-
        ExtractFilesFromYamlGroup(yaml.groups[[i]], masterdir)

      exp.data <- data.Files[[1]]
      guide.table <- data.Files[[2]]

      MakeGuideTablePlots(exp.data, guide.table, outdir)

      if (SNPTableFlag) {
        tmr.table <- data.Files[[3]]
        cy.table <- data.Files[[4]]

        SNPDir <- paste(outdir, "/SNPChannelData", sep = "")
        if (!dir.exists(SNPDir)) {
          dir.create(SNPDir)
        }

        MakeStackedCellCountPlot(exp.data, tmr.table, SNPDir)
        MakeTotalCellCountPie(exp.data, tmr.table, SNPDir)
        MakeStackedCellCountPlot(exp.data, cy.table, SNPDir)
        MakeTotalCellCountPie(exp.data, cy.table, SNPDir)

      }

    }

  }
