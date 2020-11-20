server <- function(input, output) {

  volumes <- c(Home = fs::path_home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "plinkpath", roots = volumes)
  shinyFiles::shinyDirChoose(input, "workpath", roots = volumes)
  shinyFiles::shinyFileChoose(input, "vcfpath", roots = volumes)
  shinyFiles::shinyFileChoose(input, "csvpath", roots = volumes)

  output$plinkpath <- shiny::renderPrint({
    if (is.integer(input$plinkpath)) {
      cat("Previously downloaded & extracted binaries from https://www.cog-genomics.org/plink")
    } else {
      shinyFiles::parseDirPath(volumes, input$plinkpath)
    }})

  output$workpath <- shiny::renderPrint({
    if (is.integer(input$workpath)) {
      cat("Plink files will be written here.")
    } else {
      shinyFiles::parseDirPath(volumes, input$workpath)
    }})

  output$vcfpath <- shiny::renderPrint({
    if (is.integer(input$vcfpath)) {
      cat("VCF must be multisample and uncompressed.")
    } else {
      vpath <- shinyFiles::parseFilePaths(volumes, input$vcfpath)
      vpath$datapath
    }})

  output$csvpath <- shiny::renderPrint({
    if (is.integer(input$csvpath)) {
      cat("CSV table with phenotypes as columns.")
    } else {
      cpath <- shinyFiles::parseFilePaths(volumes, input$csvpath)
      cpath$datapath
    }})


  data <- shiny::reactive({
    path <- shinyFiles::parseFilePaths(volumes, input$csvpath)
    if (is.null(path$datapath)) { return() }
    else { read.csv(path$datapath) }
  })

  output$phenotypes <- DT::renderDataTable({
    if (is.integer(input$csvpath)) { data.frame(Waiting = "Load CSV please.") }
    else { data() }
  })


  output$sexcol_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("sexcolname", "Sex column:", choices = names(data()))
    } else { return() }
  })

  output$femalecode_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("femalecodeval", "Female encoding value:",
                  choices = unique(data()[input$sexcolname]))
    } else { return() }
  })

  output$malecode_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("malecodeval", "Male encoding value:",
                  choices = unique(data()[input$sexcolname]))
    } else { return() }
  })


  output$phencol_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("phencolname", "Column to Test:",
                         choices = names(data()))
    } else { return() }
  })


  output$seltest_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("plinkmode", "Plink Test:",
             choices = c("Association", "Model", "Linear", "Logistic"))
    } else { return() }
  })


  output$selcov_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath,
        input$plinkmode != "Association", input$plinkmode != "Model")
    if (!is.integer(input$csvpath)) {
      shiny::selectInput("covariables", "Linear/ Logistic Covariables:",
                         multiple = TRUE, choices = names(data()))
    } else { return() }
  })


  output$runbutton_ui <- shiny::renderUI({
    shiny::req(input$plinkpath, input$workpath, input$vcfpath, input$csvpath)
    if (!is.integer(input$csvpath)) {
      shiny::actionButton('execute', label = 'Run test!')
    } else { return() }
  })


  output$cmd_text <- shiny::eventReactive(input$execute, {
    ppath <- shinyFiles::parseDirPath(volumes, input$plinkpath)
    wpath <- shinyFiles::parseDirPath(volumes, input$workpath)
    vpath <- shinyFiles::parseFilePaths(volumes, input$vcfpath)
    if (Sys.info()["sysname"] == "Windows ") {
      command <- paste0(ppath, "/plink.exe")
    } else {
      command <- paste0(ppath, "/plink")
    }
    prefix <- paste0(wpath, "/plinkdata-", format(Sys.time(), "%y%m%d%H%M%S"),
                     "-", tolower(input$plinkmode))
    argsPre <- paste("--vcf", vpath$datapath,
                     "--recode --out", prefix, "\n", command, "--file",
                     prefix, "--make-just-fam --out", prefix)
    returnPre <- system2(command = command, args = argsPre,
                         stdout = TRUE, stderr = TRUE)

    fam <- read.csv(paste0(prefix, ".fam"), header = FALSE,
                    sep = " ", colClasses = "character")
    ped <- read.csv(paste0(prefix, ".ped"), header = FALSE,
                    sep = " ", colClasses = "character")
    sexcol <- gsub(input$femalecodeval, 2, data()[[input$sexcolname]])
    fam$V5 <- gsub(input$malecodeval, 1, sexcol)
    fam$V6 <- data()[[input$phencolname]]
    write.table(fam, file = paste0(prefix, ".fam"), quote = FALSE,
                na = "-9", row.names = FALSE, col.names = FALSE)
    ped[ , c(1:6)] <- fam
    write.table(ped, file = paste0(prefix, ".ped"), quote = FALSE,
                na = "-9", row.names = FALSE, col.names = FALSE)
    if (!is.null(input$covariables)) {
      covar_names <- unlist(lapply(strsplit(input$covariables, ","), trimws))
      write.table(cbind(fam$V1, fam$V2, data()[, input$covariables]),
                  file = paste0(prefix, ".cov"),
                  quote = FALSE, na = "-9", row.names = FALSE,
                  col.names = c("FID", "IID", covar_names))
    }

    if (input$plinkmode == "Association") {
      if (length(setdiff(unique(fam$V6), c(-9, 0, 1, 2))) > 0) {
        # quantitative phenotype
        argsTest <- paste("--assoc mperm=5000",
                          "--hwe 1e-8 midp --seed 9034",
                          "--file", prefix, "--out", prefix)
      } else {
        argsTest <- paste("--assoc fisher-midp mperm=5000",
                          "--hwe 1e-8 midp --seed 9034",
                          "--file", prefix, "--out", prefix)
      }
    } else if (input$plinkmode == "Model") {
      argsTest <- paste("--model fisher-midp mperm=5000",
                        "--hwe 1e-8 midp --seed 9034",
                        "--file", prefix, "--out", prefix)
    } else if (input$plinkmode == "Linear") {
      if (!is.null(input$covariables)) {
        covar <- paste("--covar", paste0(prefix, ".cov"), "--covar-name",
                       paste(covar_names, collapse = ", "))
      } else { covar <- " " }
      argsTest <- paste("--linear standard-beta genotypic interaction",
                        "--parameters 1-5, 7", covar,
                        "--file", prefix, "--out", prefix)
    } else if (input$plinkmode == "Logistic") {
      if (!is.null(input$covariables)) {
        covar <- paste("--covar", paste0(prefix, ".cov"), "--covar-name",
                       paste(covar_names, collapse = ", "))
      } else { covar <- " " }
      argsTest <- paste("--logistic --parameters 1-5, 7", covar,
                        "--file", prefix, "--out", prefix)
    }

    returnTest <- system2(command = command, args = argsTest,
                         stdout = TRUE, stderr = TRUE)

    return(paste(gsub("", "", c(returnPre, returnTest)), collapse = '\n'))
  })


    output$result_ui <- shiny::renderUI({
      shiny::req(input$workpath)
      if (!is.integer(input$workpath)) {
        wpath <- shinyFiles::parseDirPath(volumes, input$workpath)
        shiny::selectInput("result_file", "Choose result file:",
            list.files(file.path(wpath), pattern =
            "\\.(?:qassoc|assoc|model|linear|logistic|best|fisher|mperm)"),
          width = "100%")
    } else { return() }
    })

    output$results <- DT::renderDataTable({
      shiny::req(input$workpath, input$result_file)
      if (!is.integer(input$workpath) & !is.null(input$result_file)) {
        wpath <- shinyFiles::parseDirPath(volumes, input$workpath)
        read.table(paste0(wpath, "/", input$result_file), header = TRUE)
     } else { return() }
    })

}
