UFA_workflow <- function(spreadsheet) {
  ##
  tryCatch(stop(""), error = function(e) {NULL}) # To clear cache of error messages
  ##
  tryCatch(gc(), error = function(e) {NULL}, warning = function(w) {NULL})
  tryCatch(closeAllConnections(), error = function(e) {NULL}, warning = function(w) {NULL})
  ##
  on.exit({
    sillyErrors <- c("",
                     "subscript out of bounds",
                     "StopIteration", # `StopIteration` is generated by `stopCluster`
                     "cannot shut down device 1 (the null device)") # error by `dev.off()`
    ##
    errorMessages <- geterrmessage()
    xSillyErrors <- which(!(errorMessages %in% sillyErrors))
    if (length(xSillyErrors) > 0) {
      IPA_logRecorder(errorMessages[xSillyErrors])
      IPA_logRecorder("Stopped IDSL.UFA workflow!")
    }
    ##
    if (exists('.logIPA')) {
      rm(.logIPA, envir = .GlobalEnv)
    }
    ##
    tryCatch(gc(), error = function(e) {NULL}, warning = function(w) {NULL})
    tryCatch(closeAllConnections(), error = function(e) {NULL}, warning = function(w) {NULL})
    ##
  })
  ##
  ##############################################################################
  ##
  listPARAM <- UFA_xlsxAnalyzer(spreadsheet)
  if (!is.null(listPARAM)) {
    PARAM <- listPARAM[["PARAM"]]
    ############################################################################
    PARAM0001 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0001'), 2])
    if (PARAM0001 == "yes") {
      PARAM0002 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2])
      if (PARAM0002 == "yes") {
        PARAM_ECS <- listPARAM[["PARAM_ECS"]]
        UFA_enumerated_chemical_space(PARAM_ECS)
      }
      ##
      PARAM0003 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0003'), 2])
      if (PARAM0003 == "yes") {
        PARAM_FormSource <- listPARAM[["PARAM_FormSource"]]
        UFA_formula_source(PARAM_FormSource)
      }
    }
    ##
    PARAM0005 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0005'), 2])
    PARAM0006 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    ##
    if ((PARAM0005 == "yes") | (PARAM0006 == "yes")) {
      ##
      ##########################################################################
      ## To create log record for IDSL.UFA
      initiation_time_UFA <- Sys.time()
      timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
      input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]
      output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
      .logIPA <- NULL
      .logIPA <<- paste0(output_path, "/logUFA_performance.txt")
      IPA_logRecorder(paste0(rep("", 100), collapse = "="))
      IPA_logRecorder("Type <<< citation('IDSL.UFA') >>> for citing this R package in publications.")
      IPA_logRecorder(paste0("mzML/mzXML/netCDF:  ", input_path_hrms))
      IPA_logRecorder(paste0("OUTPUT:  ", output_path))
      IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
      IPA_logRecorder("Initiated IDSL.UFA workflow!")
      IPA_logRecorder(paste0(as.character(initiation_time_UFA), " ", timeZone))
      IPA_logRecorder("", allowedPrinting = FALSE)
      IPA_logRecorder("", allowedPrinting = FALSE)
      IPA_logRecorder(paste0(PARAM[, 1], "\t", PARAM[, 2]),  allowedPrinting = FALSE)
      IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
      ##
      ##########################################################################
      ##
    }
    ##
    ############################################################################
    ##
    if (PARAM0005 == "yes") {
      NPT <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2])
      ##
      addressIPDB <- PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]
      IPA_logRecorder(paste0("Loading the Isotopic Profiles DataBase (IPDB) from `", addressIPDB, "`!"))
      IPDB <- IDSL.IPA::loadRdata(addressIPDB)
      ## Temp
      if (is.null(IPDB[["MolecularFormula"]])) {
        stop(IPA_logRecorder("The selected IPDB is not consistent with this version of IDSL.UFA!"))
      }
      ##
      samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
      if (tolower(samples_string) == "all") {
        file_name_hrms <- dir(path = input_path_hrms)
        file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
      } else {
        file_name_hrms <- strsplit(samples_string, ";")[[1]]
      }
      LHRMS <- length(file_name_hrms)
      if (LHRMS == 0) {
        stop(IPA_logRecorder("ERROR!!! EMPTY HRMS FOLDER!!!"))
      }
      ##
      inputPathPeaklist <- PARAM[which(PARAM[, 1] == 'PARAM0011'), 2]
      peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
      peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
      L_PL <- length(peaklistFileNames)
      ##
      peaklist2HRMS <- gsub("^peaklist_|.Rdata$", "", peaklistFileNames)
      matchPeaklist2HRMS <- file_name_hrms[file_name_hrms %in% peaklist2HRMS]
      if (length(matchPeaklist2HRMS) != L_PL) {
        ndHRMSfiles <- setdiff(file_name_hrms, matchPeaklist2HRMS)
        print("Error!!! peaklist files are not available for the following HRMS file(s):")
        for (i in ndHRMSfiles) {
          print(i)
        }
        stop()
      }
      ##
      indexedIPApeaksCheck <- FALSE
      if (LHRMS == 1) {
        PARAM0013 <- PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]
        if (tolower(PARAM0013) != "all") {
          indexedIPApeaksCheck <- TRUE
          selectedIPApeaks <- tryCatch(eval(parse(text = paste0("c(", PARAM0013, ")"))), error = function(e){NULL})
          if (is.null(selectedIPApeaks)) {
            indexedIPApeaksCheck <- FALSE
          }
        }
      }
      ##
      output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
      if (!dir.exists(output_path_annotated_mf_tables)) {
        dir.create(output_path_annotated_mf_tables, recursive = TRUE)
      }
      ##
      parallelizationMode <- PARAM[which(PARAM[, 1] == 'PARAM0015'), 2]
      ##
      RTtolerance <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0016'), 2]), warning = function(w){NA})
      if (is.na(RTtolerance)) {
        retentionTimeCheck <- FALSE
      } else {
        retentionTimeCheck <- TRUE
        ##
        if (is.null(IPDB[["Retention Time"]])) {
          stop(IPA_logRecorder("Retention time match was selected in `PARAM0016` but <<< 'Retention Time' >>> values were not detected in the IPDB!"))
        }
        ##
        correctedRetentionTimeCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2]) == "yes") {TRUE} else {FALSE}
        if (correctedRetentionTimeCheck) {
          peak_alignment_folder <- PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]
          addressListCorrectedRTpeaklists <- paste0(peak_alignment_folder, "/listCorrectedRTpeaklists.Rdata")
          if (file.exists(addressListCorrectedRTpeaklists)) {
            listCorrectedRTpeaklists <- IDSL.IPA::loadRdata(addressListCorrectedRTpeaklists)
            IPA_logRecorder(paste0("Corrected retention times from `", addressListCorrectedRTpeaklists, "` are used for retention time match!"))
          }
        }
      }
      massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])
      maxNEME <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0019'), 2])
      minPCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])
      minNDCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])
      minRCS <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])
      scoreCoefficients <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])), error = function(e) {rep(1, 5)})
      maxAllowedNumberHits <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])
      ##
      if (maxAllowedNumberHits > 0) {
        dev.offCheck <- TRUE
        while (dev.offCheck) {
          dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
        }
        ##
        exportSpectraCheck <- TRUE
        outputProfileSpectra <- paste0(output_path, "/UFA_spectra")
        if (!dir.exists(outputProfileSpectra)) {
          dir.create(outputProfileSpectra, recursive = TRUE)
        }
        IPA_logRecorder("UFA_spectra comparison plots with theoretical isotopic profiles are stored in the `UFA_spectra` folder!")
      } else {
        exportSpectraCheck <- FALSE
        exportSpectraParameters <- NULL
      }
      ##
      libraryCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0025'), 2]) == "yes") {TRUE} else {FALSE}
      if (libraryCheck) {
        molecularFormulaLibrarySearchCheck <- TRUE
        IonPathways <- tryCatch(eval(parse(text = PARAM[which(PARAM[, 1] == 'PARAM0026'), 2])), error = function(e){"[M]"})
        IPA_logRecorder(paste0("Ionization pathway(s) to match with molecular formulas library (`PARAM0026`) include: ", paste0(IonPathways, collapse = ", ")))
        ##
        MFlibraryPath <- PARAM[which(PARAM[, 1] == 'PARAM0027'), 2]
        MFlibrary <- IDSL.IPA::loadRdata(MFlibraryPath)
      } else {
        molecularFormulaLibrarySearchCheck <- FALSE
      }
      ##
      IPA_logRecorder("Initiated molecular formula annotation on individual peaklists!")
      IPA_logRecorder("Individual annotated peaklist tables are stored in `.Rdata` and `.csv` formats in the `annotated_mf_tables` folder!")
      ##
      call_molecular_formula_annotator <- function(i) {
        peaklistFileName <- paste0("peaklist_", file_name_hrms[i], ".Rdata")
        peaklist <- IDSL.IPA::loadRdata(paste0(inputPathPeaklist, "/", peaklistFileName))
        ##
        if (!indexedIPApeaksCheck) {
          selectedIPApeaks <- seq(1, nrow(peaklist), 1)
        } else {
          selectedIPApeaks <- eval(parse(text = paste0("c(", PARAM0013, ")")))
        }
        ##
        outputer <- IDSL.IPA::IPA_MSdeconvoluter(inputHRMSfolderPath = input_path_hrms, MSfileName = file_name_hrms[i], MSlevel = 1)
        spectraList <- outputer[["spectraList"]]
        msPolarity <- outputer[["msPolarity"]]
        outputer <- NULL
        ##
        if (retentionTimeCheck) {
          if (correctedRetentionTimeCheck) {
            correctedRTpeaklist <- listCorrectedRTpeaklists[[peaklistFileName]]
          } else {
            correctedRTpeaklist <- NULL
          }
        } else {
          correctedRTpeaklist <- NULL
        }
        ##
        if (exportSpectraCheck) {
          exportSpectraParameters <- c(maxAllowedNumberHits, file_name_hrms[i], msPolarity, outputProfileSpectra)
        } else {
          exportSpectraParameters <- NULL
        }
        ##
        MolecularFormulaAnnotationTable <- molecular_formula_annotator(IPDB, spectraList, peaklist, selectedIPApeaks, massAccuracy,
                                                                       maxNEME, minPCS, minNDCS, minRCS, scoreCoefficients, RTtolerance,
                                                                       correctedRTpeaklist, exportSpectraParameters, number_processing_threads = NPT)
        ##
        if (!is.null(MolecularFormulaAnnotationTable)) {
          if (molecularFormulaLibrarySearchCheck) {
            MolecularFormulaAnnotationTable <- molecular_formula_library_search(MolecularFormulaAnnotationTable, MFlibrary, IonPathways, number_processing_threads = NPT)
          }
          ##
          save(MolecularFormulaAnnotationTable, file = paste0(output_path_annotated_mf_tables, "/MolecularFormulaAnnotationTable_", file_name_hrms[i], ".Rdata"))
          write.csv(MolecularFormulaAnnotationTable, file = paste0(output_path_annotated_mf_tables, "/MolecularFormulaAnnotationTable_", file_name_hrms[i], ".csv"), row.names = TRUE)
        }
        ##
        return()
      }
      ##
      ##########################################################################
      ##
      if (NPT == 1 | parallelizationMode == "peakmode") {
        ##
        progressBARboundaries <- txtProgressBar(min = 0, max = L_PL, initial = 0, style = 3)
        for (i in 1:L_PL) {
          null_variable <- tryCatch(call_molecular_formula_annotator(i),
                                    error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
          ##
          setTxtProgressBar(progressBARboundaries, i)
        }
        close(progressBARboundaries)
        ##
      } else if (parallelizationMode == "samplemode") {
        ##
        NPT0 <- NPT
        NPT <- 1
        ##
        osType <- Sys.info()[['sysname']]
        ##
        ########################################################################
        ##
        if (osType == "Windows") {
          ##
          clust <- makeCluster(NPT0)
          clusterExport(clust, setdiff(ls(), c("clust", "L_PL")), envir = environment())
          ##
          null_variable <- parLapply(clust, 1:L_PL, function (i) {
            tryCatch(call_molecular_formula_annotator(i),
                     error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
          })
          ##
          stopCluster(clust)
          ##
          ######################################################################
          ##
        } else {
          ##
          null_variable <- mclapply(1:L_PL, function (i) {
            tryCatch(call_molecular_formula_annotator(i),
                     error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
          }, mc.cores = NPT0)
          ##
          closeAllConnections()
          ##
        }
      }
      IPA_logRecorder("Sucessfully completed molecular formula annotation on individual peaklists!")
      IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
      IPDB <- NULL
      ##
      gc()
      closeAllConnections()
    }
    ##
    ############################################################################
    ##
    if (PARAM0006 == "yes") {
      aligned_molecular_formula_annotator(PARAM)
    }
    ##
    ############################################################################
    ##
    if ((PARAM0005 == "yes") | (PARAM0006 == "yes")) {
      completion_time_UFA <- Sys.time()
      IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
      required_time <- completion_time_UFA - initiation_time_UFA
      IPA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
      IPA_logRecorder(paste0(as.character(completion_time_UFA), " ", timeZone), allowedPrinting = FALSE)
      IPA_logRecorder("", allowedPrinting = FALSE)
      IPA_logRecorder("", allowedPrinting = FALSE)
      IPA_logRecorder("Completed IDSL.UFA workflow!")
      IPA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
    }
    ##
    ############################################################################
    ##
    PARAM0007 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0007'), 2])
    if (PARAM0007 == "yes") {
      PARAM_ScoreFunc <- listPARAM[["PARAM_ScoreFunc"]]
      UFA_score_function_optimization(PARAM_ScoreFunc)
    }
    ##
    ############################################################################
    ##
  }
}
