aligned_molecular_formula_annotator <- function(PARAM) {
  print("Initiated creating the aligned molecular formula annotated table!")
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == "PARAM0009"), 2])
  ##
  input_path_peak_Xcol <- PARAM[which(PARAM[, 1] == 'PARAM0025'), 2]
  peak_Xcol <- loadRdata(input_path_peak_Xcol)
  L_peaks <- dim(peak_Xcol)[1]
  L_samples <- dim(peak_Xcol)[2] - 2
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0014'), 2]
  output_path_annotated_mf_tables <- paste0(output_path, "/annotated_mf_tables")
  mf_table_list <- dir(path = output_path_annotated_mf_tables, pattern = ".Rdata")
  if (length(mf_table_list) != L_samples) {
    AnPL <- gsub("MolecularFormulaAnnotationTable_", "", mf_table_list)
    ColPL <- colnames(peak_Xcol)[3:(L_samples + 2)]
    MissedPL <- setdiff(ColPL, AnPL)
    print("Error!!! The following MolecularFormulaAnnotationTables are not avialable:")
    for (i in 1:length(MissedPL)) {
      print(MissedPL[i])
    }
    stop()
  }
  ##
  input_path_peak_property <- PARAM[which(PARAM[, 1] == 'PARAM0026'), 2]
  peak_property <- loadRdata(input_path_peak_property)
  if (dim(peak_property)[1] != L_peaks | dim(peak_property)[2] != (L_samples + 2)) {
    stop("Error!!! aligned peak property table and indexed peak table are not in the same size!!!")
  }
  ##
  N_top_candidates <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])
  N_candidate <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
  MF_Zcol <- array(rep(0, N_top_candidates*L_peaks*L_samples), dim = c(N_top_candidates, L_peaks, L_samples))
  ##
  adjust_freq_rank <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
  N_candidate3 <- N_candidate*3
  ##
  call_calculating_median_ranks <- function(j) {
    ID_freq_Rank <- rep(0, N_candidate3)
    ##
    entire_IDs <- MF_Zcol[, j, ]
    ID_rank <- do.call(rbind, lapply(1:L_samples, function(q) {
      cbind(entire_IDs[, q], seq(1, N_top_candidates))
    }))
    x_non0 <- which(ID_rank[, 1] != 0)
    if (length(x_non0) > 0) {
      ID_rank <- matrix(ID_rank[x_non0, ], ncol = 2)
      ##
      t_IDs <- sort(table(ID_rank[, 1]), decreasing = TRUE)
      max_k <- min(c(N_candidate, length(t_IDs)))
      t_freq <- as.numeric(t_IDs[1:max_k])
      t_id <- as.numeric(names(t_IDs[1:max_k]))
      ID_freq_Rank3 <- do.call(rbind, lapply(1:max_k, function(i) {
        x_id_t <- which(ID_rank[, 1] == t_id[i])
        med_rank <- median(ID_rank[x_id_t, 2]) # A median is calculated for the rank of candidate compounds across samples
        c(t_id[i], t_freq[i], med_rank)
      }))
      ID_freq_Rank3 <- matrix(ID_freq_Rank3[order(ID_freq_Rank3[, 3], decreasing = FALSE), ], ncol = 3)
      ID_freq_Rank3 <- matrix(ID_freq_Rank3[order(ID_freq_Rank3[, 2], decreasing = TRUE), ], ncol = 3)
      ##
      for (i in 1:max_k) {
        ID_freq_Rank[3*i - 2] <- ID_freq_Rank3[i, 1]
        ID_freq_Rank[3*i - 1] <- ID_freq_Rank3[i, 2]
        ID_freq_Rank[3*i] <- ID_freq_Rank3[i, 3]
      }
    }
    return(ID_freq_Rank)
  }
  ##
  call_creating_aligned_table <- function(i) {
    molf_IDs <- M_IDs[, (3*i - 2)]
    x_non0 <- which(molf_IDs != 0)
    molf <- rep(NA, L_peaks)
    if (length(x_non0) > 0) {
      matched_IDs_vec_db <- matrix(MolVecList_DB[molf_IDs[x_non0], ], ncol = L_Elements)
      molf_hill <- hill_molecular_formula_printer(Elements, matched_IDs_vec_db)
      molf[x_non0] <- molf_hill
    }
    sub_table <- cbind(molf, M_IDs[, (3*i - 1)], M_IDs[, 3*i])
    return(sub_table)
  }
  ##
  call_mz_rt_freq_median_peak_property <- function(i) {
    m_h <- 0
    x_h <- which(peak_property[i, 3:(L_samples + 2)] != 0)
    freq_h <- length(x_h)
    if (freq_h > 0) {
      m_h <- median(peak_property[i, (x_h + 2)])
    }
    c(freq_h, m_h)
  }
  ##
  if (tolower(adjust_freq_rank) == "yes") {
    ##
    call_freq_rank_table <- function(k) {
      freq <- M_IDs[, (3*k - 1)]
      rank <- M_IDs[, (3*k)]
      sqrt(freq)/rank           # To adjust ranking and frequencies
    }
    ##
    call_M_IDs2 <- function(k) {
      x_0 <- which(freq_rank_table[k, ] > 0)
      L_x_0 <- length(x_0)
      if (L_x_0 > 1) {
        order_rank <- order(freq_rank_table[k, x_0], decreasing = TRUE)
        M_ID_ordered <- do.call(cbind, lapply(order_rank, function (i) {
          cbind(M_IDs[k, (3*i - 2)], M_IDs[k, (3*i - 1)], M_IDs[k, 3*i])
        }))
        if (L_x_0 < N_candidate) {
          M_ID_ordered <- c(M_ID_ordered, rep(0, (N_candidate3 - 3*L_x_0)))
        }
      } else {
        M_ID_ordered <- M_IDs[k, ]
      }
      return(M_ID_ordered)
    }
  }
  ##
  print("Initiated matching peak IDs!")
  progressBARboundaries <- txtProgressBar(min = 1, max = L_samples, initial = 1, style = 3)
  ##
  osType <- Sys.info()[['sysname']]
  if (osType == "Windows") {
    clust <- makeCluster(number_processing_threads)
    registerDoSNOW(clust)
    ##
    for (i in 1:L_samples) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      peak_table_id <- peak_Xcol[, (i + 2)]
      MolecularFormulaAnnotationTable <- loadRdata(paste0(output_path_annotated_mf_tables, "/", mf_table_list[i]))
      matched_peak_ids <- as.numeric(MolecularFormulaAnnotationTable[, 1])
      matched_mf_ids <- as.numeric(MolecularFormulaAnnotationTable[, 2])
      x_peak_ids <- which(peak_table_id%in%unique(matched_peak_ids) == TRUE)
      L_x_peak_ids <- length(x_peak_ids)
      ##
      if (L_x_peak_ids > 0) {
        mjx <- foreach(j = x_peak_ids, .verbose = FALSE) %dopar% {
          x_j <- which(matched_peak_ids == peak_table_id[j])
          matched_mf_ids[x_j]
        }
        ##
        for (j in 1:L_x_peak_ids) {
          x_mf_id <- mjx[[j]]
          max_k <- min(c(N_top_candidates, length(x_mf_id)))
          for (k in 1:max_k) {
            MF_Zcol[k, x_peak_ids[j], i] <- x_mf_id[k]
          }
        }
      }
    }
    close(progressBARboundaries)
    print("Completed matching peak IDs!")
    ##
    print("Initiated calculating median ranks!")
    M_IDs <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
      call_calculating_median_ranks(k)
    }
    MF_Zcol <- 0
    print("Completed calculating median ranks!")
    ##
    if (tolower(adjust_freq_rank) == "yes") {
      print("Initiated adjusting frequencies and ranks!")
      ##
      freq_rank_table <- foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
        call_freq_rank_table(k)
      }
      #
      freq_rank_table[is.nan(freq_rank_table)] <- 0
      ##
      M_IDs <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_M_IDs2(k)
      }
      print("Completed adjusting frequencies and ranks!")
    }
    ##
    address_sav_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
    print("Loading the isotopic profiles database!")
    IPDB <- loadRdata(address_sav_IPDB)
    MolVecList0 <- IPDB[[2]]
    IPDB <- 0
    Elements <- MolVecList0[[1]]
    MolVecList_DB <- MolVecList0[[2]]
    L_Elements <- length(Elements)
    print("Initiated creating the aligned table!")
    aligned_molecular_formula <- foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
      call_creating_aligned_table(k)
    }
    print("Completed creating the aligned table!")
    ##
    print("Initiated processing the peak property table!")
    ##
    IPA_Xcol <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
      c(peak_Xcol[k, 1:2], (length(which(peak_Xcol[k, ] > 0)) - 2))
    }
    ##
    mz_rt_freq_median_peak_property <- foreach(k = 1:L_peaks, .combine = 'rbind', .verbose = FALSE) %dopar% {
      call_mz_rt_freq_median_peak_property(k)
    }
    ##
    title_mat <-  foreach(k = 1:N_candidate, .combine = 'cbind', .verbose = FALSE) %dopar% {
      cbind(paste0("IonFormula_", k), paste0("Frequency_", k), paste0("MedianRank_", k))
    }
    ##
    stopCluster(clust)
  }
  if (osType == "Linux") {
    ##
    for (i in 1:L_samples) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      peak_table_id <- peak_Xcol[, (i + 2)]
      MolecularFormulaAnnotationTable <- loadRdata(paste0(output_path_annotated_mf_tables, "/", mf_table_list[i]))
      matched_peak_ids <- as.numeric(MolecularFormulaAnnotationTable[, 1])
      matched_mf_ids <- as.numeric(MolecularFormulaAnnotationTable[, 2])
      x_peak_ids <- which(peak_table_id%in%unique(matched_peak_ids) == TRUE)
      L_x_peak_ids <- length(x_peak_ids)
      ##
      if (L_x_peak_ids > 0) {
        mjx <- mclapply(x_peak_ids, function(j) {
          x_j <- which(matched_peak_ids == peak_table_id[j])
          matched_mf_ids[x_j]
        }, mc.cores = number_processing_threads)
        ##
        for (j in 1:L_x_peak_ids) {
          x_mf_id <- mjx[[j]]
          max_k <- min(c(N_top_candidates, length(x_mf_id)))
          for (k in 1:max_k) {
            MF_Zcol[k, x_peak_ids[j], i] <- x_mf_id[k]
          }
        }
      }
    }
    close(progressBARboundaries)
    print("Completed matching peak IDs!")
    ##
    print("Initiated calculating median ranks!")
    M_IDs <- do.call(rbind, mclapply(1:L_peaks, function (k) {
      call_calculating_median_ranks(k)
    }, mc.cores = number_processing_threads))
    MF_Zcol <- 0
    print("Completed calculating median ranks!")
    ##
    if (tolower(adjust_freq_rank) == "yes") {
      print("Initiated adjusting frequencies and ranks!")
      ##
      freq_rank_table <- do.call(cbind, mclapply(1:N_candidate, function (k) {
        call_freq_rank_table(k)
      }, mc.cores = number_processing_threads))
      #
      freq_rank_table[is.nan(freq_rank_table)] <- 0
      ##
      M_IDs <- do.call(rbind, mclapply(1:L_peaks, function (k) {
        call_M_IDs2(k)
      }, mc.cores = number_processing_threads))
      print("Completed adjusting frequencies and ranks!")
    }
    ##
    address_sav_IPDB <- PARAM[which(PARAM[, 1] == "PARAM0004"), 2]
    print("Loading the isotopic profiles database!")
    IPDB <- loadRdata(address_sav_IPDB)
    MolVecList0 <- IPDB[[2]]
    IPDB <- 0
    Elements <- MolVecList0[[1]]
    MolVecList_DB <- MolVecList0[[2]]
    L_Elements <- length(Elements)
    print("Initiated creating the aligned table!")
    aligned_molecular_formula <- do.call(cbind, mclapply(1:N_candidate, function (k) {
      call_creating_aligned_table(k)
    }, mc.cores = number_processing_threads))
    print("Completed creating the aligned table!")
    ##
    print("Initiated processing the peak property table!")
    ##
    IPA_Xcol <- do.call(rbind, mclapply(1:L_peaks, function (k) {
      c(peak_Xcol[k, 1:2], (length(which(peak_Xcol[k, ] > 0)) - 2))
    }, mc.cores = number_processing_threads))
    ##
    mz_rt_freq_median_peak_property <- do.call(rbind, mclapply(1:L_peaks, function (k) {
      call_mz_rt_freq_median_peak_property(k)
    }, mc.cores = number_processing_threads))
    ##
    title_mat <- do.call(cbind, mclapply(1:N_candidate, function(k) {
      cbind(paste0("IonFormula_", k), paste0("Frequency_", k), paste0("MedianRank_", k))
    }, mc.cores = number_processing_threads))
    ##
    closeAllConnections()
  }
  ##
  aligned_molecular_formula <- data.frame(cbind(IPA_Xcol, mz_rt_freq_median_peak_property, aligned_molecular_formula))
  ##
  ppn1 <- strsplit(input_path_peak_property, "/")[[1]]
  ppn <- ppn1[length(ppn1)]
  peak_property_name <- gsub(".Rdata", "", ppn)
  #
  title_mat <- cbind("m/z", "RT", "IPA detection frequency", paste0(peak_property_name, " frequency"), paste0("median ", peak_property_name), title_mat)
  colnames(aligned_molecular_formula) <- title_mat
  rownames(aligned_molecular_formula) <- c()
  print("Completed processing the peak property table!")
  output_path_aligned_table <- paste0(output_path, "/aligned_molecular_formula_table")
  if (!dir.exists(output_path_aligned_table)) {
    dir.create(output_path_aligned_table)
  }
  print("Initiated saving the aligned molecular formula table!")
  save(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.Rdata"))
  write.csv(aligned_molecular_formula, file = paste0(output_path_aligned_table, "/aligned_molecular_formula.csv"))
  print("Completed saving the aligned molecular formula table!")
}
