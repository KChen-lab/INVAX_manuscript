###########################################################
############### Modified tradeSeq functions ###############
###########################################################


assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
                           FUN = "/",
                           STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        stats::rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}

get_offset <- function(offset, counts) {
  if (is.null(offset)) {
    nf <- try(edgeR::calcNormFactors(counts), silent = TRUE)
    if (is(nf, "try-error")) {
      message("TMM normalization failed. Will use unnormalized library sizes",
              "as offset.\n")
      nf <- rep(1,ncol(counts))
    }
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
    if(any(libSize == 0)){
      message("Some library sizes are zero. Offsetting these to 1.\n")
      offset[libSize == 0] <- 0
    }
  }
  return(offset)
}

findKnots <- function(nknots, pseudotime, wSamp) {
  # Easier to recreate them all here than to pass them on
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  
  # Get the times for the knots
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }
  
  knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- stats::quantile(t1[l1 == 1],
                                probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (max(which(dupId)) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }
  
  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- vapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    }, FUN.VALUE = 1)
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }
  
  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- min(tAll)
  knots[nknots] <- max(tAll)
  
  knotList <- lapply(seq_len(ncol(pseudotime)), function(i){
    knots
  })
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))
  
  return(knotList)
}

fitGAM <- function(counts, U = NULL, pseudotime, cellWeights,
                   conditions,
                   genes = seq_len(nrow(counts)),
                   weights = NULL, offset = NULL, nknots = 6, verbose = TRUE,
                   parallel = FALSE, BPPARAM = BiocParallel::bpparam(),
                   aic = FALSE, control = mgcv::gam.control(), sce = TRUE,
                   family = "nb", gcv = FALSE){
  
  if(!is.null(conditions)){
    message("Fitting lineages with multiple conditions. This method has ",
            "been tested on a couple of datasets, but is still in an ",
            "experimental phase.")
  }
  
  if (is(genes, "character")) {
    if (!all(genes %in% rownames(counts))) {
      stop("The genes ID is not present in the models object.")
    }
    if(any(duplicated(genes))){
      stop("The genes vector contains duplicates.")
    }
    id <- match(genes, rownames(counts))
  } else {
    id <- genes
  }
  
  # if (parallel) {
  #   BiocParallel::register(BPPARAM)
  #   if (verbose) {
  #     # update progress bar 40 times
  #     BPPARAM$tasks = as.integer(40)
  #     # show progress bar
  #     BPPARAM$progressbar = TRUE
  #   }
  # }
  
  
  # Convert pseudotime and weights to matrices if need be
  # if (is.null(dim(pseudotime))) {
  #   pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
  # }
  if (is.null(dim(cellWeights))) {
    cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
  }
  
  # .checks(pseudotime, cellWeights, U, counts, conditions, family)
  
  wSamp <- assignCells(cellWeights)
  
  # define pseudotime for each lineage
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  
  # offset
  offset <- get_offset(offset, counts)
  
  # fit model
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
  }
  
  ## Get the knots
  knotList <- findKnots(nknots, pseudotime, wSamp)
  
  ## fit NB GAM
  ### Actually fit the model ----
  teller <- 0
  converged <- rep(TRUE, length(genes))
  counts_to_Gam <- function(y) {
    teller <<- teller + 1
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if (!is.null(weights)) weights <- weights[teller,]
    if (!is.null(dim(offset))) offset <- offset[teller,]
    if(is.null(conditions)){
      smoothForm <- stats::as.formula(
        paste0("y ~ ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
               }, FUN.VALUE = "formula"),
               collapse = "+"))
      )
    } else {
      for(jj in seq_len(ncol(pseudotime))){
        for(kk in seq_len(nlevels(conditions))){
          # three levels doesn't work. split it up and loop over both conditions and pseudotime
          # to get a condition-and-lineage-specific smoother. Also in formula.
          lCurrent <- get(paste0("l", jj))
          id1 <- which(lCurrent == 1)
          lCurrent[id1] <- ifelse(conditions[id1] == levels(conditions)[kk], 1, 0)
          assign(paste0("l", jj, "_", kk), lCurrent)
        }
      }
      smoothForm <- stats::as.formula(
        paste0("y ~ ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste(vapply(seq_len(nlevels(conditions)), function(kk){
                   paste0("s(t", ii, ", by=l", ii, "_", kk,
                          ", bs='cr', id=1, k=nknots)")
                 }, FUN.VALUE = "formula"),
                 collapse = "+")
               }, FUN.VALUE = "formula"),
               collapse="+"))
      )
    }
    # fit smoother, catch errors and warnings
    s <- mgcv::s
    m <- suppressWarnings(try(withCallingHandlers({
      mgcv::gam(smoothForm, family = family, knots = knotList, weights = weights,
                control = control)},
      error = function(e){ #if errors: return try-error class
        converged[teller] <<- FALSE
        return(structure("Fitting errored",
                         class = c("try-error", "character")))
      },
      warning = function(w){ #if warning: set converged to FALSE
        converged[teller] <<- FALSE
      }), silent=TRUE))
    return(m)
  }
  
  #### fit models
  # if (parallel) {
  #   gamList <- BiocParallel::bplapply(
  #     as.data.frame(t(as.matrix(counts)[id, ])),
  #     counts_to_Gam, BPPARAM = BPPARAM
  #   )
  # } else {
  if (verbose) {
    gamList <- pbapply::pblapply(
      # as.data.frame(t(as.matrix(counts)[id, ])),
      as.data.frame(t(as.matrix(counts))),
      counts_to_Gam
    )
  } else {
    gamList <- lapply(
      # as.data.frame(t(as.matrix(counts)[id, ])),
      as.data.frame(t(as.matrix(counts))),
      counts_to_Gam
    )
  }
  # }
  
  ### output
  if (aic) { # only return AIC
    # return(unlist(lapply(gamList, function(x){
    #   if (class(x)[1] == "try-error") return(NA)
    #   x$aic
    # })))
    aicVals <- unlist(lapply(gamList, function(x){
      if (class(x)[1] == "try-error") return(NA)
      x$aic
    }))
    if (gcv) {
      gcvVals <- unlist(lapply(gamList, function(x){
        if (class(x)[1] == "try-error") return(NA)
        x$gcv.ubre
      }))
      return(list(aicVals, gcvVals))
    } else return(aicVals)
  }
  
  # if (sce) { #tidy output: also return X
  # tidy smoother regression coefficients
  betaAll <- lapply(gamList, function(m) {
    if (is(m, "try-error")) {
      beta <- NA
    } else {
      beta <- matrix(stats::coef(m), ncol = 1)
      rownames(beta) <- names(stats::coef(m))
    }
    return(beta)
  })
  betaAllDf <- data.frame(t(do.call(cbind,betaAll)))
  rownames(betaAllDf) <- rownames(counts)[id]
  
  # list of variance covariance matrices
  SigmaAll <- lapply(gamList, function(m) {
    if (is(m, "try-error")) {
      Sigma <- NA
    } else {
      Sigma <- m$Vp
    }
    return(Sigma)
  })
  
  # Get X, dm and knotPoints
  element <- min(which(!is.na(SigmaAll)))
  m <- gamList[[element]]
  X <- stats::predict(m, type = "lpmatrix")
  dm <- m$model[, -1]
  knotPoints <- m$smooth[[1]]$xp
  
  # # return output
  # return(list(beta = betaAllDf,
  #             Sigma = SigmaAll,
  #             X = X,
  #             dm = dm,
  #             knotPoints = knotPoints,
  #             converged = converged)
  # )
  
  gamOutput <- list(beta = betaAllDf,
                    Sigma = SigmaAll,
                    X = X,
                    dm = dm,
                    knotPoints = knotPoints,
                    converged = converged)
  
  # return SingleCellExperiment object
  sc <- SingleCellExperiment(assays = list(counts = counts))
  # slingshot info
  SummarizedExperiment::colData(sc)$crv <- S4Vectors::DataFrame(
    pseudotime = pseudotime,
    cellWeights = cellWeights)
  # tradeSeq gene-level info
  df <- tibble::enframe(gamOutput$Sigma, value = "Sigma")
  df$beta <- tibble::tibble(beta = gamOutput$beta)
  df$converged <- gamOutput$converged
  suppressWarnings(rownames(df) <- rownames(counts))
  SummarizedExperiment::rowData(sc)$tradeSeq <- df
  # tradeSeq cell-level info
  if(is.null(conditions)){
    SummarizedExperiment::colData(sc)$tradeSeq <-
      tibble::tibble(X = gamOutput$X, dm = gamOutput$dm)
  } else {
    SummarizedExperiment::colData(sc)$tradeSeq <-
      tibble::tibble(X = gamOutput$X, dm = gamOutput$dm,
                     conditions = conditions)
  }
  # metadata: tradeSeq knots
  S4Vectors::metadata(sc)$tradeSeq <-
    list(knots = gamOutput$knotPoints)
  # return(sc)
  # } else {
  #   return(gamList)
  # }
  return(list(gamList = gamList,
              sc = sc))
}


find_conditions <- function(models) {
  if (is(models, "list")) {
    stop("Argument models must be a SingleCellExperiment object.",
         "Please refit the smoothers as such, see ?fitGAM.")
  } else if (is(models, "SingleCellExperiment")) {
    condPresent <- suppressWarnings({
      !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
    })
    if(!condPresent) stop("The models were not fitted with multiple conditions.")
    conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditions)
    if (nConditions == 1) stop("This can only be run with multiple conditions.")
  }
  return(list("conditions" = conditions, "nConditions" = nConditions))
}


construct_contrast_matrix_conditions <- function(models, X, conditions, 
                                                 nConditions, nLineages, 
                                                 nKnots, knots) {
  # do statistical test for every model through eigenvalue decomposition
  # construct pairwise contrast matrix
  # within-lineage between-condition DE
  # get linear predictor, condition specific
  combsPerCurve <- utils::combn(nConditions, 2)
  nComparisonsPerCurve <- ncol(combsPerCurve)
  ## construct generic contrast matrix to be filled in for all lineages
  Lsub <- matrix(0, nrow = nknots(models) * nConditions, 
                 ncol = nKnots * nComparisonsPerCurve)
  for (jj in seq_len(nComparisonsPerCurve)) {
    comp <- combsPerCurve[, jj]
    comp1ID <- ((comp[1] - 1) * nknots(models) + knots[1]):(
      (comp[1] - 1) * nknots(models) + knots[2])
    comp2ID <- ((comp[2] - 1) * nknots(models) + knots[1]):(
      (comp[2] - 1) * nknots(models) + knots[2])
    for (kk in seq_len(length(comp1ID))) {
      Lsub[c(comp1ID[kk], comp2ID[kk]), (jj - 1) * nKnots + kk] <- c(1, -1)
    }
  }
  # fill in contrast matrix for each lineage
  LWithin <- matrix(0, nrow = ncol(X), 
                    ncol = nLineages * nKnots * nComparisonsPerCurve)
  rownames(LWithin) <- colnames(X)
  # fill in contrast matrix
  # some helpers to identify coefficients for each lineage/condition
  smoothCoefs <- grep(x = colnames(X), pattern = "^s(t[1-9]+)*")
  smoothCoefNames <- colnames(X)[smoothCoefs]
  textAfterSplit <- unlist(lapply(strsplit(smoothCoefNames, split = ":l"), "[[", 2))
  # fill in with generic contrast
  for (jj in seq_len(nLineages)) {
    curvID <- substr(textAfterSplit, 1, 1) == jj
    limits <-  ((jj - 1) * nKnots * nComparisonsPerCurve + 1):(
      jj * nKnots * nComparisonsPerCurve)
    LWithin[smoothCoefs[curvID], limits] <- Lsub
  }
  colnames(LWithin) <- rep(paste0("lineage", seq_len(nLineages)),
                           each = nKnots * nComparisonsPerCurve)
  return(LWithin)
}


conditionTest <- function(models, global = TRUE, pairwise = FALSE, 
                          lineages = FALSE, knots = NULL, l2fc = 0, 
                          eigenThresh = 1e-2){
  # Ensure the models are ok for the conditionTest
  if (!(global | pairwise | lineages)) stop("One of global, pairwise or lineages must be true")
  checks <- find_conditions(models)
  conditions <- checks$conditions; nConditions <- checks$nConditions
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  if (is.null(knots)) {
    nKnots <- nknots(models)
    knots <- c(1, nKnots)
  } else {
    if (length(knots) != 2 | !(is.numeric(knots)) | any(knots > nknots(models))) {
      stop("knots must consists of 2 integers below the number of knots")
    }
    nKnots <- knots[2] - knots[1] + 1
  }
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData),
                                        pattern = "pseudotime")]
  # note that nCurves = nLineages * nConditions
  nCurves <- length(grep(x = colnames(dm), pattern = "l[(1-9)+]"))
  nLineages <- nCurves / nConditions
  if (nLineages == 1 & lineages == TRUE) {
    message("Only one lineage; skipping single-lineage comparison.")
    lineages <- FALSE
  }
  if (nConditions == 2 & pairwise == TRUE) {
    message("Only two conditions; skipping pairwise comparison.")
    pairwise <- FALSE
  }
  # get predictor matrix for every lineage.
  L <- construct_contrast_matrix_conditions(models, X, conditions, nConditions,
                                            nLineages, nKnots, knots)
  # perform Wald test and calculate p-value
  waldResultsOmnibus <- .allWaldStatGAMFC(models, L, l2fc, eigenThresh)
  
  # perform pairwise comparisons
  if (lineages & (!pairwise)) {
    # within-lineage DE between conditions
    # loop over list of within-lineage DE contrasts
    for (jj in seq_len(nLineages)) {
      LLin <- L[, colnames(L) == paste0("lineage", jj), drop = FALSE]
      waldResults <- .allWaldStatGAMFC(models, LLin, l2fc, eigenThresh)
      colnames(waldResults) <- c(
        paste0("waldStat_lineage", jj),
        paste0("df_lineage", jj),
        paste0("pvalue_lineage", jj)
      )
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldResAllLineages <- waldResults
      if (jj > 1) waldResAllLineages <- cbind(waldResAllLineages, waldResults)
    } 
  }
  
  if ((!lineages) & pairwise) {
    # all-lineage DE between pairs conditions
    # loop over list of pairs of conditions DE contrasts
    combs <- utils::combn(nConditions, 2)
    for (jj in seq_len(ncol(combs))) {
      conds <- combs[,jj]
      conds_L <- c(
        grep(x = rownames(L), pattern = paste0("_", conds[1], "\\.")),
        grep(x = rownames(L), pattern = paste0("_", conds[2], "\\."))
      )
      Lpair <- L
      Lpair[-conds_L, ] <- 0
      Lpair <- Lpair[, (colSums(Lpair) == 0) & (colSums(abs(Lpair)) != 0)]
      waldResults <- .allWaldStatGAMFC(models, Lpair, l2fc, eigenThresh)
      colnames(waldResults) <- c(
        paste0("waldStat_conds", paste0(conds, collapse = "vs")),
        paste0("df_conds", paste0(conds, collapse = "vs")),
        paste0("pvalue_conds", paste0(conds, collapse = "vs"))
      )
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldResAllPairs <- waldResults
      if (jj > 1) waldResAllPairs <- cbind(waldResAllPairs, waldResults)
    } 
  }
  
  if (lineages & pairwise) {
    for (ll in seq_len(nLineages)) {
      LLin <- L[, colnames(L) == paste0("lineage", ll), drop = FALSE]
      combs <- utils::combn(nConditions, 2)
      for (jj in seq_len(ncol(combs))) {
        conds <- combs[, jj]
        conds_L <- c(
          grep(x = rownames(LLin), pattern = paste0("_", conds[1], "\\.")),
          grep(x = rownames(LLin), pattern = paste0("_", conds[2], "\\."))
        )
        Lpair <- LLin; Lpair[-conds_L, ] <- 0
        Lpair <- Lpair[, (colSums(Lpair) == 0) & (colSums(abs(Lpair)) != 0)]
        waldResults <- .allWaldStatGAMFC(models, Lpair, l2fc, eigenThresh)
        colnames(waldResults) <- c(
          paste0("waldStat_lineage", ll, "_conds", paste0(conds, collapse = "vs")),
          paste0("df_lineage", ll, "_conds", paste0(conds, collapse = "vs")),
          paste0("pvalue_lineage", ll, "_conds", paste0(conds, collapse = "vs"))
        )
        waldResults <- as.data.frame(waldResults)
        if (jj == 1 & ll == 1) {
          waldResAll <- waldResults
        } else {
          waldResAll <- cbind(waldResAll, waldResults)
        }
      }
    }
  }
  
  # return output
  if ((global) & (!lineages) & (!pairwise)) return(waldResultsOmnibus)
  if ((!global) & (lineages) & (!pairwise)) return(waldResAllLineages)
  if ((!global) & (!lineages) & (pairwise)) return(waldResAllPairs)
  if ((global) & (lineages) & (!pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAllLineages))
  }
  if ((global) & (!lineages) & (pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAllPairs))
  }
  if ((!global) & (lineages) & (pairwise)) {
    return(cbind(waldResAll))
  }
  if ((global) & (lineages) & (pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAll))
  }
}


nknots <- function(models) {
  return(length(S4Vectors::metadata(models)$tradeSeq$knots))
}

.allWaldStatGAMFC <- function(models, L, l2fc, eigenThresh = 1e-2) {
  betaAll <- rowData(models)$tradeSeq$beta[[1]]
  sigmaAll <- rowData(models)$tradeSeq$Sigma
  waldRes <- lapply(seq_len(nrow(models)), function(ii){
    beta <- t(betaAll[ii,])
    Sigma <- sigmaAll[[ii]]
    if(any(is.na(beta)) | any(is.na(Sigma))) {
      return(c(NA, NA))
    } else {
      return(getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh))  
    }
  })
  names(waldRes) <- rownames(models)
  #tidy output
  waldResults <- do.call(rbind, waldRes)
  pval <- 1 - stats::pchisq(waldResults[, 1], df = waldResults[, 2])
  waldResults <- cbind(waldResults, pval)
  colnames(waldResults) <- c("waldStat", "df", "pvalue")
  waldRes <- as.data.frame(waldResults)
  return(waldRes)
}

getEigenStatGAMFC <- function(beta, Sigma, L, l2fc, eigenThresh = 1e-2){
  estFC <- t(L) %*% beta
  logFCCutoff <- log(2^l2fc) # log2 to log scale
  est <- sign(estFC)*pmax(0, abs(estFC) - logFCCutoff) # zero or remainder
  sigma <- t(L) %*% Sigma %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eSigma$values / eSigma$values[1] > eigenThresh), silent = TRUE)
  if (is(r,"try-error")) {
    return(c(NA, NA))
  }
  if (r == 1) return(c(NA, NA)) # CHECK
  halfCovInv <- eSigma$vectors[, seq_len(r)] %*% (diag(1 / sqrt(eSigma$values[seq_len(r)])))
  halfStat <- t(est) %*% halfCovInv
  stat <- crossprod(t(halfStat))
  return(c(stat, r))
}
