#theo.gcca, to use kernels when needed

library(GUniFrac)
library(RGCCA)

theo.gcca <- function(blocks,
                      connection = 1 - diag(length(blocks)),
                      useKernel = rep(0, length(blocks)),
                      kernel = rep(0, length(blocks)),
                      scale = TRUE,
                      tau = rep(1, length(blocks)),
                      ncomp = rep(1, length(blocks)),
                      scheme = "centroid",
                      init = "svd",
                      bias = TRUE,
                      tol = 1e-16,
                      verbose = TRUE,
                      na.rm = TRUE,
                      quiet = FALSE,
                      superblock = FALSE)

{
  if (mode(scheme) != "function") {
    if (verbose) {
      cat("Computation of the RGCCA block components based on the",
          scheme,
          "scheme \n")
    }
  }
  if (mode(scheme) == "function" & verbose) {
    cat("Computation of the RGCCA block components based on the g scheme \n")
  }
  
  if (!is.numeric(tau) & verbose) {
    cat("Optimal Shrinkage intensity parameters are estimated \n")
  } else {
    if (is.numeric(tau) & verbose) {
      cat("Shrinkage intensity parameters are chosen manually \n")
    }
  }
  
  if (scale == TRUE) {
    blocks = lapply(blocks, scale)
    # blocks = blocks / sqrt(blocks)
    # kernel = lapply(kernel, scale)
    }
  
  # ndefl number of deflation per block
  ndefl <- ncomp - 1
  N <- max(ndefl)
  J <- length(blocks)
  pjs <- sapply(blocks, NCOL)
  nb_ind <- NROW(blocks[[1]])
  AVE_X <- list()
  AVE_outer <- rep(NA, max(ncomp))
  
  Y <- NULL
  P <- a <- astar <- list()
  crit <- list()
  AVE_inner <- rep(NA, max(ncomp))
  
  # Whether primal or dual
  primal_dual <- rep("primal", J)
  primal_dual[which(nb_ind < pjs)] <- "dual"
  primal_dual[which(useKernel == TRUE)] <- "dual"
  
  for (b in seq_len(J)) {
    if (primal_dual[b] == "dual"){
      a[[b]] <- matrix(NA, nrow(blocks[[b]]), N + 1)
      Y[[b]] <- matrix(NA, nb_ind, N + 1)
    }
    else {
      a[[b]] <- matrix(NA, pjs[[b]], N + 1)
      Y[[b]] <- matrix(NA, nb_ind, N + 1)      
    }
  }
  
  if (!superblock) {
    for (b in seq_len(J))
      astar[[b]] <- matrix(NA, pjs[b], N + 1)
  } else {
    astar <- matrix(NA, pjs[J], N + 1)
  }
    
  # Save computed shrinkage parameter in a new variable
  computed_tau <- tau
  if (is.vector(tau))
    computed_tau <- matrix(NA, nrow = N + 1, J)
  
  # First component block
  if (is.vector(tau)) {
    rgcca_result <- theogccak(
      blocks,
      connection,
      useKernel = useKernel,
      kernel = kernel,
      primal_dual = primal_dual,
      tau = tau,
      scheme = scheme,
      init = init,
      bias = bias,
      tol = tol,
      verbose = verbose,
      na.rm = na.rm
    )
  } else {
    rgcca_result <- theogccak(
      blocks,
      connection,
      useKernel = useKernel,
      kernel = kernel,
      primal_dual = primal_dual,
      tau = tau[1,],
      scheme = scheme,
      init = init,
      bias = bias,
      tol = tol,
      verbose = verbose,
      na.rm = na.rm
    )
  }
  computed_tau[1,] <- rgcca_result$tau
  
  for (b in seq_len(J))
    Y[[b]][, 1] <- rgcca_result$Y[, b, drop = FALSE]
  for (b in seq_len(J))
    a[[b]][, 1] <- rgcca_result$a[[b]]
  
  ifelse(!superblock,
         astar <- a,
         astar[, 1] <- a[[J]][, 1, drop = FALSE])
  
  AVE_inner[1] <- rgcca_result$AVE_inner
  crit[[1]] <- rgcca_result$crit
  
  if (N > 0) {
    R <- blocks
    
    if (!superblock) {
      for (b in seq_len(J))
        P[[b]] <- matrix(NA, pjs[b], N)
    } else {
      P <- matrix(NA, pjs[J], N)
    }
    
    for (n in 2:(N + 1)) {
      if (verbose) {
        cat(
          paste0(
            "Computation of the RGCCA block components #",
            n,
            " is under
                 progress...\n"
          )
        )
      }
      
      if (!superblock) {
        defl_result <- defl_select(rgcca_result$Y, R,
                                   ndefl, n - 1, J,
                                   na.rm = na.rm)
        R <- defl_result$resdefl
        for (b in seq_len(J))
          P[[b]][, n - 1] <- defl_result$pdefl[[b]]
      } else {
        defl_result <- deflation(R[[J]], rgcca_result$Y[, J])
        R[[J]] <- defl_result$R
        P[, n - 1] <- defl_result$p
        cumsum_pjs <- cumsum(pjs)[seq_len(J - 1)]
        inf_pjs <- c(0, cumsum_pjs[seq_len(J - 2)]) + 1
        for (j in seq_len(J - 1)) {
          R[[j]] <- R[[J]][, inf_pjs[j]:cumsum_pjs[j], drop = FALSE]
          rownames(R[[j]]) <- rownames(R[[j]])
          colnames(R[[j]]) <-
            colnames(R[[J]])[inf_pjs[j]:cumsum_pjs[j]]
        }
      }
      
      if (is.vector(tau)) {
        rgcca_result <- theogccak(
          R,
          connection,
          useKernel = useKernel,
          kernel = kernel,
          primal_dual = primal_dual,
          tau = tau,
          scheme = scheme,
          init = init,
          bias = bias,
          tol = tol,
          verbose = verbose,
          na.rm = na.rm
        )
      } else {
        rgcca_result <- theogccak(
          R,
          connection,
          useKernel = useKernel,
          kernel = kernel,
          primal_dual = primal_dual,
          tau = tau[n,],
          scheme = scheme,
          init = init,
          bias = bias,
          tol = tol,
          verbose = verbose,
          na.rm = na.rm
        )
      }
      
      computed_tau[n,] <- rgcca_result$tau
      
      AVE_inner[n] <- rgcca_result$AVE_inner
      crit[[n]] <- rgcca_result$crit
      
      for (b in seq_len(J))
        Y[[b]][, n] <- rgcca_result$Y[, b]
      for (b in seq_len(J))
        a[[b]][, n] <- rgcca_result$a[[b]]
      
      if (!superblock) {
        for (b in seq_len(J)) {
          astar[[b]][, n] <- rgcca_result$a[[b]] -
            astar[[b]][, (1:(n - 1)), drop = F] %*%
            drop(t(a[[b]][, n]) %*% P[[b]][, 1:(n - 1), drop = F])
        }
      } else {
        astar[, n] <- rgcca_result$a[[J]] -
          astar[, (1:(n - 1)), drop = F] %*%
          drop(t(a[[J]][, n]) %*% P[, 1:(n - 1), drop = F])
      }
    }
  }
  
  for (j in seq_len(J)) {
    if (useKernel[j]){
      AVE_X[[j]] <- apply(cor(kernel[[j]], Y[[j]], use = "pairwise.complete.obs") ^
                            2, 2, mean)
    }
    else{
      AVE_X[[j]] <- apply(cor(blocks[[j]], Y[[j]], use = "pairwise.complete.obs") ^
                            2, 2, mean)
    }
  }
  
  outer <- matrix(unlist(AVE_X), nrow = max(ncomp))
  
  for (j in seq_len(max(ncomp))) {
    AVE_outer[j] <- sum(pjs * outer[j,]) / sum(pjs)
  }
  
  AVE_X <- shave(AVE_X, ncomp)
  
  AVE <-
    list(AVE_X = AVE_X,
         AVE_outer = AVE_outer,
         AVE_inner = AVE_inner)
  
  if (N == 0) {
    crit <- unlist(crit)
    computed_tau <- as.vector(computed_tau)
  }
  
  out <- list(
    Y = Y,
    a = a,
    astar = astar,
    tau = computed_tau,
    crit = crit,
    primal_dual = primal_dual,
    AVE = AVE
  )
  
  class(out) <- "rgccad"
  
  return(out)
}








theogccak <-
  function(A,
           C,
           useKernel = rep(0, length(A)),
           kernel = rep(0, length(A)),
           primal_dual = primal_dual,
           tau = rep(1, length(A)),
           scheme = "centroid",
           verbose = FALSE,
           init = "svd",
           bias = TRUE,
           tol = 1e-08,
           na.rm = TRUE)
    {
    if (mode(scheme) != "function") {
      if (scheme == "horst") {
        g <- function(x)
          x
        ctrl <- FALSE
      }
      if (scheme == "factorial") {
        g <- function(x)
          x ^ 2
        ctrl <- TRUE
      }
      if (scheme == "centroid") {
        g <- function(x)
          abs(x)
        ctrl <- TRUE
      }
    } else {
      # check for parity of g
      g <- scheme
      ctrl <- !any(g(-5:5) != g(5:-5))
    }
    
    dg <- Deriv::Deriv(g, env = parent.frame())
    
    J <- length(A) # number of blocks
    n <- NROW(A[[1]]) # number of individuals
    pjs <- sapply(A, NCOL) # number of variables per block
    Y <- Z <- matrix(0, n, J)
    
    if (!is.numeric(tau)) {
      # From Schafer and Strimmer, 2005
      tau <- sapply(A, tau.estimate, na.rm = na.rm)
    }
    
    A <- lapply(A, as.matrix)
    a <- alpha <- M <- Minv <- K <- list()
    
    # Test for primal or dual for each block
    which.primal <- which(primal_dual == "primal")
    which.dual <- which(primal_dual == "dual")
    
    # Initialisation by SVD
    if (init == "svd") {
      for (j in which.primal) {
        a[[j]] <- initsvd(A[[j]])
      }
      for (j in which.dual) {
        if (useKernel[j]){
          alpha[[j]] <- initsvd(kernel[[j]])
          K[[j]] <- kernel[[j]]
        }
        else {
          alpha[[j]] <- initsvd(A[[j]])
          K[[j]] <- pm(A[[j]], t(A[[j]]), na.rm = na.rm)
        }
      }
    } else if (init == "random") {
      for (j in which.primal) {
        a[[j]] <- rnorm(pjs[j]) # random initialisation
      }
      
      for (j in which.dual) {
        if (useKernel[j]){
          alpha[[j]] <- initsvd(kernel[[j]])
          K[[j]] <- kernel[[j]]
        }
        else {
          alpha[[j]] <- initsvd(A[[j]])
          K[[j]] <- pm(A[[j]], t(A[[j]]), na.rm = na.rm)
        }
      }
    } else {
      stop_rgcca("init should be either random or by SVD.")
    }
    
    N <- ifelse(bias, n, n - 1)
    for (j in which.primal) {
      ifelse(tau[j] == 1,
             yes = {
               a[[j]] <- drop(1 / sqrt(t(a[[j]]) %*% a[[j]])) * a[[j]]
               Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
             },
             no = {
               M[[j]] <- ginv(tau[j] * diag(pjs[j]) + ((1 - tau[j])) * 1 / N *
                                (pm(t(A[[j]]), A[[j]], na.rm = na.rm)))
               a[[j]] <-
                 drop(1 / sqrt(t(a[[j]]) %*% M[[j]] %*% a[[j]])) *
                 (M[[j]] %*% a[[j]])
               Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
             })
    }
    for (j in which.dual) {
      ifelse(tau[j] == 1,
             yes = {
               alpha[[j]] <- drop(1 / sqrt(t(alpha[[j]]) %*% K[[j]] %*%
                                             alpha[[j]])) * alpha[[j]]
               a[[j]] <- alpha[[j]]
               Y[, j] <- pm(K[[j]], a[[j]], na.rm = na.rm)
             },
             no = {
               M[[j]] <- tau[j] * diag(n) + ((1 - tau[j])) * 1 / N * K[[j]]
               Minv[[j]] <- ginv(M[[j]])
               alpha[[j]] <- drop(1 / sqrt(t(alpha[[j]]) %*%
                                             M[[j]] %*% K[[j]] %*% alpha[[j]])) * alpha[[j]]
               a[[j]] <- alpha[[j]] # BEWARE THAT a[[j]] ARE ACTUALLY alpha[[j]]
               Y[, j] <- pm(K[[j]], a[[j]], na.rm = na.rm)
             })
    }
    
    iter <- 1
    n_iter_max <- 1000L
    crit <- numeric(n_iter_max)
    crit_old <- sum(C * g(cov2(Y, bias = bias)))
    a_old <- a
    
    repeat {
      for (j in which.primal) {
        dgx <- dg(cov2(Y[, j], Y, bias = bias))
        CbyCovj <- drop(C[j,] * dgx)
        if (tau[j] == 1) {
          Z[, j] <- Y %*% CbyCovj
          Az <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
          a[[j]] <- drop(1 / sqrt(crossprod(Az))) * Az
          Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
        } else {
          Z[, j] <- Y %*% CbyCovj
          Az <- pm(t(A[[j]]), Z[, j], na.rm = TRUE)
          a[[j]] <-
            drop(1 / sqrt(t(Az) %*% M[[j]] %*% Az)) * (M[[j]] %*% Az)
          Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
        }
      }
      
      for (j in which.dual) {
        dgx <- dg(cov2(Y[, j], Y, bias = bias))
        CbyCovj <- drop(C[j,] * dgx)
        ifelse(tau[j] == 1,
               yes = {
                 Z[, j] <- Y %*% CbyCovj
                 alpha[[j]] <-
                   drop(1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Z[, j])) * Z[, j]
                 a[[j]] <- alpha[[j]]
                 Y[, j] <- pm(K[[j]], a[[j]], na.rm = na.rm)
               },
               no = {
                 Z[, j] <- Y %*% CbyCovj
                 alpha[[j]] <- drop(1 / sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*% Z[, j])) * (Minv[[j]] %*% Z[, j])
                 
                 a[[j]] <- alpha[[j]]
                 Y[, j] <- pm(K[[j]], a[[j]], na.rm = na.rm)
               })
      }
      
      crit[iter] <- sum(C * g(cov2(Y, bias = bias)))
      if (verbose & (iter %% 1) == 0) {
        cat(
          " Iter: ",
          formatC(iter, width = 3, format = "d"),
          " Fit:",
          formatC(
            crit[iter],
            digits = 8,
            width = 10,
            format = "f"
          ),
          " Dif: ",
          formatC(
            crit[iter] - crit_old,
            digits = 8,
            width = 10,
            format = "f"
          ),
          "\n"
        )
      }
      
      stopping_criteria <- c(drop(crossprod(unlist(a, F, F) - unlist(a_old, F, F))),
                             abs(crit[iter] - crit_old))
      
      if (any(stopping_criteria < tol) | (iter > 1000)) {
        break
      }
      crit_old <- crit[iter]
      a_old <- a
      iter <- iter + 1
    }
    
    for (j in seq_len(J)) {
      if (ctrl & a[[j]][1] < 0) {
        a[[j]] <- -a[[j]]
        if (primal_dual[j] == "dual") {
          Y[, j] <- pm(K[[j]], a[[j]], na.rm = na.rm)
        }
        else {
          Y[, j] <- pm(A[[j]], a[[j]], na.rm = na.rm)
        }
      }
    }
    
    crit <- crit[which(crit != 0)]
    
    if (iter > n_iter_max) {
      warning("The RGCCA algorithm did not converge after",
              n_iter_max,
              " iterations.")
    }
    if (iter < n_iter_max & verbose) {
      cat("The RGCCA algorithm converged to a stationary point after",
          iter - 1,
          "iterations \n")
    }
    if (verbose) {
      plot(crit, xlab = "iteration", ylab = "criteria")
    }
    
    AVEinner <- sum(C * cor(Y) ^ 2 / 2) / (sum(C) / 2)
    result <- list(
      Y = Y,
      a = a,
      crit = crit,
      AVE_inner = AVEinner,
      tau = tau
    )
    return(result)
  }
