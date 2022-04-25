# A set of useful fucnctions

make.symmetric = function(a, lower.tri = TRUE) {
  if (lower.tri) {
    ind <- upper.tri(a)
    a[ind] <- t(a)[ind]
  } else {
    ind <- lower.tri(a)
    a[ind] <- t(a)[ind]
  }
  a
}


dist2kern.exp = function(distmatrix) {
  N = dim(distmatrix)[1]
  I_N = diag(N)
  un_N = replicate(N, 1)
  D2 = distmatrix * distmatrix
  K_1 = -0.5 * (I_N - (un_N %*% t(un_N)) / N) %*% D2 %*% (I_N - (un_N %*% t(un_N)) /
                                                            N)
  K_1 = make.symmetric(K_1)
  colnames(K_1) = rownames(K_1) = colnames(distmatrix)
  return(exp(K_1 * 4))
}


dist2kern.exp2 = function(distmatrix, number) {
  N = dim(as.matrix(distmatrix))[1]
  I_N = diag(N)
  un_N = replicate(N, 1)
  D2 = distmatrix * distmatrix
  K_1 = -0.5 * (I_N - (un_N %*% t(un_N)) / N) %*% D2 %*% (I_N - (un_N %*% t(un_N)) /
                                                            N)
  K_1 = make.symmetric(K_1)
  colnames(K_1) = rownames(K_1) = colnames(distmatrix)
  return(exp(K_1 * number))
}


dist2kern.trunc = function(distmatrix) {
  N = dim(as.matrix(distmatrix))[1]
  I_N = diag(N)
  un_N = replicate(N, 1)
  D2 = distmatrix * distmatrix
  K_1 = -0.5 * (I_N - (un_N %*% t(un_N)) / N) %*% D2 %*% (I_N - (un_N %*% t(un_N)) /
                                                            N)
  K_1 = make.symmetric(K_1)
  eK <- eigen(K_1, symmetric = TRUE)
  K_1 <- eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  colnames(K_1) = rownames(K_1) = colnames(distmatrix)
  return(K_1)
}


center.train = function(K_train) {
  N = dim(as.matrix(K_train))[1]
  I_N = diag(N)
  un_N = replicate(N, 1)
  K_center = (I_N - (un_N %*% t(un_N)) / N) %*% K_train %*% (I_N - (un_N %*% t(un_N)) /
                                                               N)
  K_center = make.symmetric(K_center)
  eK <- eigen(K_center, symmetric = TRUE)
  K_center <- eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  return(K_center)
}


center.test = function(K_test, K_train) {
  Nt = dim(as.matrix(K_test))[1]
  N = dim(as.matrix(K_train))[2]
  
  I_N = diag(N)
  un_Nt = as.matrix(replicate(Nt, 1))
  un_N = as.matrix(replicate(N, 1))
  K_center = (K_test - ((un_Nt %*% t(un_N)) / N) %*% K_train) %*% (I_N - (un_N %*% t(un_N)) /
                                                                     N)
  
  K_center = make.symmetric(K_center)
  eK <- eigen(K_center, symmetric = TRUE)
  K_center <- eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  return(K_center)
}