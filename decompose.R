set.seed(12345)
ymat <- matrix(rnorm(12), 4, 3)
rmat <- crossprod(matrix(rnorm(40), 10, 4)) / 10
cmat <- crossprod(matrix(rnorm(30), 10, 3)) / 10

decompose <- function(xmat, rmat, cmat) {
  n <- nrow(xmat)
  m <- ncol(xmat)
  pmat <- matrix(colSums(rmat), n, n, byrow = TRUE) / sum(rmat)
  qmat <- matrix(rowSums(cmat), m, m) / sum(cmat)
  px <- pmat %*% xmat
  xq <- xmat %*% qmat
  dd <- rep(list(0),4) 
  dd[[1]] <- px %*% qmat
  dd[[2]] <- px - dd[[1]]
  dd[[3]] <- xq - dd[[1]]
  dd[[4]] <- xmat - dd[[1]] - dd[[2]] - dd[[3]]
  ip <- matrix(0, 4, 4)
  for (i in 1:4) {
    for (j in 1:4) {
      ip[i, j] <- kroneckerIP(dd[[i]], dd[[j]], rmat, cmat)
    }
  }
  ssq <- c(kroneckerIP(xmat, xmat, rmat, cmat), sum(diag(ip)))
  return(list(dd = dd, ip = ip, ssq = ssq))
}

kroneckerIP <- function(xmat, ymat, rmat, cmat) {
  return(sum(ymat * (rmat %*% xmat %*% cmat)))
}