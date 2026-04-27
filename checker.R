
set.seed(12345)
amat <- array(rnorm(60), c(4, 3, 5))
rmat <- crossprod(matrix(rnorm(40), 10, 4)) / 10
cmat <- crossprod(matrix(rnorm(30), 10, 3)) / 10
smat <- crossprod(matrix(rnorm(50), 10, 5)) / 10

oo <- function(am = amat,
               rm = rmat,
               cm = cmat,
               sm = smat) {
  avec <- as.vector(am)
  adim <- dim(am)
  n <- adim[1]
  m <- adim[2]
  l <- adim[3]
  kp <-  kronecker(sm, kronecker(cm, rm))
  s1 <- sum(avec * (kp %*% avec))
  s2 <- 0.0
  for (i in 1:n) {
    for (j in 1:m) {
      for (k in 1:l) {
        for (ii in 1:n) {
          for (jj in 1:m) {
            for (kk in 1:l) {
              s2 <- s2 + (rmat[i, ii] * cmat[j, jj] *
                smat[k, kk] * amat[i, j, k] * amat[ii, jj, kk])
            }
          }
        }
      }
    }
  }
  return(c(s1, s2))
}

xmat <- matrix(rnorm(12), 4, 3)

v <- matrix(rnorm(6),3,2)
u <- matrix(rnorm(8),4,2)
v[, 1] <- v[, 1] / sqrt(sum(v[, 1] * (cmat %*% v[, 1])))
v[, 2] <- v[, 2] - sum(v[, 2] * (cmat %*% v[, 1])) * v[, 1]
v[, 2] <- v[, 2] / sqrt(sum(v[, 2] * (cmat %*% v[, 2])))
u[, 1] <- u[, 1] / sqrt(sum(u[, 1] * (rmat %*% u[, 1])))
u[, 2] <- u[, 2] - sum(u[, 2] * (rmat %*% u[, 1])) * u[, 1]
u[, 2] <- u[, 2] / sqrt(sum(u[, 2] * (rmat %*% u[, 2])))
p <- u %*% crossprod(u, rmat)
q <- tcrossprod(cmat %*% v, v)
rxc <- rmat %*% xmat %*% cmat
pm <- p %*% xmat %*% q
pa <- p %*% xmat - pm
pb <- xmat %*% q - pm
ph <- xmat - pa - pb - pm

loss <- function(x, y, rm = rmat, cm = cmat) {
  return(sum(y * (rm %*% x %*% cm)))
}