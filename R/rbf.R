# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- kernlab::rbfdot(sigma = 1/l)
  return(sigf * kernlab::kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- kernlab::rbfdot(sigma = 1/l)
  return(sigf * kernlab::kernelMatrix(rbf, x = x, y = y))
}
