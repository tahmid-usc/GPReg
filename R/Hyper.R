#require(mvtnorm)
#require(kernlab)

Hyper <- function(x, y,...) {

  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    - mvtnorm::dmvnorm(x = y, sigma = k + theta[3] * diag(n), log = T)
  }

  hyp <- optim(par=rep(1, 3), fn = marlik, method = 'BFGS',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par^2)

}
