predict.GP <- function(gp, x) {

  if(gp$ker == 'rbf') source("R/rbf.R")
  if(gp$ker == 'laplace') source("R/laplace.R")

  testx <- x
  x <- x/gp$xmax
  x <- as.matrix(x)
  theta <- gp$theta
  n <- dim(gp$x)[1]
  nx <- dim(x)[1]
  #x <- t((t(as.matrix(x)) - gp$mx) / gp$sdx)

   kx <- ker2(x = x, y = gp$tx, l = theta[1], sigf = theta[2])
   kxx <- ker(x = x, l = theta[1], sigf = theta[2]) + theta[3] * diag(nx)
  # k <- ker(x = gp$x, l = theta[1], sigf = theta[2]) + theta[3] * diag(n)
  # kinv <- chol2inv(chol(k))
  mu <- kx %*% (gp$kinv %*% as.matrix(gp$ty))
  mu <- mu * gp$sdy + gp$my
  sigma <- kxx - kx %*% (gp$kinv %*% t(kx))
  diag(sigma)[diag(sigma)<0] <- 0
  ll <- mu - 1.96 * (sqrt(diag(sigma)) * gp$sdy)
  ul <- mu + 1.96 * (sqrt(diag(sigma)) * gp$sdy)

  pred <- list(mu = mu, sigma = sigma, ll = ll, ul = ul, x = gp$x, y = gp$y, testx = testx)
  class(pred) <- 'prediction_GP'
  return(pred)

}




plot.prediction_GP <- function(fitted,...) {

  ggdata <- data.frame(fitted$testx, fitted$mu, fitted$ll, fitted$ul)
  ggplot(data = ggdata, aes(x = fitted$testx, y = fitted$mu)) + geom_line() +
    geom_smooth(aes(ymin = fitted$ll, ymax = fitted$ul), stat = 'identity', size = 2) +
    xlab('X') + ylab('Y')

}
