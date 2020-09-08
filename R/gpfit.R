gpfit <- function(x, y, ker = 'rbf',...) {

  if(ker == 'rbf') source("R/rbf.R")
  if(ker == 'laplace') source("R/laplace.R")

  xmax <- max(x)
  tx <- x/xmax
  tx <- as.matrix(tx)
  n <- dim(tx)[1]

  #mx <- apply(x, 2, mean)
  #sdx <- apply(x, 2, sd)
  my <- mean(y)
  sdy <- sd(y)

  #x <- t((t(as.matrix(x)) - mx) / sdx)
  ty <- (y- my)/sdy

  theta <- Hyper(tx, ty)

  k <- ker(x = tx, l = theta[1], sigf = theta[2]) + theta[3] * diag(n)
  kinv <- chol2inv(chol(k))

  gpfit <- list(x = x, y=y, theta = theta, ker = ker, my = my, sdy = sdy, xmax = xmax, kinv = kinv,tx = tx, ty = ty)

  class(gpfit) <- 'GP'
  return(gpfit)
}




plot.GP <- function(gp,...) {
  fitted <- predict(gp, gp$x)

  ggdata <- data.frame(fit$x, fit$y, fitted$mu, fitted$ll, fitted$ul)
  ggplot(data = ggdata, aes(x = fit$x, y = fitted$mu)) + geom_line() +
    geom_point(aes(y = fit$y), size = 2) +
    geom_smooth(aes(ymin = fitted$ll, ymax = fitted$ul), stat = 'identity', size = 2) +
    xlab('X') + ylab('Y')
}



