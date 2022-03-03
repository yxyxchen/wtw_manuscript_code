set.seed(1)
n <- 1e3
dat <- data.frame(
  x = 1:n,
  y = sin(seq(0, 5*pi, length.out = n)) + rnorm(n=n, mean = 0, sd=0.1)
)

approxData <- data.frame(
  with(dat, 
       approx(x, y, xout = seq(1, n, by = 10), method = "linear")
  ),
  method = "approx()"
)

splineData <- data.frame(
  with(dat, 
       spline(x, y, xout = seq(1, n, by = 10))
  ),
  method = "spline()"
)

smoothData <- data.frame(
  x = 1:n,
  y = as.vector(smooth(dat$y)),
  method = "smooth()"
)

loessData <- data.frame(
  x = 1:n,
  y = predict(loess(y~x, dat, span = 0.1)),
  method = "loess()"
)

library(ggplot2)
ggplot(rbind(approxData, splineData, smoothData, loessData), aes(x, y)) + 
  geom_point(dat = dat, aes(x, y), alpha = 0.2, col = "red") +
  geom_line(col = "blue") +
  facet_wrap(~method) +
  ggtitle("Interpolation and smoothing functions in R") +
  theme_bw(16)
