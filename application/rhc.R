library(effectbounds)

set.seed(10016)

data("RHC", package = "ATbounds")

X <- setdiff(colnames(RHC), c("survival", "RHC"))
A <- "RHC"
Y <- "survival"

thresholds <- 10^seq(-5, log10(0.05), 0.2)
smoothness <- c(0.001, 0.01)

bounds <- ate_bounds(
  RHC, X, A, Y, 
  smoothness = smoothness, 
  thresholds = thresholds
)

summary(bounds)

pdf("application/rhc_plot.pdf", width = 8, height = 6)
title <- "Population Average Treatment Effect of \nRight Heart Catheterization on 30-day Mortality"
ylim = c(-0.2, 0.10)
plot(bounds, 
     main = title,
     ylim = ylim,
     bounds_color = "darkred", smoothness = 0.01)
dev.off()

pdf("application/rhc_plot_sensitivity.pdf", width = 8, height = 6)
title <- "Population Average Treatment Effect of \nRight Heart Catheterization on 30-day Mortality"
ylim = c(-0.2, 0.10)
plot(bounds, 
     main = title,
     ylim = ylim,
     bounds_color = "darkred", smoothness = 0.001)
dev.off()
