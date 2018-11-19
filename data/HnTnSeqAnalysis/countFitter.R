setup <- function(){
  library(MASS)
  library(pscl)
  library(ggplot2)
  library(car)
  setwd("~/Documents/berkeley/SavageLab/scripts")
  library(readr)
  hnea_genes <- read_csv("~/Documents/berkeley/SavageLab/scripts/hnea_genes.csv")
}

plot_nbinom <- function(x, siz, m) {
  return(dnbinom(round(x),size = siz, mu = m))
}
plot_Prob_nbinom <- function(x, siz, m) {
  return(pnbinom(round(x),size = siz, mu = m))
}

rootogram <- function(obj, max = NULL, ...) {
  y <- model.response(model.frame(obj))
  tab <- table(factor(y, levels = 0:max(y)))
  tab2 <- colSums(predprob(obj))
  if(is.null(max)) {
    max <- if(all(tab2 >= 1)) max(y) else max(ceiling(mean(y)), min(which(tab2 < 1)) - 1)
  }
  max <- min(max, length(tab) - 1) + 1
  obsrvd <- sqrt(tab[1:max])
  expctd <- sqrt(tab2[1:max])
  res <- obsrvd - expctd
  x <- barplot(obsrvd, offset = -res, xlab = "Count", ylab = "sqrt(Frequency)")
  lines(x, expctd, col = 2, type = ifelse(max > 25, "l", "b"), pch = 19)
  abline(h = 0)
  invisible(cbind(observed = tab, expected = tab2))
}

plotZINB <- function(zinb){
  x<- zinb$y
  p <- ggplot() +
    geom_histogram(aes(x, ..density..), binwidth = 5, colour = "black", 
                   fill = "white") +
    stat_function(geom = "line", data=data.frame(zinb$y), fun = plot_nbinom,
                  args = c(zinb$theta,unname(exp(coef(zinb)[1]))),
                  colour = "red", lwd = 1.5) +
    ylab("Density")
  print(p)
  rootogram(zinb)
}

findEssentialsZinb <- function(data){
  dat<-data.frame(data)
  colnames(dat) <- "X1"
  zinb<-zeroinfl(formula = X1 ~ 1 | 1, data = dat, dist = "negbin", 
                 y = TRUE)
  plotZINB(zinb)
  return(zinb)
}

qqAndMeanDiff <- function(values, fitt){
  temp <- cbind(values,fitt)
  forMeanDiff <- data.frame(temp)
  forMeanDiff$mean <- (forMeanDiff$X1+forMeanDiff$X2)/2
  forMeanDiff$diff <- (forMeanDiff$X1-forMeanDiff$X2)
  MDReg <- lm(formula = mean~diff, data = forMeanDiff)
  print("Mean Diff regression Summary:")
  print(summary(MDReg))
  aveDiff <-mean(forMeanDiff$diff)
  down95Diff <- mean(forMeanDiff$diff) - (1.96 * sd(forMeanDiff$diff))
  up95Diff <- mean(forMeanDiff$diff) + (1.96 * sd(forMeanDiff$diff))
  p <- ggplot(forMeanDiff, aes(x = mean, y = diff)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = aveDiff, colour = "blue", size = 0.5) +
    geom_hline(yintercept = down95Diff, colour = "red", size = 0.5) +
    geom_hline(yintercept = up95Diff, colour = "red", size = 0.5) +
    ylab("Diff. Between Measures") +
    xlab("Average Measure") +
    abline(MDReg, col = "green")
  
  print(p)
  
  QQReg <- lm(formula = X1~X2, data = forMeanDiff)
  print("QQ regression Summary:")
  print(summary(QQReg))
  qqplot(forMeanDiff$X1,forMeanDiff$X2)
  abline(QQReg, col = "green")
  
}

distTest <- function(data){
  dat <-data+1
  qqp(dat, "norm")
  qqp(dat, "lnorm")
  nbinom <- fitdistr(dat, "Negative Binomial")
  qqp(dat, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  poisson <- fitdistr(dat, "Poisson")
  qqp(dat, "pois", poisson$estimate)
  gamma <- fitdistr(dat, "gamma")
  qqp(dat, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
}
# 
# plot_mix_comps <- function(x, alpha, beta, lambda) {
#   return(lambda * dgamma(x, alpha, shape = beta))
# }
# 
# plotModel <- function(mixmdl){
#   x<- data.frame(mixmdl$x)
#   x1 <- mixmdl$x
#   p <- ggplot() +
#     geom_histogram(aes(x1, ..density..), binwidth = 0.01, colour = "black", 
#                    fill = "white") +
#     stat_function(geom = "line", data=x, fun = plot_mix_comps,
#                   args = list(mixmdl$gamma.pars[1,1], mixmdl$gamma.pars[2,1], lam = mixmdl$lambda[1]),
#                   colour = "red", lwd = 1.5) +
#     stat_function(geom = "line", data=x, fun = plot_mix_comps,
#                   args = list(mixmdl$gamma.pars[1,2], mixmdl$gamma.pars[2,2], lam = mixmdl$lambda[2]),
#                   colour = "blue", lwd = 1) +
#     ylab("Density")
#   print(p)
# }
# 
# plotProbs<- function(model){
#   hist(model$x,breaks = seq(1,1.3,by=0.01))
#   points(model$x,model$posterior[,2]*500,col='red')
#   points(model$x,model$posterior[,1]*500)
# }
# 
# findEssentials <- function(data){
#   if(min(data)>0){
#     dat <- data
#   } else{
#     dat <- data+1
#   }
#   mixmdl <- gammamixEM(dat,maxit = 10000,k=2)
#   plotProbs(mixmdl)
#   plotModel(mixmdl)
#   return(mixmdl)
# }
# 
# findSplit <- function(data){
#   tt<- density(data)
#   plot(tt)
#   s<-data.frame()
#   point <- (which(diff(sign(diff(tt$y)))==2)+1)[1]
#   points(tt$x[point],tt$y[point],col='red')
#   small <- data[data<=tt$x[point]]
#   big <- data[data>tt$x[point]]
#   return(list("point" = point, "big" = big, "small" = small))
# }
# 
# findEssentialSeeded <- function(data){
#   #if(min(data)>0){
#   #  dat <- data
#   #} else{
#   #  dat <- data+1
#   #}
#   dat<-data+1
#   split<-findSplit(dat)
#   normSmall <- fitdistr(split$small,'gamma')
#   normBig <- fitdistr(split$big,'gamma')
#   mixmdl <- gammamixEM(dat,alpha = c(normSmall$estimate[1],normBig$estimate[1]),beta = c(1/normSmall$estimate[2],1/normBig$estimate[2]),maxit = 10000,k=2)
#   plotProbs(mixmdl)
#   plotModel(mixmdl)
#   return(mixmdl)
# }
# plot_norm_mix_comps <- function(x, mu, sigma, lam) {
#   lam * dnorm(x, mu, sigma)
# }
# 
# plotNormModel <- function(mixmdl){
#   x<- data.frame(mixmdl$x)
#   x1 <- mixmdl$x
#   p <- ggplot() +
#     geom_histogram(aes(x1, ..density..), binwidth = 0.01, colour = "black", 
#                    fill = "white") +
#     stat_function(geom = "line", data=x, fun = plot_norm_mix_comps,
#                   args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
#                   colour = "red", lwd = 1.5) +
#     stat_function(geom = "line", data=x, fun = plot_norm_mix_comps,
#                   args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
#                   colour = "blue", lwd = 1) +
#     ylab("Density")
#   print(p)
# }
# 
# findEssentialsNorm <- function(data){
#   dat <- data+1
#   split<-findSplit(dat)
#   normSmall <- fitdistr(split$small,'normal')
#   normBig <- fitdistr(split$big,'normal')
#   mixmdl <- normalmixEM(dat,mu = c(normSmall$estimate[1],normBig$estimate[1]),sigma = c(normSmall$estimate[2],normBig$estimate[2]),maxit = 10000,k=2)
#   #plotProbs(mixmdl)
#   plotNormModel(mixmdl)
#   return(mixmdl)
# }

