#/*
# * Dependencies and functions for N_invest.R analysis
# *
# *
# */
######################
# AUXILLIARY FUNCTIONS
######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figs', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# REQUIRED LIBRARIES
####################
library(rstan)
library(loo)
library(corrplot)
library(lubridate)
library(MASS)
library(coda)
library(MCMCpack)
library(MCMCglmm)
library(lme4) # remember to detatch("package:nlme") because of conflicts
#library(boot)
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(bayesboot)

################
# STAN OPTIONS
################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################
# SUMMARIZE MCMC RESULTS
########################
MCMCsum <- function(x) {
   data.frame(Mean=mean(x, na.rm=TRUE) , Median=median(x, na.rm=TRUE), t(quantile(x,probs=c(0,0.2,0.25,0.5,0.75,0.8,1),na.rm=TRUE)), HPDinterval(as.mcmc(x)))
}

aMCMCsum <- function(x) {
   c(mean(x, na.rm=TRUE) , median(x, na.rm=TRUE), t(quantile(x, probs=c(0,0.2,0.25,0.5,0.75,0.8,1),na.rm=TRUE)), HPDinterval(as.mcmc(x))[1:2])
}


looDiffSE  <-  function(x) {
        n  <-  length(x)
        sqrt(n * var(x))
}


#######################
# SUMMARIZE LOO RESULTS
#######################
#' Summarize loo model comparison results
#'
#' @title Make a table to summarize LOO model comparisons.
#' @param looDiff Result table from baked in function compare() from the 'loo' package.
#' @return A Results table in matrix form with pairwise elpd_loo differences, s.e., and p-values.
#' @author Colin Olito.
#' @export
makeLooTable <- function(looDiff) {
    # Variables, Containers
    index            <-  combn(nrow(looDiff),2)
    if(ncol(index) > 10)
      index  <-  index[,1:(max(index)-1)]
    nameIndex        <-  row.names(looDiff)
    elpd_pair        <-  rep(0, length = ncol(index))
    selooDiff        <-  rep(0, length = ncol(index))
    bayesBootDiff    <-  rep(0, length = ncol(index))
    bayesBootSE      <-  rep(0, length = ncol(index))
    bayesBootSE95    <-  rep(0, length = ncol(index))
    pDiff            <-  rep(0, length = ncol(index))
    pDiffBootSE      <-  rep(0, length = ncol(index))
    pDiffBootSE95    <-  rep(0, length = ncol(index))
    names            <-  rep(0, length = ncol(index))
    looList          <-  lapply(row.names(looDiff), get)
    
    # Check sample sizes
    for (i in 1:ncol(index)) {
        if(length(looList[[index[1,i]]]$pointwise[,"elpd_loo"]) != 
            length(looList[[index[2,i]]]$pointwise[,"elpd_loo"]))
            stop("sample sizes (chain lengths) differ between models")
    }
    n  <-  length(looList[[index[1,1]]]$pointwise[,"elpd_loo"])

   # Pairwise elpd_loo differences for all models
    for (i in 1:ncol(index)) {
        elpd_pair[i]  <-  looDiff[index[1,i], 3] - looDiff[index[2,i], 3]
    }

    # Standard Error of pairwise elpd_loo differences
    for (i in 1:ncol(index)) {
        selooDiff[i]  <-  sqrt(n * var(looList[[index[1,i]]]$pointwise[,"elpd_loo"]  - 
                                       looList[[index[2,i]]]$pointwise[,"elpd_loo"]))
    }

    # P-values for pairwise elpd_loo differences
        pDiff  <-  as.numeric(rounded(2*pnorm(-abs((elpd_pair - 0)/selooDiff)), 3))

    # Bayesian Bootstrap for differences
    for (i in 1:ncol(index)) {
      diffs         <-  bayesboot((looList[[index[1,i]]]$pointwise[,"elpd_loo"] - 
                                   looList[[index[2,i]]]$pointwise[,"elpd_loo"]), weighted.mean, use.weights=TRUE)
      bayesBootDiff[i]  <- sum(diffs < 0)/nrow(diffs)
      rm(diffs)
    }

    # Bayesian Bootstrap for s.e.
    for (i in 1:ncol(index)) {
      seSumm  <-  summary(bayesboot((looList[[index[1,i]]]$pointwise[,"elpd_loo"] - 
                                                    looList[[index[2,i]]]$pointwise[,"elpd_loo"]), R2=n, statistic=looDiffSE))
        bayesBootSE[i]    <-  seSumm$value[1]
        bayesBootSE95[i]  <-  seSumm$value[4]
    }

    # P-values for pairwise elpd_loo differences
        pDiffBootSE    <-  as.numeric(rounded(2*pnorm(-abs((elpd_pair - 0)/bayesBootSE)), 3))
        pDiffBootSE95  <-  as.numeric(rounded(2*pnorm(-abs((elpd_pair - 0)/bayesBootSE95)), 3))

    # Combine results into table
    LooDiff  <-  cbind(elpd_pair, selooDiff, pDiff, pDiffBootSE, pDiffBootSE95, bayesBootDiff, bayesBootSE)
    for(i in 1:ncol(index)){
        names[i]  <-  paste0(unlist(strsplit(nameIndex[index[1,i]],split='L'))[1]," v. ",
                             unlist(strsplit(nameIndex[index[2,i]],split='L'))[1])
    }
    row.names(LooDiff)  <-  names

    # Return result
    (LooDiff)
}

########################
# INVERSE LOGIT FUNCTION
########################
inv_logit <- function(u) {
    1 / (1 + exp(-u)); 
}

####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
whiteGrid  <-  function(...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}


#' Quick calcualtion of pvalue from stan sample data
#'
#' @title Quick calcualtion of pvalue from stan sample data
#' @param x relevant column of stan sample
#' @export
pval  <-  function(x) length(x[x < 0])/length(x)

#' Makes density plot for contrasts 
#'
#' @title Makes density plot for contrasts 
#' @param Dens density object that needs plotting
#' @param name name for plot
#' @export
plotContr  <-  function(Dens, name="title") {
  plot(NA, xlab=expression(paste(Delta)), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, xlim=c(min(Dens$x), (max(Dens$x)+0.4*(max(Dens$x) - min(Dens$x)))), ylim=c(0, (max(Dens$y)+0.05*(max(Dens$y) - min(Dens$y)))), yaxs='i')
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  polygon(c(Dens$x), c(Dens$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
  abline(v=0, lwd=2,col=2)
  axis(1, cex.axis=0.9)
  axis(2, cex.axis=0.9, las=1)
}
