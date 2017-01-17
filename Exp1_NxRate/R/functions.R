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

################
# STAN OPTIONS
################
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################
# SUMMARIZE MCMC RESULTS
########################
MCMCsum <- function(x) {
   data.frame(Mean=mean(x, na.rm=TRUE) , Median=median(x, na.rm=TRUE), t(quantile(x,na.rm=TRUE)), HPDinterval(as.mcmc(x)))
}

#######################
# SUMMARIZE LOO RESULTS
#######################
#' Summarize loo model comparison results
#'
#' @title Make a table to summarize LOO model comparisons.
#' @param looDiff Result table from baked in function compare() from the 'loo' package.
#' @param looList A list of the loo objects for each of the models being compared.
#'                CRITICAL: order of these loo objects must match the ranked row.names from
#'                          the associated looDiff object!!!!
#' @return A Results table in matrix form with pairwise elpd_loo differences, s.e., and p-values.
#' @author Colin Olito.

#' @export
makeLooTable <- function(looDiff, looList) {
    # Variables, Containers
    index      <-  combn(nrow(looDiff),2)
    nameIndex  <-  row.names(looDiff)
    elpd_pair  <-  rep(0, length = ncol(index))
    selooDiff  <-  rep(0, length = ncol(index))
    pDiff      <-  rep(0, length = ncol(index))
    names      <-  rep(0, length = ncol(index))
    
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
#    for (i in 1:ncol(index)) {
        pDiff  <-  as.numeric(rounded(2*pnorm(-abs((elpd_pair - 0)/selooDiff)), 3))
#    }

    # Combine results into table
    LooDiff  <-  cbind(elpd_pair, selooDiff, pDiff)
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