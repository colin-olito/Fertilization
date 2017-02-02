###############
# DEPENDENCIES
###############
# All Dependencies are called from the R/functions.R file


#########################################################
# IMPORT DATA SETS AND STAN FILES NECESSARY FOR PLOTTING
#########################################################


###############################
# ANCILLARY PLOTTING FUNCTIONS
###############################



Fast5.ci   <-  function(x) {
    dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), 
                    to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
    inv_logit(x[1] + (x[2] * dummyX))
                  }

Fast55.ci   <-  function(x) {
    dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), 
                    to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
    inv_logit((x[1] + x[4]) + (x[2] + x[6]) * dummyX)
                  }

Fast.ci   <-  function(x) {
    dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                    to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
    inv_logit((x[1] + (x[4])/2) + (x[2] + (x[6])/2) * dummyX)
                  }

Slow.ci   <-  function(x) {
    dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                    to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
    inv_logit((x[1] + x[3] + (0.5*(x[7]))) + (x[2] + x[5] + (0.5*(x[8]))) * dummyX)
                  }

Fast5.plots  <-  function(betas, allBetas, gammas, Z, ..., data) {
    dummyX       <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), 
                    to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
    dummyXRaw    <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), 
                    to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
    y            <-  inv_logit(betas[1] + (betas[2] * dummyX))
    preds        <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,Fast5.ci)))
    CIs          <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
 	realYs       <-  (data$nFert - data$nControlFert)/data$nEggs 
 	yHats        <-  inv_logit(betas[1] + (betas[2] * data$nSperm_z))
    rawResids    <-  yHats - realYs
    adjResids    <-  inv_logit(betas[1] + (betas[2] * data$nSperm_z) + Z %*% gammas) - realYs 
    yAdj         <-  yHats  + adjResids
    list("y"          =  y,
         "yAdj"       =  yAdj,
         "yReal"      =  realYs,
         "xReal"      =  data$nSperm,
         "yHats"      =  yHats,
         "x"          =  dummyX,
         "xRaw"       =  dummyXRaw,
         "rawResids"  =  rawResids,
         "adjResids"  =  adjResids,
         "preds"      =  preds,
         "CIs"        =  CIs
        )
}

Fast55.plots  <-  function(betas, allBetas, gammas, Z, data) {
  dummyX       <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  dummyXRaw    <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  y            <-  inv_logit((betas[1] + betas[4]) + (betas[2] + ((betas[6])/2)) * dummyX)
  preds        <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,Fast55.ci)))
  CIs          <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
 	realYs     <-  ((data$nFert - data$nControlFert)/data$nEggs) 
 	yHats      <-  inv_logit((betas[1] + betas[4]) + (betas[2] + betas[6]) * data$nSperm_z)
    rawResids  <-  yHats - realYs
    adjResids  <-  inv_logit((betas[1] + betas[4]) + (betas[2] + betas[6]) * data$nSperm_z + Z %*% gammas) - realYs
    yAdj       <-  yHats  + adjResids
    list("y"          =  y,
         "yAdj"       =  yAdj,
         "yReal"      =  realYs,
         "xReal"      =  data$nSperm,
         "yHats"      =  yHats,
         "x"          =  dummyX,
         "xRaw"       =  dummyXRaw,
         "rawResids"  =  rawResids,
         "adjResids"  =  adjResids,
         "preds"      =  preds,
         "CIs"        =  CIs
        )
}

Fast.plots  <-  function(betas, allBetas, gammas, Z, data) {
    dummyX     <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                    to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
    dummyXRaw  <-  seq(from = min(data$nSperm[data$Rate == "Slow"]), 
                    to   = max(data$nSperm[data$Rate == "Slow"]), length=500)
    y          <-  inv_logit((betas[1] + (betas[4])/2) + (betas[2] + (betas[6])/2) * dummyX)
    preds      <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas), 1, Fast.ci)))
    CIs        <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
 	realYs     <-  ((data$nFert - data$nControlFert)/data$nEggs) 
 	yHats      <-  inv_logit((betas[1] + (betas[4])/2) + (betas[2] + (betas[6])/2) * data$nSperm_z)
    rawResids  <-  yHats - realYs
    adjResids  <-  inv_logit((betas[1] + (betas[4])/2) + ((betas[2] + (betas[6])/2) * data$nSperm_z) + Z %*% gammas) - realYs 
    yAdj       <-  yHats  + adjResids
    list("y"          =  y,
         "yAdj"       =  yAdj,
         "yReal"      =  realYs,
         "xReal"      =  data$nSperm,
         "yHats"      =  yHats,
         "x"          =  dummyX,
         "xRaw"       =  dummyXRaw,
         "rawResids"  =  rawResids,
         "adjResids"  =  adjResids,
         "preds"      =  preds,
         "CIs"        =  CIs
        )
}

Slow.plots  <-  function(betas, allBetas, gammas, Z, data) {
    dummyX     <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                    to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
    dummyXRaw  <-  seq(from = min(data$nSperm[data$Rate == "Slow"]), 
                    to   = max(data$nSperm[data$Rate == "Slow"]), length=500)
    y       <-  inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * dummyX)
    preds   <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,Slow.ci)))
    CIs     <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
 	realYs     <-  ((data$nFert - data$nControlFert)/data$nEggs) 
 	yHats      <-  inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * data$nSperm_z)
    rawResids  <-  yHats - realYs
    adjResids  <-  inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * data$nSperm_z + Z %*% gammas) - realYs 
    yAdj       <-  yHats  + adjResids
    list("y"          =  y,
         "yAdj"       =  yAdj,
         "yReal"      =  realYs,
         "xReal"      =  data$nSperm,
         "yHats"      =  yHats,
         "x"          =  dummyX,
         "xRaw"       =  dummyXRaw,
         "rawResids"  =  rawResids,
         "adjResids"  =  adjResids,
         "preds"      =  preds,
         "CIs"        =  CIs
        )
}

Invest.plots  <-  function(m2.summ, data) {
    dummyX     <-  seq(from = min(data$nSperm_z), to   = max(data$nSperm_z), length=500)
    dummyXRaw  <-  seq(from = min(data$nSperm), to   = max(data$nSperm), length=500)
    yHats       <-  inv_logit(m2.summ$Mean[1] + m2.summ$Mean[2] * dummyX)
    	x1  <-  seq(from = min(data$nSperm_z[data$Run == 1]), to = max(data$nSperm_z[data$Run == 1]), length=500)
    	x2  <-  seq(from = min(data$nSperm_z[data$Run == 2]), to = max(data$nSperm_z[data$Run == 2]), length=500)
    	x3  <-  seq(from = min(data$nSperm_z[data$Run == 3]), to = max(data$nSperm_z[data$Run == 3]), length=500)
    	x4  <-  seq(from = min(data$nSperm_z[data$Run == 4]), to = max(data$nSperm_z[data$Run == 4]), length=500)
    	x5  <-  seq(from = min(data$nSperm_z[data$Run == 5]), to = max(data$nSperm_z[data$Run == 5]), length=500)
    	x6  <-  seq(from = min(data$nSperm_z[data$Run == 6]), to = max(data$nSperm_z[data$Run == 6]), length=500)
    	x7  <-  seq(from = min(data$nSperm_z[data$Run == 7]), to = max(data$nSperm_z[data$Run == 7]), length=500)
    	x8  <-  seq(from = min(data$nSperm_z[data$Run == 8]), to = max(data$nSperm_z[data$Run == 8]), length=500)
    runs  <-  list(
                   Run1  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[3])  + m2.summ$Mean[2] * x1),
                   Run2  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[4])  + m2.summ$Mean[2] * x2),
                   Run3  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[5])  + m2.summ$Mean[2] * x3),
                   Run4  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[6])  + m2.summ$Mean[2] * x4),
                   Run5  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[7])  + m2.summ$Mean[2] * x5),
                   Run6  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[8])  + m2.summ$Mean[2] * x6),
                   Run7  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[9])  + m2.summ$Mean[2] * x7),
                   Run8  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[10]) + m2.summ$Mean[2] * x8),
                   x1    <-  seq(from = min(data$nSperm[data$Run == 1]), to = max(data$nSperm[data$Run == 1]), length=500),
                   x2    <-  seq(from = min(data$nSperm[data$Run == 2]), to = max(data$nSperm[data$Run == 2]), length=500),
                   x3    <-  seq(from = min(data$nSperm[data$Run == 3]), to = max(data$nSperm[data$Run == 3]), length=500),
                   x4    <-  seq(from = min(data$nSperm[data$Run == 4]), to = max(data$nSperm[data$Run == 4]), length=500),
                   x5    <-  seq(from = min(data$nSperm[data$Run == 5]), to = max(data$nSperm[data$Run == 5]), length=500),
                   x6    <-  seq(from = min(data$nSperm[data$Run == 6]), to = max(data$nSperm[data$Run == 6]), length=500),
                   x7    <-  seq(from = min(data$nSperm[data$Run == 7]), to = max(data$nSperm[data$Run == 7]), length=500),
                   x8    <-  seq(from = min(data$nSperm[data$Run == 8]), to = max(data$nSperm[data$Run == 8]), length=500)
                  )
    list("yHats"  =  yHats,
    	 "runs"   =  runs,
    	 "xReal"  =  data$nSperm,
    	 "x"      =  dummyX,
    	 "xRaw"   =  dummyXRaw
    	 )
}


###############################
# OUTPUT PLOTTING FUNCTIONS
###############################

#' Plot main regression from N_invest experiment
#'
#' @title Plot main regression from N_invest experiment
#' @param stanfit A stanfit object created from reading the relevant
#'                stan sample files.
#' @param dataLoc A character string giving the name of the data set
#'                used to fit the model (N_invest_master.csv).
#' @return A pdf of the final regression plot for this analysis.
#' @author Colin Olito.
N_investPlot  <-  function(stanfit = NIm2, data = NinvData) {

    ####################
    # Process Stanfit
    m2.df    <-  as.data.frame(extract(NIm2))
    m2.summ  <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]

    #######################
    # Make plotting object
    Inv.plt  <-  Invest.plots(m2.summ, data)

    ############
    # Make Plot
    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    plot((data$nFert/data$nEggs) ~ data$nSperm_z, 
        xlab='', ylab=substitute(''), 
        type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = data)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()

    # plot all regression lines from MCMC chains
    # apply(m2.df, 1, function(x, data, nSperm_z){
    #     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
    #     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
    #     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
    # }, data=data, nSperm_z=data$nSperm_z)
    # plot run-specific regression lines
     for(i in 1:8) {
       lines(Inv.plt$runs[[i]] ~ Inv.plt$runs[[i+8]],col='grey75', lwd=3)
     }
    # plot main regression line
    lines(Inv.plt$yHats ~ Inv.plt$xRaw, col='black', lwd=3)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=16, col=transparentColor('dodgerblue1', 0.7), cex=1.1)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=1,  col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Fertilization rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.02, 1.05, 'A', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
}







#' Plot main regression from NxRate experiment
#'
#' @title Plot main regression from NxRate experiment
#' @param stanfit A stanfit object created from reading the relevant
#'                stan sample files.
#' @param data A character string giving the name of the data set
#'             used to fit the model (NxRate_master.csv).
#' @return A pdf of the final regression plot for this analysis.
#' @author Colin Olito.
NxRate_Plot  <-  function(df = m12.df, summ = m12.summ, Zmat = m12Z, data, ...) {

    #######################
    # Extract Coefficients
    m12.betas      <-  summ$Mean[1:8]
    m12.allBetas   <-  df[1:8]
    m12.gammas     <-  summ$Mean[9:58]
    m12.allGammas  <-  df[9:58]

    ###################################################
    # Create plotting objects for each regression line
    Fast5.plt   <-  Fast5.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)
    Fast55.plt  <-  Fast55.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)
    Slow.plt    <-  Slow.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)

    ################
    # Make the plot
    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
        xlab='', ylab='', 
        type='n', axes=FALSE, ylim=c(0,1), xlim=c(min(data$nSperm),max(data$nSperm)))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    # plot regression lines
#    polygon(x=c(Slow.plt$xRaw, rev(Slow.plt$xRaw)), 
#            y=c(Slow.plt$CIs$lower, rev(Slow.plt$CIs$upper)), 
#            col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
#    polygon(x=c(Fast5.plt$xRaw, rev(Fast5.plt$xRaw)), 
#            y=c(Fast5.plt$CIs$lower, rev(Fast5.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
#    polygon(x=c(Fast55.plt$xRaw, rev(Fast55.plt$xRaw)), 
#            y=c(Fast55.plt$CIs$lower, rev(Fast55.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    lines(Slow.plt$y ~ Slow.plt$xRaw, col='orangered1', lwd=3)
    lines(Fast5.plt$y ~ Fast5.plt$xRaw, col='dodgerblue1', lwd=3)
    lines(Fast55.plt$y ~ Fast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
    points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=16, 
            col=transparentColor('orangered1', 0.7), cex=1.1)
    points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=1, 
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points(Fast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    points(Fast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ Fast55.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.2),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Adjusted fertilization rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.02, 1.05, 'B', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        legend(
              x       =  usr[2]*0.3,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(Fast:~5~cm)),
                          expression(paste(Fast:~55~cm)),
                          expression(paste(Slow))),
              pch     =  c(21,21,21),
              pt.bg   =  c(transparentColor('dodgerblue1',0.7),transparentColor('dodgerblue1',0.2),transparentColor('orangered1',0.7)),
              col     =  c('dodgerblue4','dodgerblue4','orangered4'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
}








###############################
# PLOT FUNCTIONS - FOR FIGURES
###############################

#' Plot main regression results from both N_invest & NxRate experiments
#'
#' @title Plot main regression results from both N_invest & NxRate experiments
#' @param stanfits A list of stanfit objects for the best models from the N_invest
#'                 and NxRate analyses
#' @param NinvData The NinvData data frame.
#' @param data The NxRate data frame (data).
#' @return A plot of the final regression for this analysis.
#' @author Colin Olito.
regressionPlot  <-  function(Ninv.fit = NIm2, NRate.df = m12.df, NRate.summ = m12.summ, Zmat = m12Z, NinvData, data) {

    #######################
    # Extract coefficients
    NIm2.df     <-  as.data.frame(extract(Ninv.fit))[,-1]
    NIm2.summ   <-  plyr:::adply(as.matrix(NIm2.df),2,MCMCsum)
    m12.betas      <-  NRate.summ$Mean[1:8]
    m12.allBetas   <-  NRate.df[1:8]
    m12.gammas     <-  NRate.summ$Mean[9:58]
    m12.allGammas  <-  NRate.df[9:58]

    ###################################################
    # Create plotting objects for each regression line
    Inv.plt     <-  Invest.plots(NIm2.summ, data=NinvData)
    Fast5.plt   <-  Fast5.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)
    Fast55.plt  <-  Fast55.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)
    Slow.plt    <-  Slow.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)

    ###################################################
    # MAKE PLOTS

    # Set plot layout
    layout.mat <- matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    ############
    # N_invPlot
    par(omi=rep(0.5, 4), mar = c(4,4,1,2), bty='o', xaxt='s', yaxt='s')
    plot((nFert/nEggs) ~ nSperm_z, xlab='', ylab=substitute(''), type='n', axes=FALSE, 
    	  ylim=c(0,1), xlim=c(min(nSperm),max(nSperm)), data = NinvData)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    # plot all regression lines from MCMC chains
    # apply(m2.df, 1, function(x, data, nSperm_z){
    #     xrange  <-  seq(min(data$nSperm), max(data$nSperm), length.out=100)
    #     xrange2  <-  seq(min(nSperm_z), max(nSperm_z), length.out=100)
    #     lines(xrange, inv_logit(x['beta.1'] + x['beta.2'] * xrange2), col=transparentColor('grey68',0.1))
    # }, data=NinvData, nSperm_z=NinvData$nSperm_z)
    # plot run-specific regression lines
     for(i in 1:8) {
       lines(Inv.plt$runs[[i]] ~ Inv.plt$runs[[i+8]],col='grey75', lwd=3)
     }
    # plot main regression line
    lines(Inv.plt$yHats ~ Inv.plt$xRaw, col='black', lwd=3)
    points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=16, col=transparentColor('dodgerblue1', 0.7), cex=1.1)
    points((NinvData$nFert/NinvData$nEggs) ~ NinvData$nSperm, pch=1,  col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Fertilization rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
#    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.02, 1.05, 'A', cex=1.5, adj=c(0.5, 0.5), xpd=NA)


    #############
    # NxRate Plot
    plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
        xlab='', ylab='', 
        type='n', axes=FALSE, ylim=c(-0.05,1), xlim=c(min(data$nSperm),max(data$nSperm)))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    # plot regression lines
#    polygon(x=c(Slow.plt$xRaw, rev(Slow.plt$xRaw)), 
#            y=c(Slow.plt$CIs$lower, rev(Slow.plt$CIs$upper)), 
#            col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
#    polygon(x=c(Fast5.plt$xRaw, rev(Fast5.plt$xRaw)), 
#            y=c(Fast5.plt$CIs$lower, rev(Fast5.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
#    polygon(x=c(Fast55.plt$xRaw, rev(Fast55.plt$xRaw)), 
#            y=c(Fast55.plt$CIs$lower, rev(Fast55.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    lines(Slow.plt$y ~ Slow.plt$xRaw, col='orangered1', lwd=3)
    lines(Fast5.plt$y ~ Fast5.plt$xRaw, col='dodgerblue1', lwd=3)
    lines(Fast55.plt$y ~ Fast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
    points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=16, 
            col=transparentColor('orangered1', 0.7), cex=1.1)
    points(Slow.plt$yAdj[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=1, 
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points(Fast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ Fast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    points(Fast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ Fast55.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.2),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Adjusted fertilization rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.02, 1.05, 'B', cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        legend(
              x       =  usr[2]*0.225,
              y       =  usr[4],
              legend  =  c(
                          expression(paste('')),
                          expression(paste('')),
                          expression(paste(''))),
              lty     =  c(1,2,1),
              lwd     =  c(2,2,2),
              seg.len =  3,
              col     =  c(
              	           'dodgerblue1',
              	           'dodgerblue1',
              	           'orangered1'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.398,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(~~~Fast:~5~cm)),
                          expression(paste(~~~Fast:~55~cm)),
                          expression(paste(~~~Slow))),
              pch     =  c(21,21,21),
              pt.bg   =  c(
              	           transparentColor('dodgerblue1',0.7),
              	           transparentColor('dodgerblue1',0.2),
              	           transparentColor('orangered1',0.7)),
              col     =  c(
              	           'dodgerblue4',
              	           'dodgerblue4',
              	           'orangered4'),
              pt.cex  =  1.25,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
}




#' Plot gamete fertilization success for pooled Fast & Slow regression results.
#'
#' @title Plot gamete fertilization success for pooled Fast & Slow regression results.
#' @param stanfit Stanfit object for the best model from the NxRate analysis
#' @param data The NxRate data frame (data).
#' @return A plot of the Fast & Slow regression lines with credible intervals, transformed
#'         into per-capita fertilization success.
#' @author Colin Olito.
perGameteFertPlot  <-  function(df = m12.df, summ = m12.summ, Zmat = m12Z, data, ...) {

    ####################
    #######################
    # Extract Coefficients
    m12.betas      <-  summ$Mean[1:8]
    m12.allBetas   <-  df[1:8]
    m12.gammas     <-  summ$Mean[9:58]
    m12.allGammas  <-  df[9:58]

    ###################################################
    # Create plotting objects for each regression line
    Fast.plt  <-  Fast.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)
    Slow.plt  <-  Slow.plots(m12.betas, m12.allBetas, m12.gammas, Z=Zmat, data = data)

    datPerGamete  <-  (((data$nFert - data$nControlFert)/data$nEggs) / data$nSperm)
    yRange        <-  range(datPerGamete)[2] - range(datPerGamete)[1]
    ################
    # Make the plot
    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    plot(datPerGamete ~ data$nSperm_z, 
        xlab='', ylab='', type='n', axes=FALSE, 
        ylim=c((min(datPerGamete) - (yRange*0.05)), (max(datPerGamete) + (yRange*0.15))), 
        xlim=c(min(data$nSperm),max(data$nSperm)))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    # plot regression lines
    polygon(x=c(Slow.plt$xRaw, rev(Slow.plt$xRaw)), 
            y=c((Slow.plt$CIs$X20 / Slow.plt$xRaw), rev((Slow.plt$CIs$X80 / Slow.plt$xRaw))), 
            col=transparentColor('orangered1', 0.05), border=transparentColor('orangered4',0.2))
    polygon(x=c(Fast.plt$xRaw, rev(Fast.plt$xRaw)), 
            y=c((Fast.plt$CIs$X20 / Fast.plt$xRaw), rev((Fast.plt$CIs$X80 / Fast.plt$xRaw))), 
            col=transparentColor('dodgerblue1', 0.05), border=transparentColor('dodgerblue4',0.2))
    lines((Slow.plt$y / Slow.plt$xRaw) ~ Slow.plt$xRaw, col='orangered1', lwd=3)
    lines((Fast.plt$y / Fast.plt$xRaw) ~ Fast.plt$xRaw, col='dodgerblue1', lwd=3, lty=1)
    points((Slow.plt$yAdj / Slow.plt$xReal)[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=16, 
            col=transparentColor('orangered1', 0.7), cex=1.1)
    points((Slow.plt$yAdj / Slow.plt$xReal)[data$Rate == "Slow"] ~ Slow.plt$xReal[data$Rate == "Slow"], pch=1, 
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points((Fast.plt$yAdj / Fast.plt$xReal)[data$Rate == "Fast"] ~ Fast.plt$xReal[data$Rate == "Fast"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7), col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.175, 0.5, expression(paste("Adjusted per-gamete fertilization rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
        legend(
              x       =  usr[2]*0.93,
              y       =  usr[4],
              legend  =  c(
                          expression(paste('')),
                          expression(paste(''))),
              lty     =  c(1,1),
              lwd     =  2,
              seg.len = 3,
              col     =  c(
              	           'dodgerblue1',
              	           'orangered1'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
        legend(
              x       =  usr[2]*0.9825,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(~~~Fast)),
                          expression(paste(~~~Slow))),
              pch     =  c(21,21),
              pt.bg   =  c(
              	           transparentColor('dodgerblue1',0.7),
              	           transparentColor('orangered1',0.7)),
              col     =  c(
              	           'dodgerblue4',
              	           'orangered4'),
              pt.cex  =  1.25,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
    )
}






#########################################
# PLOT FUNCTIONS - SUPPLEMENTARY FIGURES
#########################################



#' Multi-panel plot of Fixed Effects a-priori contrasts for the NxRate analysis.
#'
#' @title Multi-panel plot of Fixed Effects a-priori contrasts for the NxRate analysis.
#' @param df data frame with chosen model samples.
#' @param summ Corresponding MCMC summary of df.
#' @return A multipanel density plot showing simple contrasts.
#' @author Colin Olito.
coefContrastPlots  <-  function(df = m12.df, summ = m12.summ, ...) {

    #######################
    # Extract Coefficients
    m12.allBetas   <-  as.matrix(m12.df[1:8])

    # for model m12
    b0Fast    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2))
    b0Slow    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))))
    b0Fast5   <-  inv_logit(m12.allBetas[,1])
    b0Fast55  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]))
    b0Slow5   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]))
    b0Slow55  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]))
    b1Fast    <-  inv_logit((m12.allBetas[,2] + (m12.allBetas[,6])/2))
    b1Slow    <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))))
    b1Fast5   <-  inv_logit(m12.allBetas[,2])
    b1Fast55  <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,6]))
    b1Slow5   <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5]))
    b1Slow55  <-  inv_logit((m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]))
 
    #  Contrasts
    Contrasts  <-  list(
                        c1  =  b1Slow   - b1Fast,
                        c2  =  b1Fast5  - b1Fast55,
                        c3  =  b1Slow5  - b1Slow55,
                        c4  =  b1Fast5  - b1Slow5,
                        c5  =  b1Fast55 - b1Slow55,
                        c6  =  b1Fast55 - b1Slow,
                        c7  =  b1Fast5  - b1Slow
                        )
 
    #  P-values
    pVals  <-  list(
                    p1  =  pval(Contrasts[[1]]),
                    p2  =  pval(Contrasts[[2]]),
                    p3  =  pval(Contrasts[[3]]),
                    p4  =  pval(Contrasts[[4]]),
                    p5  =  pval(Contrasts[[5]]),
                    p6  =  pval(Contrasts[[6]]),
                    p7  =  pval(Contrasts[[7]])
                    )
    ###################################################
    # MAKE PLOTS

    # Set plot layout
    layout.mat <- matrix(c(1:8), nrow=2, ncol=4, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    par(omi=rep(0.5, 4), mar = c(4,4,3,3), bty='o', xaxt='s', yaxt='s')
    for(i in 1:7) {
        Dens  <-  density(Contrasts[[i]])
        plot(NA, xlab=expression(paste(Delta)), type='n', axes=FALSE, ylab='', cex.lab=1.2, 
             xlim=c(min(Dens$x), (max(Dens$x)+0.4*(max(Dens$x) - min(Dens$x)))), 
             ylim=c(0, (max(Dens$y)+0.05*(max(Dens$y) - min(Dens$y)))), yaxs='i')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        whiteGrid()
        box()
        polygon(c(Dens$x), c(Dens$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
        abline(v=0, lwd=2,col=2, xpd=FALSE)
        axis(1, cex.axis=0.9)
        axis(2, cex.axis=0.9, las=1)
        proportionalLabel(0.65, 0.9, substitute(p~"="~val, list(val = rounded(pVals[[i]], 2))), cex=1.5, adj=c(0, 0.5), font=3, xpd=NA)
        
        #  Plot Labels
        if(i == 1) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Slow")]~-~beta["1,"~italic("Fast")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        }
        if(i == 2) 
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Fast,5cm")]~-~beta["1,"~italic("Fast,55cm")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        if(i == 3) 
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Slow,5cm")]~-~beta["1,"~italic("Slow,55cm")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        if(i == 4) 
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Fast,5cm")]~-~beta["1,"~italic("Slow,5cm")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        if(i == 5) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Fast,55cm")]~-~beta["1,"~italic("Slow,55cm")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        }
        if(i == 6) 
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Fast,55cm")]~-~beta["1,"~italic("Slow")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        if(i == 7) 
            proportionalLabel(0.5, 1.1, expression(paste(beta["1,"~italic("Fast,5cm")]~-~beta["1,"~italic("Slow")])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
    }
}





#' Multi-panel plot of simple contrasts for the NxRate regression lines
#'
#' @title Multi-panel plot of simple contrasts for the NxRate regression lines
#' @param df data frame with chosen model samples
#' @param summ Corresponding MCMC summary of df
#' @return A multipanel density plot showing simple contrasts
#' @author Colin Olito.
simpleContrastPlots  <-  function(df = m12.df, summ = m12.summ, ...) {

    #######################
    # Extract Coefficients
    m12.allBetas   <-  as.matrix(m12.df[1:8])
    m12.allGammas  <-  as.matrix(m12.df[9:58])

    ############################
    # Calculate predicted lines

    #  Simple Contrasts at -1, 0, 1, & 2 standard deviations of nSperm
    m12Fast5.neg2   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (-2))
    m12Fast55.neg2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (-2))
    m12Slow5.neg2   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (-2))
    m12Slow55.neg2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (-2))
    m12Fast.neg2    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (-2))
    m12Slow.neg2    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (-2))

    m12Fast5.neg1   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (-1))
    m12Fast55.neg1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (-1))
    m12Slow5.neg1   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (-1))
    m12Slow55.neg1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (-1))
    m12Fast.neg1    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (-1))
    m12Slow.neg1    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (-1))

    m12Fast5.0   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (0))
    m12Fast55.0  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (0))
    m12Slow5.0   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (0))
    m12Slow55.0  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (0))
    m12Fast.0    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (0))
    m12Slow.0    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (0))

    m12Fast5.1   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (1))
    m12Fast55.1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (1))
    m12Slow5.1   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (1))
    m12Slow55.1  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (1))
    m12Fast.1    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (1))
    m12Slow.1    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (1))

    m12Fast5.2   <-  inv_logit(m12.allBetas[,1] + m12.allBetas[,2] * (2))
    m12Fast55.2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,4]) + (m12.allBetas[,2] + m12.allBetas[,6]) * (2))
    m12Slow5.2   <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3]) + (m12.allBetas[,2] + m12.allBetas[,5]) * (2))
    m12Slow55.2  <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + m12.allBetas[,7]) + (m12.allBetas[,2] + m12.allBetas[,5] + m12.allBetas[,8]) * (2))
    m12Fast.2    <-  inv_logit((m12.allBetas[,1] + (m12.allBetas[,4])/2) + (m12.allBetas[,2] + (m12.allBetas[,6])/2) * (2))
    m12Slow.2    <-  inv_logit((m12.allBetas[,1] + m12.allBetas[,3] + (0.5*(m12.allBetas[,7]))) + (m12.allBetas[,2] + m12.allBetas[,5] + (0.5*(m12.allBetas[,8]))) * (2))


    # Perform contrasts
    simpContr  <-  list(
      cSimp1   =  m12Fast5.neg1 - m12Fast55.neg1,
      cSimp2   =  m12Slow5.neg1 - m12Slow55.neg1,
      cSimp3   =  m12Fast5.neg1 - m12Slow5.neg1,
      cSimp4   =  m12Fast55.neg1 - m12Slow55.neg1,
      cSimp5   =  m12Fast.neg1 - m12Slow.neg1,
      cSimp6   =  m12Fast5.0 - m12Fast55.0,
      cSimp7   =  m12Slow5.0 - m12Slow55.0,
      cSimp8   =  m12Fast5.0 - m12Slow5.0,
      cSimp9   =  m12Fast55.0 - m12Slow55.0,
      cSimp10  =  m12Fast.0 - m12Slow.0,
      cSimp11  =  m12Fast5.1 - m12Fast55.1,
      cSimp12  =  m12Slow5.1 - m12Slow55.1,
      cSimp13  =  m12Fast5.1 - m12Slow5.1,
      cSimp14  =  m12Fast55.1 - m12Slow55.1,
      cSimp15  =  m12Fast.1 - m12Slow.1,
      cSimp16  =  m12Fast5.2 - m12Fast55.2,
      cSimp17  =  m12Slow5.2 - m12Slow55.2,
      cSimp18  =  m12Fast5.2 - m12Slow5.2,
      cSimp19  =  m12Fast55.2 - m12Slow55.2,
      cSimp20  =  m12Fast.2 - m12Slow.2
    )

    # Calcuate P-values
    simpPVals  <-  list(
    pVal1   =  pval(simpContr[[1]]),
    pVal2   =  pval(simpContr[[2]]),
    pVal3   =  pval(simpContr[[3]]),
    pVal4   =  pval(simpContr[[4]]),
    pVal5   =  pval(simpContr[[5]]),
    pVal6   =  pval(simpContr[[6]]),
    pVal7   =  pval(simpContr[[7]]),
    pVal8   =  pval(simpContr[[8]]),
    pVal9   =  pval(simpContr[[9]]),
    pVal10  =  pval(simpContr[[10]]),
    pVal11  =  pval(simpContr[[11]]),
    pVal12  =  pval(simpContr[[12]]),
    pVal13  =  pval(simpContr[[13]]),
    pVal14  =  pval(simpContr[[14]]),
    pVal15  =  pval(simpContr[[15]]),
    pVal16  =  pval(simpContr[[16]]),
    pVal17  =  pval(simpContr[[17]]),
    pVal18  =  pval(simpContr[[18]]),
    pVal19  =  pval(simpContr[[19]]),
    pVal20  =  pval(simpContr[[20]])
    )

    ###################################################
    # MAKE PLOTS

    # Set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    par(omi=rep(0.5, 4), mar = c(4,4,1,2), bty='o', xaxt='s', yaxt='s')
    for(i in 1:20) {
        Dens  <-  density(simpContr[[i]])
        plot(NA, xlab=expression(paste(Delta)), type='n', axes=FALSE, ylab='', cex.lab=1.2, 
             xlim=c(min(Dens$x), (max(Dens$x)+0.4*(max(Dens$x) - min(Dens$x)))), 
             ylim=c(0, (max(Dens$y)+0.05*(max(Dens$y) - min(Dens$y)))), yaxs='i')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        whiteGrid()
        box()
        polygon(c(Dens$x), c(Dens$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
        abline(v=0, lwd=2,col=2, xpd=FALSE)
        axis(1, cex.axis=0.9)
        axis(2, cex.axis=0.9, las=1)
        proportionalLabel(0.65, 0.9, substitute(p~"="~val, list(val = rounded(simpPVals[[i]], 2))), cex=1.5, adj=c(0, 0.5), font=3, xpd=NA)
        
        #  Plot Labels
        if(i == 1) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Fast.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(-0.4, 0.5, expression(paste('x = -1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        }
        if(i == 2) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow.5 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 3) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.5 - Slow.5')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 4) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.2, expression(paste('Distribution of Fast.55 - Slow.55')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 5) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(0.5, 1.2, expression(paste('Distribution of Slow - Fast')), xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 6) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(-0.4, 0.5, expression(paste('x = 0',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 11) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(-0.4, 0.5, expression(paste('x = 1',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            }
        if(i == 16) {
            proportionalLabel(-0.25, 0.5, expression(paste('Density')), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
            proportionalLabel(-0.4, 0.5, expression(paste('x = 2',sigma)), srt=90, xpd=NA, adj=c(0.5, 0.5), font=3, cex=1.5)
        }
    }
}

