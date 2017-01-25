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


m3aFast5.ci   <-  function(x) {
  dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
  inv_logit((x[1] + (x[2] + (x[6])/2) * dummyX))
                  }
m3aFast55.ci   <-  function(x) {
  dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  inv_logit((x[1] + x[4]) + (x[2] + ((x[6])/2)) * dummyX)
                  }
m3aSlow.ci   <-  function(x) {
  dummyX  <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                  to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
  inv_logit((x[1] + x[3] + (0.5*(x[7]))) + (x[2] + x[5] + (0.5*(x[8]))) * dummyX)
                  }

m3aFast5.plots  <-  function(betas, allBetas, gammas, data) {
    dummyX       <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), 
                    to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
    dummyXRaw    <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), 
                    to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
    y            <-  inv_logit((betas[1] + (betas[2] + (betas[6])/2) * dummyX))
    preds        <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,m3aFast5.ci)))
    CIs          <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
    Z            <-  model.matrix(~ -1 + data$Run + data$Run:data$nSperm_z)
 	realYs       <-  (data$nFert - data$nControlFert)/data$nEggs 
 	yHats        <-  inv_logit((betas[1] + (betas[2] + (betas[6])/2) * data$nSperm_z))
    yAdj         <-  yHats  + inv_logit((betas[1] + (betas[2] + (betas[6])/2) * data$nSperm_z) + Z %*% gammas) - realYs 
    list("y"      =  y,
         "yAdj"   =  yAdj,
         "xReal"  =  data$nSperm,
         "x"      =  dummyX,
         "xRaw"   =  dummyXRaw,
         "preds"  =  preds,
         "CIs"    =  CIs
        )
}

m3aFast55.plots  <-  function(betas, allBetas, gammas, data) {
  dummyX       <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  dummyXRaw    <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  y            <-  inv_logit((betas[1] + betas[4]) + (betas[2] + ((betas[6])/2)) * dummyX)
  preds        <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,m3aFast55.ci)))
  CIs          <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
    Z          <-  model.matrix(~ -1 + data$Run + data$Run:data$nSperm_z)
 	realYs     <-  ((data$nFert - data$nControlFert)/data$nEggs) 
 	yHats      <-  inv_logit((betas[1] + betas[4]) + (betas[2] + ((betas[6])/2)) * data$nSperm_z)
    yAdj       <-  yHats  + inv_logit((betas[1] + betas[4]) + (betas[2] + ((betas[6])/2)) * data$nSperm_z + Z %*% gammas) - realYs 
    list("y"      =  y,
         "yAdj"   =  yAdj,
         "xReal"  =  data$nSperm,
         "x"      =  dummyX,
         "xRaw"   =  dummyXRaw,
         "preds"  =  preds,
         "CIs"    =  CIs
        )
}

m3aSlow.plots  <-  function(betas, allBetas, gammas, data) {
    dummyX     <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                    to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
    dummyXRaw  <-  seq(from = min(data$nSperm[data$Rate == "Slow"]), 
                    to   = max(data$nSperm[data$Rate == "Slow"]), length=500)
    y       <-  inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * dummyX)
    preds   <-  unname(as.matrix(plyr:::aaply(as.matrix(allBetas),1,m3aSlow.ci)))
    CIs     <-  plyr:::adply(preds, 2, MCMCsum)
    # Calculate adjusted y-values
    Z          <-  model.matrix(~ -1 + data$Run + data$Run:data$nSperm_z)
 	realYs     <-  ((data$nFert - data$nControlFert)/data$nEggs) 
 	yHats      <-  inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * data$nSperm_z)
    yAdj       <-  yHats  + inv_logit((betas[1] + betas[3] + (0.5*(betas[7]))) + (betas[2] + betas[5] + (0.5*(betas[8]))) * data$nSperm_z + Z %*% gammas) - realYs 
    list("y"      =  y,
         "yAdj"   =  yAdj,
         "xReal"  =  data$nSperm,
         "x"      =  dummyX,
         "xRaw"   =  dummyXRaw,
         "preds"  =  preds,
         "CIs"    =  CIs
        )
}


m2Invest.plots  <-  function(m2.summ, data) {
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
    m2.plt  <-  m2Invest.plots(m2.summ, data)

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
       lines(m2.plt$runs[[i]] ~ m2.plt$runs[[i+8]],col='grey75', lwd=3)
     }
    # plot main regression line
    lines(m2.plt$yHats ~ m2.plt$xRaw, col='black', lwd=3)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=16, col=transparentColor('dodgerblue1', 0.7), cex=1.1)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=1,  col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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
NxRate_Plot  <-  function(stanfit = m3a, data) {

    ####################
    # Process Stanfit
    m3a.df         <-  as.data.frame(extract(m3a))
    m3a.summ       <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)[-1,]
    m3a.allBetas   <-  m3a.df[2:9]
    m3a.betas      <-  m3a.summ$Mean[1:8]
    m3a.allBetas   <-  m3a.df[2:9]
    m3a.gammas     <-  m3a.summ$Mean[9:28]
    m3a.allGammas  <-  m3a.df[10:29]

    ###################################################
    # Create plotting objects for each regression line
    m3aFast5.plt   <-  m3aFast5.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)
    m3aFast55.plt  <-  m3aFast55.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)
    m3aSlow.plt    <-  m3aSlow.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)

    ################
    # Make the plot
    par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
    plot(((data$nFert - data$nControlFert)/data$nEggs) ~ data$nSperm_z, 
        xlab='', ylab='', 
        type='n', axes=FALSE, ylim=c(-0.05,1), xlim=c(min(data$nSperm),max(data$nSperm)))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    # plot regression lines
#    polygon(x=c(m3aSlow.plt$xRaw, rev(m3aSlow.plt$xRaw)), 
#            y=c(m3aSlow.plt$CIs$lower, rev(m3aSlow.plt$CIs$upper)), 
#            col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
#    polygon(x=c(m3aFast5.plt$xRaw, rev(m3aFast5.plt$xRaw)), 
#            y=c(m3aFast5.plt$CIs$lower, rev(m3aFast5.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
#    polygon(x=c(m3aFast55.plt$xRaw, rev(m3aFast55.plt$xRaw)), 
#            y=c(m3aFast55.plt$CIs$lower, rev(m3aFast55.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    lines(m3aSlow.plt$y ~ m3aSlow.plt$xRaw, col='orangered1', lwd=3)
    lines(m3aFast5.plt$y ~ m3aFast5.plt$xRaw, col='dodgerblue1', lwd=3)
    lines(m3aFast55.plt$y ~ m3aFast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
    points(m3aSlow.plt$yAdj[data$Rate == "Slow"] ~ m3aSlow.plt$xReal[data$Rate == "Slow"], pch=16, 
            col=transparentColor('orangered1', 0.7), cex=1.1)
    points(m3aSlow.plt$yAdj[data$Rate == "Slow"] ~ m3aSlow.plt$xReal[data$Rate == "Slow"], pch=1, 
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points(m3aFast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ m3aFast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    points(m3aFast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ m3aFast55.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.2),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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
regressionPlot  <-  function(stanfits = list(NIm2, m3a), NinvData, data) {

    ####################
    # Process Stanfits
    m2.df          <-  as.data.frame(extract(NIm2))
    m2.summ        <-  plyr:::adply(as.matrix(m2.df),2,MCMCsum)[-1,]
    m3a.df         <-  as.data.frame(extract(m3a))
    m3a.summ       <-  plyr:::adply(as.matrix(m3a.df),2,MCMCsum)[-1,]
    m3a.allBetas   <-  m3a.df[2:9]
    m3a.betas      <-  m3a.summ$Mean[1:8]
    m3a.allBetas   <-  m3a.df[2:9]
    m3a.gammas     <-  m3a.summ$Mean[9:28]
    m3a.allGammas  <-  m3a.df[10:29]

    ###################################################
    # Create plotting objects for each regression line
    m2.plt         <-  m2Invest.plots(m2.summ, data)
    m3aFast5.plt   <-  m3aFast5.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)
    m3aFast55.plt  <-  m3aFast55.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)
    m3aSlow.plt    <-  m3aSlow.plots(m3a.betas, m3a.allBetas, m3a.gammas, data)


    # Set plot layout
    layout.mat <- matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    ############
    # N_invPlot
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
       lines(m2.plt$runs[[i]] ~ m2.plt$runs[[i+8]],col='grey75', lwd=3)
     }
    # plot main regression line
    lines(m2.plt$yHats ~ m2.plt$xRaw, col='black', lwd=3)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=16, col=transparentColor('dodgerblue1', 0.7), cex=1.1)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=1,  col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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
#    polygon(x=c(m3aSlow.plt$xRaw, rev(m3aSlow.plt$xRaw)), 
#            y=c(m3aSlow.plt$CIs$lower, rev(m3aSlow.plt$CIs$upper)), 
#            col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
#    polygon(x=c(m3aFast5.plt$xRaw, rev(m3aFast5.plt$xRaw)), 
#            y=c(m3aFast5.plt$CIs$lower, rev(m3aFast5.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
#    polygon(x=c(m3aFast55.plt$xRaw, rev(m3aFast55.plt$xRaw)), 
#            y=c(m3aFast55.plt$CIs$lower, rev(m3aFast55.plt$CIs$upper)), 
#            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    lines(m3aSlow.plt$y ~ m3aSlow.plt$xRaw, col='orangered1', lwd=3)
    lines(m3aFast5.plt$y ~ m3aFast5.plt$xRaw, col='dodgerblue1', lwd=3)
    lines(m3aFast55.plt$y ~ m3aFast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
    points(m3aSlow.plt$yAdj[data$Rate == "Slow"] ~ m3aSlow.plt$xReal[data$Rate == "Slow"], pch=16, 
            col=transparentColor('orangered1', 0.7), cex=1.1)
    points(m3aSlow.plt$yAdj[data$Rate == "Slow"] ~ m3aSlow.plt$xReal[data$Rate == "Slow"], pch=1, 
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points(m3aFast5.plt$yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ m3aFast5.plt$xReal[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    points(m3aFast55.plt$yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ m3aFast55.plt$xReal[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.2),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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
