###############
# DEPENDENCIES
###############
# All Dependencies are called from the R/functions.R file



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

m3aFast5.plots  <-  function(betas, allBetas, data) {
  dummyX       <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
  dummyXRaw    <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), 
                  to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "5"]), length=500)
  y            <-  inv_logit((m3a.betas[1] + (m3a.betas[2] + (m3a.betas[6])/2) * dummyX))
  preds        <-  unname(as.matrix(plyr:::aaply(as.matrix(m3a.allBetas),1,m3aFast5.ci)))
  CIs          <-  plyr:::adply(preds, 2, MCMCsum)
  list("y"      =  y,
       "x"      =  dummyX,
       "xRaw"   =  dummyXRaw,
       "preds"  =  preds,
       "CIs"    =  CIs
      )
}
m3aFast55.plots  <-  function(betas, allBetas, data) {
  dummyX     <-  seq(from = min(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm_z[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  dummyXRaw  <-  seq(from = min(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), 
                  to   = max(data$nSperm[data$Rate == "Fast" & data$EggPos == "55"]), length=500);
  y       <-  inv_logit((m3a.betas[1] + m3a.betas[4]) + (m3a.betas[2] + ((m3a.betas[6])/2)) * dummyX)
  preds   <-  unname(as.matrix(plyr:::aaply(as.matrix(m3a.allBetas),1,m3aFast55.ci)))
  CIs     <-  plyr:::adply(preds, 2, MCMCsum)
  list("y"      =  y,
       "x"      =  dummyX,
       "xRaw"   =  dummyXRaw,
       "preds"  =  preds,
       "CIs"    =  CIs
      )
}
m3aSlow.plots  <-  function(betas, allBetas, data) {
  dummyX     <-  seq(from = min(data$nSperm_z[data$Rate == "Slow"]), 
                  to   = max(data$nSperm_z[data$Rate == "Slow"]), length=500)
  dummyXRaw  <-  seq(from = min(data$nSperm[data$Rate == "Slow"]), 
                  to   = max(data$nSperm[data$Rate == "Slow"]), length=500)
  y       <-  inv_logit((m3a.betas[1] + m3a.betas[3] + (0.5*(m3a.betas[7]))) + (m3a.betas[2] + m3a.betas[5] + (0.5*(m3a.betas[8]))) * dummyX)
  preds   <-  unname(as.matrix(plyr:::aaply(as.matrix(m3a.allBetas),1,m3aSlow.ci)))
  CIs     <-  plyr:::adply(preds, 2, MCMCsum)
  list("y"      =  y,
       "x"      =  dummyX,
       "xRaw"   =  dummyXRaw,
       "preds"  =  preds,
       "CIs"    =  CIs
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
N_investPlot  <-  function(stanfit = 'N_invest_m2', dataFile = 'Ninvest_master.csv') {

    ##############
    # Import Data
    dataLoc  <-  paste0('./data/', datafile)
    data <- read.csv(dataLoc, header=TRUE, stringsAsFactors=FALSE)
    data <- data.frame(data)
        
    # Convert grouping variables to factors; Correct Dates
    data$Run       <-  factor(data$Run)
    data$Colony    <-  factor(data$Colony)
    data$N         <-  factor(data$N)
    data$Lane      <-  factor(data$Lane)
    data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
    data$Date      <-  dmy(data$Date)
    data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

    ####################
    # Import Stan Files
    stanLoc   <-  paste0('./output/StanFits/')
    csvFiles  <-  c(paste0(stanLoc, stanfit, '.csv1'),
                    paste0(stanLoc, stanfit, '.csv2'),
                    paste0(stanLoc, stanfit, '.csv3'))
    m2        <-  read_stan_csv(csvFiles, col_major = TRUE)
    rm(csvFiles)

    ####################
    # Process Stanfit
    m2.df    <-  as.data.frame(extract(m2))
    runs  <-  list(
               Run1  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[3])  + m2.summ$Mean[2] * data$nSperm_z),
               Run2  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[4])  + m2.summ$Mean[2] * data$nSperm_z),
               Run3  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[5])  + m2.summ$Mean[2] * data$nSperm_z),
               Run4  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[6])  + m2.summ$Mean[2] * data$nSperm_z),
               Run5  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[7])  + m2.summ$Mean[2] * data$nSperm_z),
               Run6  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[8])  + m2.summ$Mean[2] * data$nSperm_z),
               Run7  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[9])  + m2.summ$Mean[2] * data$nSperm_z),
               Run8  <- inv_logit((m2.summ$Mean[1] + m2.summ$Mean[10]) + m2.summ$Mean[2] * data$nSperm_z)
              )
    RegLine  <-  inv_logit(m2.summ$Mean[1] + m2.summ$Mean[2] * data$nSperm_z)

    ############
    # Make Plot
    par(omi=rep(0.3, 4))
    plot((data$nFert/data$nEggs) ~ nSperm_z, 
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
    # }, data=data, nSperm_z=nSperm_z)
    # plot run-specific regression lines
     for(i in 1:8) {
       lines(runs[[i]][data$Run == i][order(data$nSperm_z[data$Run == i])] ~ data$nSperm[data$Run == i][order(data$nSperm_z[data$Run == i])],
                       col='grey75', lwd=3)
     }
    # plot main regression line
    lines(RegLine[order(data$nSperm_z)] ~ data$nSperm[order(data$nSperm_z)],
                      col='black', lwd=3)
    points((data$nFert/data$nEggs) ~ data$nSperm, pch=21, 
            bg=transparentColor('dodgerblue3', 0.7),
            col=transparentColor('dodgerblue1', 0.7), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
    proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
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
NxRate_Plot  <-  function(stanfit = 'NxRate_m1', dataFile = 'NxRate_master.csv') {

    ##############
	# Import Data
    dataLoc  <-  paste0('./data/', dataFile)
    data <- read.csv(dataLoc, header=TRUE, stringsAsFactors=FALSE)
    data <- data.frame(data)
    
    # Convert grouping variables to factors; Correct Dates
    data$Run       <-  factor(data$Run)
    data$Colony    <-  factor(data$Colony)
    data$N         <-  factor(data$N)
    data$Rate      <-  factor(data$Rate)
    data$EggPos    <-  factor(data$EggPos)
    data$Lane      <-  factor(data$Lane)
    data$nSperm_c  <-  data$nSperm - mean(data$nSperm)
    data$Date      <-  dmy(data$Date)
    data$nSperm_z  <-  (data$nSperm - mean(data$nSperm))/sd(data$nSperm)

    ####################
    # Import Stan Files
    stanLoc   <-  paste0('./output/StanFits/')
    csvFiles  <-  c(paste0(stanLoc, stanfit, '.csv1'),
                    paste0(stanLoc, stanfit, '.csv2'),
                    paste0(stanLoc, stanfit, '.csv3'))
    m3a       <-  read_stan_csv(csvFiles, col_major = TRUE)
    rm(csvFiles)

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
    m3aFast5.plt   <-  m3aFast5.plots(m3a.betas, m3a.allBetas, data)
    m3aFast55.plt  <-  m3aFast55.plots(m3a.betas, m3a.allBetas, data)
    m3aSlow.plt    <-  m3aSlow.plots(m3a.betas, m3a.allBetas, data)

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
    polygon(x=c(m3aSlow.plt$xRaw, rev(m3aSlow.plt$xRaw)), 
            y=c(m3aSlow.plt$CIs$lower, rev(m3aSlow.plt$CIs$upper)), 
            col=transparentColor('orangered1', 0.01), border=transparentColor('orangered4',0.2))
    polygon(x=c(m3aFast5.plt$xRaw, rev(m3aFast5.plt$xRaw)), 
            y=c(m3aFast5.plt$CIs$lower, rev(m3aFast5.plt$CIs$upper)), 
            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    polygon(x=c(m3aFast55.plt$xRaw, rev(m3aFast55.plt$xRaw)), 
            y=c(m3aFast55.plt$CIs$lower, rev(m3aFast55.plt$CIs$upper)), 
            col=transparentColor('dodgerblue1', 0.01), border=transparentColor('dodgerblue4',0.2))
    lines(m3aSlow.plt$y ~ m3aSlow.plt$xRaw, col='orangered1', lwd=3)
    lines(m3aFast5.plt$y ~ m3aFast5.plt$xRaw, col='dodgerblue1', lwd=3)
    lines(m3aFast55.plt$y ~ m3aFast55.plt$xRaw, col='dodgerblue1', lwd=3, lty=2)
    points(m3aSlow_yAdj[data$Rate == "Slow"] ~ data$nSperm[data$Rate == "Slow"], pch=21, 
            bg=transparentColor('orangered1', 0.7),
            col=transparentColor('orangered4', 0.9), cex=1.1)
    points(m3aFast5_yAdj[data$Rate == "Fast" & data$EggPos == "5"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "5"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.7),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    points(m3aFast55_yAdj[data$Rate == "Fast" & data$EggPos == "55"] ~ data$nSperm[data$Rate == "Fast" & data$EggPos == "55"], pch=21, 
            bg=transparentColor('dodgerblue1', 0.2),
            col=transparentColor('dodgerblue4', 0.9), cex=1.1)
    axis(2, las=1)
    axis(1)
    proportionalLabel(-0.15, 0.5, expression(paste("Adjusted Fertilization Rate")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
    proportionalLabel(0.5, -0.15, expression(paste("Sperm Released")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)
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