#####################################################
#  Fertilization in broadcast spawners
#
#  Functions to generate Figures for: 
#    Article title goes here...
#  
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())

library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

#source('paths.R')
source('R/functions.R')
source('R/functions-figures.R')
source('./loadFigData.R')



######################
# FIGURES OF INTEREST
######################

toPdf(N_investPlot(NIm2.df, NIm2.summ, NinvData), figPath(name='N_invest_Regression.pdf'), width=7, height=7)
embed_fonts(figPath(name='N_invest_Regression.pdf'))

Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : Rate   : EggPos,
                         data = data)

toPdf(NxRate_Plot(df=m12.df, summ=m12.summ, Zmat = Z, data = data), figPath(name='NxRate_Regression.pdf'), width=7, height=7)
embed_fonts(figPath(name='NxRate_Regression.pdf'))
rm(Z)

Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z,
                         data = data)

toPdf(NxRate_Plot(df=m21BB.df, summ=m21BB.summ, Zmat = Z, data = data), figPath(name='NxRate_RegressionBB.pdf'), width=7, height=7)
embed_fonts(figPath(name='NxRate_RegressionBB.pdf'))
rm(Z)

Z       <-  model.matrix(~ -1 + Run +
                                Run : nSperm_z,
                         data = data)
toPdf(perGameteFertPlot(df = m12.df, summ=m12.summ, Zmat=Z, data=data), 
	  figPath(name='perGametePlotBB.pdf'), width=7, height=7)
embed_fonts(figPath(name='perGametePlotBB.pdf'))

######################
# FIGURES FOR THE PAPER
######################


Z       <-  model.matrix(~ -1 + Run            +
                                Run : nSperm_z +
                                Run : Rate     +
                                Run : EggPos   +
                                Run : Rate   : EggPos,
                         data = data)

toPdf(regressionPlot(Ninv.df=NIm2.df, Ninv.summ=NIm2.summ, 
	                 NRate.df=m12.df, NRate.summ=m12.summ, Zmat=Z, NinvData=NinvData, data=data), 
                     figPath(name='fertPlots.pdf'), width=7, height=14)
embed_fonts(figPath(name='fertPlots.pdf'))

toPdf(perGameteFertPlot(df = m12.df, summ=m12.summ, Zmat=Z, data=data), 
	  figPath(name='perGametePlot.pdf'), width=7, height=7)
embed_fonts(figPath(name='perGametePlot.pdf'))


########################
# SUPPLEMENTARY FIGURES
########################

toPdf(coefContrastPlots(df=m12.df, summ=m12.summ), figPath(name='NxRate_coefContrasts.pdf'), width=14, height=7)
embed_fonts(figPath(name='NxRate_coefContrasts.pdf'))

toPdf(simpleContrastPlots(df=m12.df, summ=m12.summ), figPath(name='NxRate_simpleContrasts.pdf'), width=16, height=14)
embed_fonts(figPath(name='NxRate_simpleContrasts.pdf'))

