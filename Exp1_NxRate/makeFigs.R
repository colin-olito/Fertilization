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
source('./loadData.R')



######################
# FIGURES OF INTEREST
######################

toPdf(N_investPlot(), figPath(name='N_invest_Regression.pdf'), width=7, height=7)
embed_fonts(figPath(name='N_invest_Regression.pdf'))

toPdf(NxRate_Plot(data=data), figPath(name='NxRate_Regression.pdf'), width=7, height=7)
embed_fonts(figPath(name='NxRate_Regression.pdf'))


######################
# FIGURES FOR THE PAPER
######################

toPdf(regressionPlot(NinvData=NinvData, data=data), figPath(name='fertPlots.pdf'), width=7, height=14)
embed_fonts(figPath(name='fertPlots.pdf'))

toPdf(perGameteFertPlot(data=data), figPath(name='perGametePlot.pdf'), width=7, height=7)
embed_fonts(figPath(name='perGametePlot.pdf'))


########################
# SUPPLEMENTARY FIGURES
########################

toPdf(coefContrastPlots(), figPath(name='NxRate_coefContrasts.pdf'), width=14, height=7)
embed_fonts(figPath(name='NxRate_coefContrasts.pdf'))

toPdf(simpleContrastPlots(), figPath(name='NxRate_simpleContrasts.pdf'), width=16, height=14)
embed_fonts(figPath(name='NxRate_simpleContrasts.pdf'))

