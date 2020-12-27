
##:::::::::::::::::::::::::::::::::
## Project: Receiving Sentinel 2 images and finding  various spatial indices 
## Author: Aleksandar Tendjer
## Description:
##Doing Mofrological opening and closing, opening and closing with reconstruction on 3 different images of at least 1000x1000 pixels and one multispectral index
##Bands 6 to 9 are in the Near Infrared Range (NIR)
##10 m resolution band 2, band 3, band 4 and band 8
##20 m resolution band 5, band 6, band 7, band 11 and band 12
##60 m resolution band 1, band 9 and band 10
##LINK to the image from scihub : https://scihub.copernicus.eu/dhus/odata/v1/Products('c92416ce-4477-496c-8552-d96f2fb10519')/$value
##:::::::::::::::::::::::::::::::::
#install.packages(c("sp", "rgdal", "raster", "viridis", "rasterVis"))





###########################################################################
###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           FUNCTION DEFINITION                       ###
###                                                                     ###
###########################################################################
###########################################################################
##########HELP FUNCTIONS##########





### Load packages
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
#library that can represent images for people with low eye sight
library(viridis)
library(rasterVis)
library(sen2r)

