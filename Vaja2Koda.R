
##:::::::::::::::::::::::::::::::::
## Project: Morphological opening and closing with reconstruction
## Author: Aleksandar Tendjer
## Description:
##Doing Mofrological opening and closing, opening and closing with reconstruction on 3 different images of at least 1000x1000 pixels and one multispectral index
##Bands 6 to 9 are in the Near Infrared Range (NIR)
##10 m resolution band 2, band 3, band 4 and band 8
##20 m resolution band 5, band 6, band 7, band 11 and band 12
##60 m resolution band 1, band 9 and band 10
##LINK to the 1 image from scihub : https://scihub.copernicus.eu/dhus/odata/v1/Products('c92416ce-4477-496c-8552-d96f2fb10519')/$value
##LINK to the 2 image from scihub :https://scihub.copernicus.eu/dhus/odata/v1/Products('35c3bff7-d385-48d0-9ce5-30f0e16bd2e4')/$value
##LINK to the 3 image from scihub : https://scihub.copernicus.eu/dhus/odata/v1/Products('7c25da3e-0c99-4b9b-91f2-ec53b1a35129')/$value
##:::::::::::::::::::::::::::::::::






###########################################################################
###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           FUNCTION DEFINITION                       ###
###                                                                     ###
###########################################################################
###########################################################################
##########HELP FUNCTIONS##########

lowFilter<-function(data, filter_size){
      temp = array()
    indexer = filter_size / 2
    for (i in length(data))
    {
      for (j in 1:length(data[0]))
        for (z in 1:filter_size)
          if (i + z - indexer < 0 || i + z - indexer > len(data) - 1)
          {
            for (c in 1:range(filter_size))
            temp.append(0)
          }
          else
            {
              if (j + z - indexer < 0 || j + indexer > len(data[0]) - 1)
              temp.append(0)
                else
              for (k in range(filter_size))
              temp.append(data[i + z - indexer][j + k - indexer])
            }
          temp.sort()
          data_final[i][j] = temp[len(temp) / 2]
          temp = null
    }
}
#spectral index MSI  function bands 8 and 11
msi_sentinel<-function(band8,band11){
  band8=raster::brick(band8)
  band11=raster::brick(band11)
  msi=overlay(band11, band8, fun=function(x,y){(x/y)})
  output_name= "msi.tiff"
  output_name=paste(num,output_name)
  
  #export the image to the working directory
  raster::writeRaster(msi, filename = output_name)
  
  return(msi)
}


### Load packages
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
#library that can represent images for people with low eye sight
library(viridis)
library(rasterVis)
library(sen2r)

###########################################################################
###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           DATA INPUT                                ###
###                                                                     ###
###########################################################################
### DESCRIPTION:we collect bands 8,11  out of three images from different###
### locations                                                           ###
###########################################################################
setwd("imgs/")

b8_img1<-"T34UDV_20190922T094031_B08_10m.jp2"
b11_img1<-"T34UDV_20190922T094031_B11_20m.jp2"

b8_img2<-"T14UNB_20200910T173911_B08_10m.jp2"
b11_img2<-"T14UMC_20200910T173911_B11_20m.jp2"

b8_img3<-"T10UGA_20200815T185921_B08_10m.jp2"
b11_img3<-"T10UGA_20200815T185921_B11_20m.jp2"


b8_raster3<-raster(readGDAL(b8_img3))
#make bigger resolution of pixel
b8_raster3_20m= aggregate(b8_raster3, fact=2)
rm(b8_raster3)
b11_raster3<-raster(readGDAL(b11_img3))



b8_raster2<-raster(readGDAL(b3_img2))
b8_raster3_20m= aggregate(b8_raster3, fact=2)
b11_raster2<-raster(readGDAL(b11_img2))


###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           DATA INPUT                                ###
###                                                                     ###
###########################################################################
### DESCRIPTION: We get the values of indices                               ###
###                                                                     ###
###########################################################################


msi_3=msi_sentinel(b8_raster3_20m,b11_raster3)