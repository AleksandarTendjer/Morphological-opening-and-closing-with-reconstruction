
##:::::::::::::::::::::::::::::::::
## Project: Morphological opening and closing with reconstruction
## Author: Aleksandar Tendjer
## Description:
##Doing Mofrological opening and closing, opening and closing with reconstruction on 3 different images of at least 1000x1000 pixels and one multispectral index
##n opening is an erosion followed by a dilation, while a closing is a dilation followed by an erosion.
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
##########FUNCTIONS##########
extract_matrix<-function(i,j,img,kernel){
  #extract neighbouring pixel values to calculate erosion
  x=2*kernel+1
  extract_mat=matrix(data=0,nrow=x,ncol=x)
  
  for( k in 1:x)
    for(l in 1:x)
    {
      extract_mat[k,l] = img[i-kernel+k,j-kernel+l]
    }
  return(extract_mat)
}

######krcenje/erosion#######
krcenje<-function(img,kernel){
  
  img_height = dim(img)[1]
  img_width = dim(img)[2]
  #we  have only odd numbers like 3,5,9...
  
  s=2*kernel+1
  res=matrix(0L, nrow = img_height, ncol = img_width) 

    for(i in (kernel+1):(img_height-kernel-1))
     for(j in (kernel+1):(img_height-kernel-1))
    {
      #extract neighbouring pixel values to calculate erosion
      p = extract_matrix(i,j,img,kernel)
      #find minimum
      res[i,j] = min(p) 
    }
  #return errosion
  return(res)
}
######dilatuion#######
sirjenje<-function(img,kernel){
  
  img_height = dim(img)[1]
  img_width = dim(img)[2]
  #we  have only odd numbers like 3,5,9...
  
  s=2*kernel+1
  res=matrix(0L, nrow = img_height, ncol = img_width) 
  
  for(i in (kernel+1):(img_height-kernel-1))
    for(j in (kernel+1):(img_width-kernel-1))
    {
      #extract neighbouring pixel values to calculate erosion
      p = extract_matrix(i,j,img,kernel)
      #find minimum
      res[i,j] = max(p) 
    }
  return(res)
}

#spectral index MSI  function bands 8 and 11
msi_sentinel<-function(band8,band11,num){
  band8=raster::brick(band8)
  band11=raster::brick(band11)
  msi=overlay(band11, band8, fun=function(x,y){(x/y)})
  output_name= "msi.tiff"
  output_name=paste(num,output_name,sep = "")
  
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
getwd()
setwd("imgs/")

b8_img1<-"T34UDV_20190922T094031_B08_10m.jp2"
b11_img1<-"T34UDV_20190922T094031_B11_20m.jp2"

b8_img2<-"T14UNB_20200910T173911_B08_10m.jp2"
b11_img2<-"T14UMC_20200910T173911_B11_20m.jp2"

b8_img3<-"T10UGA_20200815T185921_B08_10m.jp2"
b11_img3<-"T10UGA_20200815T185921_B11_20m.jp2"


b8_raster3<-raster(readGDAL(b8_img3))
#make bigger resolution of each pixel
b8_raster3_60m= aggregate(b8_raster3, fact=6)
rm(b8_raster3)
b11_raster3<-raster(readGDAL(b11_img3))
#make bigger resolution of each pixel
b11_raster3_60m= aggregate(b11_raster3, fact=3)
rm(b11_raster3)


b8_raster2<-raster(readGDAL(b3_img2))
b8_raster3_20m= aggregate(b8_raster3, fact=2)
b11_raster2<-raster(readGDAL(b11_img2))


###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           DATA INPUT                                ###
###                                                                     ###
###########################################################################
### DESCRIPTION: We get the values of indices and their matrices        ###
###                                                                     ###
###########################################################################

num='3'
msi_3=msi_sentinel(b8_raster3_60m,b11_raster3_60m,num)
rm(b11_raster3_60m)
rm(b8_raster3_60m)
rm(msi_3)
msi_3=raster(readGDAL("3msi.tif"))
msi_3_matrix=raster::as.matrix(msi_3)
rm(msi_3)

###########################################################################
###                                                                     ###
###                           SECTION 3:                                ###
###                           Morphological opening###
###                           and         closing                       ###
###########################################################################
### DESCRIPTION: We find the erosion and dilation of indices            ###
###                                                                     ###
###########################################################################
    
    kernel=3
    ##################morphological opening########################
    msi_3_eroded=krcenje(msi_3_matrix,3)
    rm(msi_3_matrix)
    
    msi_3_dilated=sirjenje(msi_3_eroded,3)
    rm(msi_3_eroded)
    View(msi_3_dilated)
    
    output_name="msi_3_opened.tif"
    raster::writeRaster(msi_3_dilated, filename = output_name, format="GTiff")
    rm(msi_3_dilated)
    
    
    
    ###################################################################
    ##################morphological closing########################
    msi_3_dilated=sirjenje(msi_3_matrix,3)
    rm(msi_3_matrix)
    
    msi_3_eroded=krcenje(msi_3_dilated,3)
    rm(msi_3_dilated)

    output_name="msi_3_closed.tif"
    raster::writeRaster(msi_3_eroded, filename = output_name, format="GTiff")
    rm(msi_3_eroded)

    ###################################################################
    
    