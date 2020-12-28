
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
######krcenje/erosion#######
raster2matrix<-function(input_raster){
  img_dim=dim(input_raster)
  #create array from raster
  img_array=as.array(values(input_raster))
  #new matrix
  img_matrix=matrix(as.numeric(0),nrow =img_dim[1],ncol=img_dim[2])
  k=img_dim[1]
  #go from last row to the 1st row w column goin from 1 to last 
  for(i in img_dim[1]:1){
    l=0
    for(j in 1:img_dim[2])
    {
      img_matrix[l,k]=img_array[(img_dim[1]*(i-1))+j]
      #}else
      # {
      #img_matrix[i,j]=img_array[j]
      
      #}
      l=l+1
    }
    k=k-1
  }
  return(img_matrix)
}
#if kernel is 9, than we have a 9x9 matrix. which means we need to start from 4 index 
extract_matrix<-function(i,j,img,kernel){
  #extract neighbouring pixel values to calculate erosion
  x = floor(kernel/2)+1
 # extract =matrix(0L, nrow = kernel, ncol = kernel)
#  extract[3][1]
  extract=matrix(data=0,nrow=kernel,ncol=kernel)
  
  for( k in 1:kernel)
    for(l in 1:kernel)
    {
      print(i-x+k)
      print(j-x+l)
      print(img[i-x+k,j-x+l])
      extract[k,l] = img[i-x+k,j-x+l]
    }
  return(extract)
}

krcenje<-function(img,kernel){
  
  img_height = dim(img)[1]
  img_width = dim(img)[2]
  #we  have only odd numbers like 3,5,9...
  s=floor(kernel/2)+1
  res=matrix(0L, nrow = img_height, ncol = img_height) 
  #for(i in img_height-s:img_height)
   # for(j in img_width-s:img_width)
  for(i in s:img_height-s)
     for(j in s:img_height-s)
    {
      #extract neighbouring pixel values to calculate erosion
      p = extract_matrix(i,j,img,kernel)
      #find minimum
      res[i][j] = min(p) 
    }
  #return errosion
  return(res)
}
######sirjenje/erosion#######

#spectral index MSI  function bands 8 and 11
msi_sentinel<-function(band8,band11){
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
### DESCRIPTION: We get the values of indices and their matrices        ###
###                                                                     ###
###########################################################################

num='3'
msi_3=msi_sentinel(b8_raster3_20m,b11_raster3)
rm(b11_raster3)
rm(b8_raster3_20m)
msi_3
msi_3_matrix=raster2matrix(msi_3)
###########################################################################
###                                                                     ###
###                           SECTION 3:                                ###
###                           Morphological opening###
###                           and         closing                       ###
###########################################################################
### DESCRIPTION: We find the erosion and dilation of indices            ###
###                                                                     ###
###########################################################################


    #distancce between the  central element  and the ones on the edges 
    distance=3
    kernel_len=2*distance+1
    
    extract=extract_matrix(9,9,msi_3_matrix,5)
    msi_3_matrix[9,9]
    msi_3_matrix[9,8]

  msi_skrcen=krcenje(msi_3,3)

