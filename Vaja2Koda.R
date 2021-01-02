
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
erosion<-function(img,kernel){
  
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
dilation<-function(img,kernel){
  
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

dilation_with_reconstruction<-function(img,original_img,kernel){
  img_height = dim(img)[1]
  img_width = dim(img)[2]
  #we  have only odd numbers like 3,5,9...
  
  s=2*kernel+1
  res_b4=original_img
  res=img
  while(isTRUE(all.equal(res_b4,res)))
  {
    res_b4=res
    res=matrix(0L, nrow = img_height, ncol = img_width) 
  for(i in (kernel+1):(img_height-kernel-1))
    for(j in (kernel+1):(img_width-kernel-1))
    {
      #extract neighbouring pixel values to calculate erosion
      p = extract_matrix(i,j,img,kernel)
      #find minimum
      maximum=max(p)
      if(maximum<original_img[i,j])
      {
        res[i,j] =maximum 
      }else
      {
        res[i,j] =original_img[i,j]
      }
      
    }
  }
  return(res)
}
erosion_with_reconstruction<-function(img,original_img,kernel){
  img_height = dim(img)[1]
  img_width = dim(img)[2]
  #we  have only odd numbers like 3,5,9...
  
  s=2*kernel+1
  res_b4=original_img
  res=img
  while(isTRUE(all.equal(res_b4,res)))
  {
    res_b4=res
    res=matrix(0L, nrow = img_height, ncol = img_width) 
    for(i in (kernel+1):(img_height-kernel-1))
      for(j in (kernel+1):(img_width-kernel-1))
      {
        #extract neighbouring pixel values to calculate erosion
        p = extract_matrix(i,j,img,kernel)
        #find minimum
        minimum=min(p)
        if(minimum>original_img[i,j])
        {
          res[i,j] =minimum 
        }else
        {
          res[i,j] =original_img[i,j]
        }
        
      }
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

opening_with_reconstruction<-function(img,kernel){
  
  img_eroded=erosion(img,kernel)
  
  
  img_opened=dilation_with_reconstruction(img_eroded,img,1)
  
  return(img_opened)
}
closing_with_reconstruction<-function(img,kernel){
  img_eroded=dilation(img,kernel)
  
  img_closed=erosion_with_reconstruction(img_eroded,img,1)
  
  return(img_closed)
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

b8_img1<-"T14UMC_20200910T173911_B08_10m.jp2"
b11_img1<-"T14UMC_20200910T173911_B11_20m.jp2"


b8_img2<-"T34UDV_20190922T094031_B08_10m.jp2"
b11_img2<-"T34UDV_20190922T094031_B11_20m.jp2"


b8_img3<-"T10UGA_20200815T185921_B08_10m.jp2"
b11_img3<-"T10UGA_20200815T185921_B11_20m.jp2"


b8_raster1<-raster(readGDAL(b8_img1))
#make bigger resolution of each pixel
b8_raster1_60m= aggregate(b8_raster1, fact=6)
rm(b8_raster1)
b11_raster1<-raster(readGDAL(b11_img1))
b11_raster1_60m= aggregate(b11_raster1, fact=3)
rm(b11_raster1)

b8_raster2<-raster(readGDAL(b8_img2))
b8_raster2_60m= aggregate(b8_raster2, fact=6)
rm(b8_raster2)
b11_raster2<-raster(readGDAL(b11_img2))
b11_raster2_60m= aggregate(b11_raster2, fact=3)
rm(b11_raster2)




b8_raster3<-raster(readGDAL(b8_img3))
#make bigger resolution of each pixel
b8_raster3_60m= aggregate(b8_raster3, fact=6)
rm(b8_raster3)
b11_raster3<-raster(readGDAL(b11_img3))
#make bigger resolution of each pixel
b11_raster3_60m= aggregate(b11_raster3, fact=3)
rm(b11_raster3)




###########################################################################
###                                                                     ###
###                           SECTION 1:                                ###
###                           DATA INPUT                                ###
###                                                                     ###
###########################################################################
### DESCRIPTION: We get the values of indices and their matrices        ###
###                                                                     ###
###########################################################################


num='1'

msi_1=msi_sentinel(b8_raster1_60m,b11_raster1_60m,num)
msi_1_matrix=raster::as.matrix(msi_1)

rm(msi_1)

num='2'
msi_2=msi_sentinel(b8_raster2_60m,b11_raster2_60m,num)
#msi_2=raster(readGDAL("2msi.tif"))
msi_2_matrix=raster::as.matrix(msi_2)
rm(msi_2)
rm(b11_raster2_60m)
rm(b8_raster2_60m)

num='3'
msi_3=msi_sentinel(b8_raster3_60m,b11_raster3_60m,num)
rm(b11_raster1_60m)
rm(b8_raster1_60m)

rm(msi_3)
#msi_3=raster(readGDAL("3msi.tif"))
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
    msi_1_eroded=erosion(msi_1_matrix,3)
#    rm(msi_1_matrix)
    
    msi_1_dilated=dilation(msi_1_eroded,3)
    msi_1_dilated=raster(msi_1_dilated)
    rm(msi_1_eroded)
    View(msi_1_dilated)
    
    output_name="msi_1_opened.tif"
    raster::writeRaster(msi_1_dilated, filename = output_name, format="GTiff")
    rm(msi_1_dilated)
  #2
  
      msi_2_eroded=erosion(msi_2_matrix,3)
    #    rm(msi_1_matrix)
    
    msi_2_dilated=dilation(msi_2_eroded,3)
    msi_2_dilated=raster(msi_2_dilated)
    rm(msi_2_eroded)
    View(msi_2_dilated)
    
    output_name="msi_2_opened.tif"
    raster::writeRaster(msi_2_dilated, filename = output_name, format="GTiff")
    rm(msi_2_dilated)
    

#3

    msi_3_eroded=erosion(msi_3_matrix,3)
    rm(msi_3_matrix)
    
    msi_3_dilated=dilation(msi_3_eroded,3)
    msi_3_dilated=raster(msi_3_dilated)
    rm(msi_3_eroded)
    View(msi_3_dilated)
    
    output_name="msi_3_opened.tif"
    raster::writeRaster(msi_3_dilated, filename = output_name, format="GTiff")
    rm(msi_3_dilated)
    
    
    
    ###################################################################
    ##################morphological closing########################
    #1
    msi_1_dilated=dilation(msi_1_matrix,3)
    rm(msi_1_matrix)
    
    msi_1_eroded=erosion(msi_1_dilated,3)
    msi_1_eroded=raster(msi_1_eroded)
    
    rm(msi_1_dilated)
    
    output_name="msi_1_closed.tif"
    raster::writeRaster(msi_1_eroded, filename = output_name, format="GTiff")
    rm(msi_1_eroded)
    #2
    msi_2_dilated=dilation(msi_2_matrix,3)
    #rm(msi_2_matrix)
    
    msi_2_eroded=erosion(msi_2_dilated,3)
    msi_2_eroded=raster(msi_2_eroded)
    
    rm(msi_2_dilated)
    
    output_name="msi_2_closed.tif"
    raster::writeRaster(msi_2_eroded, filename = output_name, format="GTiff")
    rm(msi_2_eroded)
    
    
    
    #3
    msi_3_dilated=dilation(msi_3_matrix,3)
    rm(msi_3_matrix)
    
    msi_3_eroded=erosion(msi_3_dilated,3)
    msi_3_eroded=raster(msi_3_eroded)
    
    rm(msi_3_dilated)

    output_name="msi_3_closed.tif"
    raster::writeRaster(msi_3_eroded, filename = output_name, format="GTiff")
    rm(msi_3_eroded)
    
    
    

    ###################################################################
    
    ###########################################################################
    ###                                                                     ###
    ###                           SECTION 4:                                ###
    ###                           Morphological opening###
    ###                           and         closing with reconstruction                      ###
    ###########################################################################
    ### DESCRIPTION: We find the erosion and dilation with reconstruction of indices            ###
    ###                                                                     ###
    ###########################################################################
    #1
    msi_1_closed_rec=closing_with_reconstruction(msi_1_matrix,3)
    msi_1_closed_rec=raster(msi_1_closed_rec)
    raster::writeRaster(msi_1_closed_rec, filename = "msi_closed_rec_1", format="GTiff")
    rm(msi_1_closed_rec)
    #opening#
    msi_1_opening_rec=opening_with_reconstruction(msi_1_matrix,3)
    msi_1_opened_rec=raster(msi_1_opening_rec)
    raster::writeRaster(msi_1_opened_rec, filename = "msi_opened_rec_1", format="GTiff")
    rm(msi_1_opening_rec)
    rm(msi_1_opened_rec)
    #2
    msi_2_closed_rec=closing_with_reconstruction(msi_2_matrix,3)
    msi_2_closed_rec=raster(msi_2_closed_rec)
    raster::writeRaster(msi_2_closed_rec, filename = "msi_closed_rec_2", format="GTiff")
    rm(msi_2_closed_rec)
    #opening#
    msi_2_opening_rec=opening_with_reconstruction(msi_2_matrix,3)
    msi_2_opened_rec=raster(msi_2_opening_rec)
    raster::writeRaster(msi_2_opened_rec, filename = "msi_opened_rec_2", format="GTiff")
    rm(msi_2_opened_rec)
    rm(msi_2_opening_rec)
    
    #3
    msi_3_closed_rec=closing_with_reconstruction(msi_3_matrix,3)
    msi_3_losed_rec=raster(msi_3_closed_rec)
    raster::writeRaster(msi_3_losed_rec, filename = "msi_closed_rec_3", format="GTiff")
    rm(msi_3_closed_rec)
    #opening#
    msi_3_opening_rec=opening_with_reconstruction(msi_3_matrix,3)
    msi_3_opened_rec=raster(msi_3_opening_rec)
    raster::writeRaster(msi_3_opened_rec, filename = "msi_opened_rec_3", format="GTiff")
    rm(msi_3_opened_rec)
    rm(msi_3_opening_rec)