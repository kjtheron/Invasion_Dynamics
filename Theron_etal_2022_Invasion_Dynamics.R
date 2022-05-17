# Raster and shapefile processing
library(rgdal)
library(rgeos)
library(raster)
library(RStoolbox)
library(sf)
library(fasterize)
library(spatialEco)

# Data manipulation and plotting
library(tidyverse)
library(ggpubr)

# Machine learning
library(caret)
library(doParallel)
library(SSDM)

# Stats
library(vcd)
library(ape)
library(ggcorrplot)
library(lme4)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(performance)
library(VSURF)
library(scutr)
library(car)

# Community
library(mvabund)
library(ade4)
library(vegan)

# Set working directory
setwd("~/Desktop/")

###########################################################################
# How wide-spread is bramble, and what is causing higher bramble abundance?
###########################################################################

#### Process satellite imagery ####
# Clip and process Sentinel 2
b2<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B02_20m.jp2")
b3<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B03_20m.jp2")
b4<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B04_20m.jp2")
b8<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B08_10m.jp2")
b8<-projectRaster(b8,b2)
b5<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B05_20m.jp2")
b6<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B06_20m.jp2")
b7<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B07_20m.jp2")
b8a<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B8A_20m.jp2")
b11<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B11_20m.jp2")
b12<-raster("Raster_Data/Sentinel/T36JTN_20200104T074311_B12_20m.jp2")
S2<-stack(b2,b3,b4,b8,b5,b6,b7,b8a,b11,b12)
rm(b2,b3,b4,b8,b5,b6,b7,b8a,b11,b12)
Study_ROI<-readOGR("Shapefiles/Plantations_Clipped.shp")
S2<-crop(S2,Study_ROI)
S2<-mask(S2,Study_ROI)
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
writeRaster(S2,filename="Products/Cliped_S2.tif",format="GTiff",overwrite=TRUE)
rm(S2)

# Clip and process PlanetScope
Estate1_a<-stack("Raster_Data/PlanetScope/Estate1/20200105_074608_1025_3B_AnalyticMS_SR.tif")
Estate1_b<-stack("Raster_Data/PlanetScope/Estate1/20200105_074609_1025_3B_AnalyticMS_SR.tif")
Estate1_c<-stack("Raster_Data/PlanetScope/Estate1/20200105_074610_1025_3B_AnalyticMS_SR.tif")
Estate1<-raster::mosaic(Estate1_a,Estate1_b,Estate1_c,fun=mean)
Estate1<-projectRaster(Estate1,datatype="INT2U",res=3,crs="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Estate1<-crop(Estate1,Study_ROI)
Estate1<-mask(Estate1,Study_ROI)
# Scale blue band
Estate1_blue<-Estate1[[1]]
Estate1_blue<-Estate1_blue-100
Estate1_green<-Estate1[[2]]
Estate1_red<-Estate1[[3]]
Estate1_nir<-Estate1[[4]]
Estate1<-stack(Estate1_blue,Estate1_green,Estate1_red,Estate1_nir)
rm(Estate1_blue,Estate1_green,Estate1_red,Estate1_nir,Estate1_a,Estate1_b,Estate1_c)
Estate2_a<-stack("Raster_Data/PlanetScope/Estate2/20200104_074833_0f3f_3B_AnalyticMS_SR.tif")
Estate2_b<-stack("Raster_Data/PlanetScope/Estate2/20200104_074834_0f3f_3B_AnalyticMS_SR.tif")
Estate2<-raster::mosaic(Estate2_a,Estate2_b,fun=mean)
Estate2<-crop(Estate2,Study_ROI)
Estate2<-mask(Estate2,Study_ROI)
rm(Estate2_a,Estate2_b)
Estate3_a<-stack("Raster_Data/PlanetScope/Estate3/20200104_074828_0f3f_3B_AnalyticMS_SR.tif")
Estate3_b<-stack("Raster_Data/PlanetScope/Estate3/20200104_074829_0f3f_3B_AnalyticMS_SR.tif")
Estate3<-raster::mosaic(Estate3_a,Estate3_b,fun=mean)
Estate3<-crop(Estate3,Study_ROI)
Estate3<-mask(Estate3,Study_ROI)
rm(Estate3_a,Estate3_b)
Estate4_a<-stack("Raster_Data/PlanetScope/Estate4/20200104_074506_0f22_3B_AnalyticMS_SR.tif")
Estate4_b<-stack("Raster_Data/PlanetScope/Estate4/20200104_074507_0f22_3B_AnalyticMS_SR.tif")
Estate4_c<-stack("Raster_Data/PlanetScope/Estate4/20200104_074508_0f22_3B_AnalyticMS_SR.tif")
Estate4<-raster::mosaic(Estate4_a,Estate4_b,Estate4_c,fun=mean)
Estate4<-crop(Estate4,Study_ROI)
Estate4<-mask(Estate4,Study_ROI)
# Scale blue band
Estate4_blue<-Estate4[[1]]
Estate4_blue<-Estate4_blue-30
Estate4_green<-Estate4[[2]]
Estate4_red<-Estate4[[3]]
Estate4_nir<-Estate4[[4]]
Estate4<-stack(Estate4_blue,Estate4_green,Estate4_red,Estate4_nir)
rm(Estate4_blue,Estate4_green,Estate4_red,Estate4_nir,Estate4_a,Estate4_b,Estate4_c)
PS<-merge(Estate2,Estate1,Estate3,Estate4,tolerance=0.2)
names(PS)<-c("Blue","Green","Red","Infrared")
writeRaster(PS,filename="Products/Cliped_PS.tif",format="GTiff",overwrite=TRUE)
rm(Estate1,Estate2,Estate3,Estate4,PS)

# Clip and process super resolution
predict_ms210<-stack("Raster_Data/Super_Resolve/predict_ms210.tif")
predict_ms210<-crop(predict_ms210,Study_ROI)
predict_ms210<-mask(predict_ms210,Study_ROI)
predict_ms220<-stack("Raster_Data/Super_Resolve/predict_ms220.tif")
predict_ms220<-crop(predict_ms220,Study_ROI)
predict_ms220<-mask(predict_ms220,Study_ROI)
SR_S2<-stack(predict_ms210,predict_ms220)
SR_S2<-subset(SR_S2,order(c(3,2,1,4,5,6,7,8,9,10)))
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
writeRaster(SR_S2,filename="Products/Cliped_SR_S2.tif",format="GTiff",overwrite=TRUE)
rm(Study_ROI,predict_ms210,predict_ms220,SR_S2)

#### Create vegetation indices ####
#https://custom-scripts.sentinel-hub.com/sentinel-2/ndvi/
#https://custom-scripts.sentinel-hub.com/sentinel-2/gndvi/
#https://custom-scripts.sentinel-hub.com/sentinel-2/nbr/#
#https://custom-scripts.sentinel-hub.com/sentinel-2/savi/
#https://custom-scripts.sentinel-hub.com/sentinel-2/evi/
#https://custom-scripts.sentinel-hub.com/sentinel-2/red_edge_position/#
#https://custom-scripts.sentinel-hub.com/sentinel-2/msi/#
#https://custom-scripts.sentinel-hub.com/sentinel-2/ndmi/#
#https://custom-scripts.sentinel-hub.com/sentinel-2/ndwi/#
#https://custom-scripts.sentinel-hub.com/sentinel-2/chl_rededge/#

# Sentinel 2
S2<-stack("Products/Cliped_S2.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
NDVI<-(S2[["Infrared"]]-S2[["Red"]])/(S2[["Infrared"]]+S2[["Red"]])
gNDVI<-(S2[["Infrared"]]-S2[["Green"]])/(S2[["Infrared"]]+S2[["Green"]])
NBR<-(S2[["Infrared"]]-S2[["SWIR2"]])/(S2[["Infrared"]]+S2[["SWIR2"]])
SAVI<-(S2[["Infrared"]]-S2[["Red"]])/(S2[["Infrared"]]+S2[["Red"]]+0.725)*(1*0.725)
EVI<-(2.5*(S2[["Infrared"]]-S2[["Red"]])/((S2[["Infrared"]]+(6*S2[["Red"]])-(7.5*S2[["Blue"]]))+1))
REPO<-(700+(40*((((S2[["Red"]]+S2[["RedEdge3"]])/2)-S2[["RedEdge1"]])/(S2[["RedEdge2"]]-S2[["RedEdge1"]]))))
MSI<-S2[["SWIR1"]]/ S2[["RedEdge4"]]
NDMI<-(S2[["Infrared"]]-S2[["SWIR1"]])/(S2[["Infrared"]]+S2[["SWIR1"]])
NDWI<-(S2[["Green"]]-S2[["RedEdge4"]])/(S2[["Green"]]+S2[["RedEdge4"]])
Chl_RE<-(S2[["RedEdge3"]]/S2[["RedEdge1"]])^(-1)
PCA<-rasterPCA(S2,nComp=1)
PCA<-PCA$map
S2<-stack(S2,NDVI,gNDVI,NBR,SAVI,EVI,REPO,MSI,NDMI,NDWI,Chl_RE,PCA)
writeRaster(S2,filename="Products/Cliped_S2_VI.tif",format="GTiff",overwrite=TRUE)
rm(S2,NDVI,gNDVI,NBR,SAVI,EVI,REPO,MSI,NDMI,NDWI,Chl_RE,PCA)

# PlanetScope
PS<-stack("Products/Cliped_PS.tif")
names(PS)<-c("Blue","Green","Red","Infrared")
NDVI<-(PS[["Infrared"]]-PS[["Red"]])/(PS[["Infrared"]]+PS[["Red"]])
gNDVI<-(PS[["Infrared"]]-PS[["Green"]])/(PS[["Infrared"]]+PS[["Green"]])
SAVI<-(PS[["Infrared"]]-PS[["Red"]])/(PS[["Infrared"]]+PS[["Red"]]+0.725)*(1*0.725)
EVI<-(2.5*(PS[["Infrared"]]-PS[["Red"]])/((PS[["Infrared"]]+(6*PS[["Red"]])-(7.5*PS[["Blue"]]))+1))
NDWI<-(PS[["Green"]]-PS[["Infrared"]])/(PS[["Green"]]+PS[["Infrared"]])
PCA<-rasterPCA(PS,nComp=1)
PCA<-PCA$map
PS<-stack(PS,NDVI,gNDVI,SAVI,EVI,NDWI,PCA)
writeRaster(PS,filename="Products/Cliped_PS_VI.tif",format="GTiff",overwrite=TRUE)
rm(PS,NDVI,gNDVI,SAVI,EVI,NDWI,PCA)

# Super-resolved
SR_S2<-stack("Products/Cliped_SR_S2.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
NDVI<-(SR_S2[["Infrared"]]-SR_S2[["Red"]])/(SR_S2[["Infrared"]]+SR_S2[["Red"]])
gNDVI<-(SR_S2[["Infrared"]]-SR_S2[["Green"]])/(SR_S2[["Infrared"]]+SR_S2[["Green"]])
NBR<-(SR_S2[["Infrared"]]-SR_S2[["SWIR2"]])/(SR_S2[["Infrared"]]+SR_S2[["SWIR2"]])
SAVI<-(SR_S2[["Infrared"]]-SR_S2[["Red"]])/(SR_S2[["Infrared"]]+SR_S2[["Red"]]+0.725)*(1*0.725)
EVI<-(2.5*(SR_S2[["Infrared"]]-SR_S2[["Red"]])/((SR_S2[["Infrared"]]+(6*SR_S2[["Red"]])-(7.5*SR_S2[["Blue"]]))+1))
REPO<-(700+(40*((((SR_S2[["Red"]]+SR_S2[["RedEdge3"]])/2)-SR_S2[["RedEdge1"]])/(SR_S2[["RedEdge2"]]-SR_S2[["RedEdge1"]]))))
MSI<-SR_S2[["SWIR1"]]/ SR_S2[["RedEdge4"]]
NDMI<-(SR_S2[["Infrared"]]-SR_S2[["SWIR1"]])/(SR_S2[["Infrared"]]+SR_S2[["SWIR1"]])
NDWI<-(SR_S2[["Green"]]-SR_S2[["RedEdge4"]])/(SR_S2[["Green"]]+SR_S2[["RedEdge4"]])
Chl_RE<-(SR_S2[["RedEdge3"]]/SR_S2[["RedEdge1"]])^(-1)
PCA<-rasterPCA(SR_S2,nComp=1)
PCA<-PCA$map
SR_S2<-stack(SR_S2,NDVI,gNDVI,NBR,SAVI,EVI,REPO,MSI,NDMI,NDWI,Chl_RE,PCA)
writeRaster(SR_S2,filename="Products/Cliped_SR_S2_VI.tif",format="GTiff",overwrite=TRUE)
rm(SR_S2,NDVI,gNDVI,NBR,SAVI,EVI,REPO,MSI,NDMI,NDWI,Chl_RE,PCA)

#### Process training data ####
# Load points
Training_Points<-read.csv("Excel_Sheets/Training_Points.csv")
# Project
Training_Points<-st_as_sf(Training_Points,coords=c("x","y"),crs=4326)
Training_Points<-st_transform(Training_Points,crs=32736)

# Extract training data for Sentinel imagery
S2<-stack("Products/Cliped_S2_VI.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
             "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
             "MSI","NDMI","NDWI","Chl_RE","PCA")
df_S2<-raster::extract(S2,Training_Points,method="simple",df=TRUE)
df_S2$ID=NULL
Class<-as.matrix(Training_Points$Class)
# Save csv file
df_S2<-as.data.frame(cbind(df_S2,"Class"=Class))
df_S2<-df_S2[complete.cases(df_S2),]
write.csv(df_S2,"Products/df_S2_VI.csv",row.names=FALSE)
rm(Class,S2,df_S2)

# Extract training data for PlanetScope imagery
PS<-stack("Products/Cliped_PS_VI.tif")
names(PS)<-c("Blue","Green","Red","Infrared",
             "NDVI","gNDVI","SAVI","EVI","NDWI","PCA")
df_PS<-raster::extract(PS,Training_Points,method="simple",df=TRUE)
df_PS$ID=NULL
Class<-as.matrix(Training_Points$Class)
# Save csv file
df_PS<-as.data.frame(cbind(df_PS,"Class"=Class))
df_PS<-df_PS[complete.cases(df_PS),]
write.csv(df_PS,"Products/df_PS_VI.csv",row.names=FALSE)
rm(Class,PS,df_PS)

# Extract training data for super resolution imagery
SR_S2<-stack("Products/Cliped_SR_S2_VI.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
                "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
                "MSI","NDMI","NDWI","Chl_RE","PCA")
df_SR_S2<-raster::extract(SR_S2,Training_Points,method="simple",df=TRUE)
df_SR_S2$ID=NULL
Class<-as.matrix(Training_Points$Class)
# Save csv file
df_SR_S2<-as.data.frame(cbind(df_SR_S2,"Class"=Class))
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]
write.csv(df_SR_S2,"Products/df_SR_S2_VI.csv",row.names=FALSE)
rm(Class,SR_S2,df_SR_S2,Training_Points)

#### Classification Sentinel 2 ####
# Classify Sentinel 2 using spectral bands only
df_S2<-read.csv("Products/df_S2_VI.csv")
df_S2<-df_S2[complete.cases(df_S2),]
df_S2<-select(df_S2,c(1,2,3,4,5,6,7,8,9,10,22))

# Drop raster variables
S2<-stack("Products/Cliped_S2_VI.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
             "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
             "MSI","NDMI","NDWI","Chl_RE","PCA")
S2<-dropLayer(S2,c(11,12,13,14,15,16,17,18,19,20,21))

# Inspect class balance
plot(df_S2$Class,main="Class Frequencies S2")
table(df_S2$Class)

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_S2$Class,p=0.7,list=FALSE)
df_S2_train<-df_S2[df_Split,]
table(df_S2_train$Class)
df_S2_test<-df_S2[-df_Split,]
table(df_S2_test$Class)
rm(df_Split,df_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_S2_train,"Bramble","Class",70)
Other<-df_S2_train[59:1374,]
df_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2,
                           method="ranger",importance="impurity",
                           data=df_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_S2_Spec_Bands.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_S2_test)
confusionMatrix(data=predict_ranger,df_S2_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_S2_Spec_Bands.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_S2_test,df_S2_train,S2)

####
# Classify Sentinel 2 using spectral bands with variable selection
df_S2<-read.csv("Products/df_S2_VI.csv")
df_S2<-df_S2[complete.cases(df_S2),]
df_S2<-select(df_S2,c(1,2,3,4,5,6,7,8,9,10,22))

# Variable selection
Var_Select<-VSURF(df_S2[,1:10],df_S2[,11])
Var_Select$varselect.pred
df_S2<-select(df_S2,c(1,2,3,5,7,9,11))
rm(Var_Select)

# Drop raster variables
S2<-stack("Products/Cliped_S2_VI.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
             "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
             "MSI","NDMI","NDWI","Chl_RE","PCA")
S2<-dropLayer(S2,c(4,6,8,10,11,12,13,14,15,16,17,18,19,20,21))

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_S2$Class,p=0.7,list=FALSE)
df_S2_train<-df_S2[df_Split,]
table(df_S2_train$Class)
df_S2_test<-df_S2[-df_Split,]
table(df_S2_test$Class)
rm(df_Split,df_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_S2_train,"Bramble","Class",70)
Other<-df_S2_train[59:1374,]
df_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+RedEdge1+RedEdge3+SWIR1,
                           method="ranger",importance="impurity",
                           data=df_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_S2_Spec_Bands_Sel.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_S2_test)
confusionMatrix(data=predict_ranger,df_S2_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_S2_Spec_Bands_Sel.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_S2_test,df_S2_train,S2)

####
# Classify Sentinel 2 using vegetation indices
df_S2<-read.csv("Products/df_S2_VI.csv")
df_S2<-df_S2[complete.cases(df_S2),]

# Load raster variables
S2<-stack("Products/Cliped_S2_VI.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
             "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
             "MSI","NDMI","NDWI","Chl_RE","PCA")

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_S2$Class,p=0.7,list=FALSE)
df_S2_train<-df_S2[df_Split,]
table(df_S2_train$Class)
df_S2_test<-df_S2[-df_Split,]
table(df_S2_test$Class)
rm(df_Split,df_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_S2_train,"Bramble","Class",70)
Other<-df_S2_train[59:1374,]
df_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2+NDVI+gNDVI+NBR+SAVI+EVI+REPO+MSI+NDMI+NDWI+Chl_RE+PCA,
                           method="ranger",importance="impurity",
                           data=df_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_S2_VI.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_S2_test)
confusionMatrix(data=predict_ranger,df_S2_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_S2_VI.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_S2_test,df_S2_train,S2)

####
# Classify Sentinel 2 using vegetation indices and variable selection
df_S2<-read.csv("Products/df_S2_VI.csv")
df_S2<-df_S2[complete.cases(df_S2),]

# Variable selection
Var_Select<-VSURF(df_S2[,1:21],df_S2[,22])
Var_Select$varselect.pred
df_S2<-select(df_S2,c(1,2,3,5,9,12,15,19,22))
rm(Var_Select)

# Drop raster variables
S2<-stack("Products/Cliped_S2_VI.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
             "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
             "MSI","NDMI","NDWI","Chl_RE","PCA")
S2<-dropLayer(S2,c(4,6,7,8,10,11,13,14,16,17,18,20,21))

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_S2$Class,p=0.7,list=FALSE)
df_S2_train<-df_S2[df_Split,]
table(df_S2_train$Class)
df_S2_test<-df_S2[-df_Split,]
table(df_S2_test$Class)
rm(df_Split,df_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_S2_train,"Bramble","Class",70)
Other<-df_S2_train[59:1374,]
df_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+RedEdge1+SWIR1+gNDVI+EVI+NDWI,
                           method="ranger",importance="impurity",
                           data=df_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_S2_VI_Sel.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_S2_test)
confusionMatrix(data=predict_ranger,df_S2_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_S2_VI_Sel.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_S2_test,df_S2_train,S2)

#### Classification PlanetScope ####
# Classify PlanetScope using spectral bands only
df_PS<-read.csv("Products/df_PS_VI.csv")
df_PS<-df_PS[complete.cases(df_PS),]
df_PS<-select(df_PS,c(1,2,3,4,11))

# Drop raster variables
PS<-stack("Products/Cliped_PS_VI.tif")
names(PS)<-c("Blue","Green","Red","Infrared",
             "NDVI","gNDVI","SAVI","EVI","NDWI","PCA")
PS<-dropLayer(PS,c(5,6,7,8,9,10))

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_PS$Class,p=0.7,list=FALSE)
df_PS_train<-df_PS[df_Split,]
table(df_PS_train$Class)
df_PS_test<-df_PS[-df_Split,]
table(df_PS_test$Class)
rm(df_Split,df_PS)

# Oversample Bramble
Bramble<-oversample_smote(df_PS_train,"Bramble","Class",70)
Other<-df_PS_train[59:1376,]
df_PS_train<-as.data.frame(rbind(Bramble,Other))
table(df_PS_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared,
                           method="ranger",importance="impurity",
                           data=df_PS_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_PS_Spec_Bands.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_PS_test)
confusionMatrix(data=predict_ranger,df_PS_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(PS,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_PS_Spec_Bands.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_PS_test,df_PS_train,PS)

####
# Classify PlanetScope using vegetation indices
df_PS<-read.csv("Products/df_PS_VI.csv")
df_PS<-df_PS[complete.cases(df_PS),]

# Load rasters
PS<-stack("Products/Cliped_PS_VI.tif")
names(PS)<-c("Blue","Green","Red","Infrared","NDVI","gNDVI","SAVI","EVI","NDWI","PCA")

# Split dataset
set.seed(691)
df_Split<-createDataPartition(df_PS$Class,p=0.7,list=FALSE)
df_PS_train<-df_PS[df_Split,]
table(df_PS_train$Class)
df_PS_test<-df_PS[-df_Split,]
table(df_PS_test$Class)
rm(df_Split,df_PS)

# Oversample Bramble
Bramble<-oversample_smote(df_PS_train,"Bramble","Class",70)
Other<-df_PS_train[59:1376,]
df_PS_train<-as.data.frame(rbind(Bramble,Other))
table(df_PS_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)  
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+NDVI+gNDVI+SAVI+EVI+NDWI+PCA,
                           method="ranger",importance="impurity",
                           data=df_PS_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_PS_VI.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_PS_test)
confusionMatrix(data=predict_ranger,df_PS_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(PS,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_PS_VI.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_PS_test,df_PS_train,PS)

#### Classification Super-resolution ####
# Classify Super-resolve using spectral bands only
df_SR_S2<-read.csv("Products/df_SR_S2_VI.csv")
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]
df_SR_S2<-select(df_SR_S2,c(1,2,3,4,5,6,7,8,9,10,22))

# Drop raster variables
SR_S2<-stack("Products/Cliped_SR_S2_VI.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
                "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
                "MSI","NDMI","NDWI","Chl_RE","PCA")
SR_S2<-dropLayer(SR_S2,c(11,12,13,14,15,16,17,18,19,20,21))

# Split dataset
set.seed(691)
df_SR_S2_Split<-createDataPartition(df_SR_S2$Class,p=0.7,list=FALSE)
df_SR_S2_train<-df_SR_S2[df_SR_S2_Split,]
table(df_SR_S2_train$Class)
df_SR_S2_test<-df_SR_S2[-df_SR_S2_Split,]
table(df_SR_S2_test$Class)
rm(df_SR_S2_Split,df_SR_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_SR_S2_train,"Bramble","Class",70)
Other<-df_SR_S2_train[59:1376,]
df_SR_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_SR_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2,
                           method="ranger",importance="impurity",
                           data=df_SR_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_SR_S2_Spec_Bands.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_SR_S2_test)
confusionMatrix(data=predict_ranger,df_SR_S2_test$Class)#,mode="prec_recall"
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(SR_S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_SR_S2_Spec_Bands.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_SR_S2_test,df_SR_S2_train,SR_S2)

####
# Classify Super-resolve using spectral bands and variable selection
df_SR_S2<-read.csv("Products/df_SR_S2_VI.csv")
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]
df_SR_S2<-select(df_SR_S2,c(1,2,3,4,5,6,7,8,9,10,22))

# Variable selection
Var_Select<-VSURF(df_SR_S2[,1:10],df_SR_S2[,11])
Var_Select$varselect.pred
df_SR_S2<-select(df_SR_S2,c(1,2,3,4,9,10,11))
rm(Var_Select)

# Drop raster variables
SR_S2<-stack("Products/Cliped_SR_S2_VI.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
                "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
                "MSI","NDMI","NDWI","Chl_RE","PCA")
SR_S2<-dropLayer(SR_S2,c(5,6,7,8,11,12,13,14,15,16,17,18,19,20,21))

# Split dataset
set.seed(691)
df_SR_S2_Split<-createDataPartition(df_SR_S2$Class,p=0.7,list=FALSE)
df_SR_S2_train<-df_SR_S2[df_SR_S2_Split,]
table(df_SR_S2_train$Class)
df_SR_S2_test<-df_SR_S2[-df_SR_S2_Split,]
table(df_SR_S2_test$Class)
rm(df_SR_S2_Split,df_SR_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_SR_S2_train,"Bramble","Class",70)
Other<-df_SR_S2_train[59:1376,]
df_SR_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_SR_S2_train$Class)
rm(Bramble,Other)
 
# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+SWIR1+SWIR2,
                           method="ranger",importance="impurity",
                           data=df_SR_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_SR_S2_Spec_Bands_Sel.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_SR_S2_test)
confusionMatrix(data=predict_ranger,df_SR_S2_test$Class)#,mode="prec_recall"
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(SR_S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_SR_S2_Spec_Bands_Sel.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_SR_S2_test,df_SR_S2_train,SR_S2)

####
# Classify Super-resolve using vegetation indices
df_SR_S2<-read.csv("Products/df_SR_S2_VI.csv")
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]

# Load raster variables
SR_S2<-stack("Products/Cliped_SR_S2_VI.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
                "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
                "MSI","NDMI","NDWI","Chl_RE","PCA")

# Split dataset
set.seed(691)
df_SR_S2_Split<-createDataPartition(df_SR_S2$Class,p=0.7,list=FALSE)
df_SR_S2_train<-df_SR_S2[df_SR_S2_Split,]
table(df_SR_S2_train$Class)
df_SR_S2_test<-df_SR_S2[-df_SR_S2_Split,]
table(df_SR_S2_test$Class)
rm(df_SR_S2_Split,df_SR_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_SR_S2_train,"Bramble","Class",70)
Other<-df_SR_S2_train[59:1376,]
df_SR_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_SR_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2+NDVI+gNDVI+NBR+SAVI+EVI+REPO+MSI+NDMI+NDWI+Chl_RE+PCA,
                           method="ranger",importance="impurity",
                           data=df_SR_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_SR_S2_VI.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_SR_S2_test)
confusionMatrix(data=predict_ranger,df_SR_S2_test$Class)#,mode="prec_recall"
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(SR_S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_SR_S2_VI.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_SR_S2_test,df_SR_S2_train,SR_S2)

####
# Classify Super-resolve using vegetation indices and variable selection
df_SR_S2<-read.csv("Products/df_SR_S2_VI.csv")
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]

# Variable selection
Var_Select<-VSURF(df_SR_S2[,1:21],df_SR_S2[,22])
Var_Select$varselect.pred
df_SR_S2<-select(df_SR_S2,c(2,3,9,10,14,15,19,22))
rm(Var_Select)

# Drop raster variables
SR_S2<-stack("Products/Cliped_SR_S2_VI.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2",
                "NDVI","gNDVI","NBR","SAVI","EVI","REPO",
                "MSI","NDMI","NDWI","Chl_RE","PCA")
SR_S2<-dropLayer(SR_S2,c(1,4,5,6,7,8,11,12,13,16,17,18,20,21))

# Split dataset
set.seed(691)
df_SR_S2_Split<-createDataPartition(df_SR_S2$Class,p=0.7,list=FALSE)
df_SR_S2_train<-df_SR_S2[df_SR_S2_Split,]
table(df_SR_S2_train$Class)
df_SR_S2_test<-df_SR_S2[-df_SR_S2_Split,]
table(df_SR_S2_test$Class)
rm(df_SR_S2_Split,df_SR_S2)

# Oversample Bramble
Bramble<-oversample_smote(df_SR_S2_train,"Bramble","Class",70)
Other<-df_SR_S2_train[59:1376,]
df_SR_S2_train<-as.data.frame(rbind(Bramble,Other))
table(df_SR_S2_train$Class)
rm(Bramble,Other)

# Classification: Random Forest (ranger)
set.seed(873)
model_ranger<-caret::train(Class~Green+Red+SWIR1+SWIR2+SAVI+EVI+NDWI,
                           method="ranger",importance="impurity",
                           data=df_SR_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_SR_S2_VI_Sel.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_SR_S2_test)
confusionMatrix(data=predict_ranger,df_SR_S2_test$Class)
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(SR_S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_SR_S2_VI_Sel.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_SR_S2_test,df_SR_S2_train,SR_S2)

#### Classification plots ####
# Plot variable importance together
model_ranger_SR_S2_VI_Sel<-readRDS("Products/model_ranger_SR_S2_VI_Sel.rds")
model_ranger_SR_S2_VI<-readRDS("Products/model_ranger_SR_S2_VI.rds")
model_ranger_SR_S2_Spec_Bands_Sel<-readRDS("Products/model_ranger_SR_S2_Spec_Bands_Sel.rds")
model_ranger_SR_S2_Spec_Bands<-readRDS("Products/model_ranger_SR_S2_Spec_Bands.rds")
model_ranger_S2_VI_Sel<-readRDS("Products/model_ranger_S2_VI_Sel.rds")
model_ranger_S2_VI<-readRDS("Products/model_ranger_S2_VI.rds")
model_ranger_S2_Spec_Bands_Sel<-readRDS("Products/model_ranger_S2_Spec_Bands_Sel.rds")
model_ranger_S2_Spec_Bands<-readRDS("Products/model_ranger_S2_Spec_Bands.rds")
model_ranger_PS_VI<-readRDS("Products/model_ranger_PS_VI.rds")
model_ranger_PS_Spec_Bands<-readRDS("Products/model_ranger_PS_Spec_Bands.rds")
a<-ggplot(varImp(model_ranger_SR_S2_VI_Sel))+ggtitle("Variable Importance Ranger SR VI selected")+labs(tag="a)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
b<-ggplot(varImp(model_ranger_SR_S2_VI))+ggtitle("Variable Importance Ranger SR VI")+labs(tag="b)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
c<-ggplot(varImp(model_ranger_SR_S2_Spec_Bands_Sel))+ggtitle("Variable Importance Ranger SR bands selected")+labs(tag="c)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
d<-ggplot(varImp(model_ranger_SR_S2_Spec_Bands))+ggtitle("Variable Importance Ranger SR bands")+labs(tag="d)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
e<-ggplot(varImp(model_ranger_S2_VI_Sel))+ggtitle("Variable Importance Ranger S2 VI selected")+labs(tag="e)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
f<-ggplot(varImp(model_ranger_S2_VI))+ggtitle("Variable Importance Ranger S2 VI")+labs(tag="f)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
g<-ggplot(varImp(model_ranger_S2_Spec_Bands_Sel))+ggtitle("Variable Importance Ranger S2 bands selected")+labs(tag="g)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
h<-ggplot(varImp(model_ranger_S2_Spec_Bands))+ggtitle("Variable Importance Ranger S2 bands")+labs(tag="h)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
i<-ggplot(varImp(model_ranger_PS_VI))+ggtitle("Variable Importance Ranger PS VI")+labs(tag="i)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
j<-ggplot(varImp(model_ranger_PS_Spec_Bands))+ggtitle("Variable Importance Ranger PS bands")+labs(tag="j)")+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=13,colour="black"),axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),panel.grid=element_blank(),panel.background=element_blank(),plot.tag=element_text(hjust=0.5))
ggarrange(a,b,c,d,e,f,g,h,i,j,ncol=4,nrow=3)
rm(a,b,c,d,e,f,g,h,i,j,model_ranger_SR_S2_VI_Sel,model_ranger_SR_S2_VI,model_ranger_SR_S2_Spec_Bands_Sel,model_ranger_SR_S2_Spec_Bands,
   model_ranger_S2_VI_Sel,model_ranger_S2_VI,model_ranger_S2_Spec_Bands_Sel,model_ranger_S2_Spec_Bands,model_ranger_PS_VI,model_ranger_PS_Spec_Bands)

# Plot spectral profiles
df<-read.csv("Products/df_SR_S2_VI.csv")
df<-df[complete.cases(df),]
ms<-aggregate(df,list(df$Class),mean)
rownames(ms)<-ms[,1]
ms<-ms[,-1]
ms<-ms[,1:10]
ms<-as.matrix(ms)
mycolor<-c('darkred','yellow','burlywood','cyan','blue','green','gray','magenta')
plot(0,ylim=c(0,4000), xlim=c(1,10),type='n',xlab="Bands",ylab="Reflectance",xaxt="n")
for (i in 1:nrow(ms)){
  lines(ms[i,],type="l",lwd=3,lty=1,col=mycolor[i])
}
title(main="Spectral profiles from Super-resolved imagery",font.main=2)
legend("topleft",rownames(ms),cex=0.8,col=mycolor,lty=1,lwd=3,bty="n")
name<-c("B02","B03","B04","B08","B05","B06","B07","B8A","B11","B12")
axis(1,at=1:10,labels=name)
rm(df,ms,i,mycolor,name)

#### Create bramble response variable ####
# Extract bramble class
Bramble_Class<-raster("Products/Map_ranger_SR_S2_VI.tif")
Bramble_Class[Bramble_Class>1]<-NA
writeRaster(Bramble_Class,filename="Products/Bramble_Class.tif",format="GTiff",overwrite=TRUE)
# Convert to vector
Vector_Bramble_Class<-rasterToPolygons(Bramble_Class)
Vector_Bramble_Class$Obj_ID<-1:nrow(Vector_Bramble_Class)
writeOGR(Vector_Bramble_Class,"Products/Shapefiles/","Vector_Bramble_Class",driver="ESRI Shapefile")
rm(Bramble_Class)
# Dissolve vector
Vector_Bramble_Class_agg<-aggregate(Vector_Bramble_Class,dissolve=TRUE)
Vector_Bramble_Class_dis<-disaggregate(Vector_Bramble_Class_agg)
Vector_Bramble_Class_buf<-buffer(Vector_Bramble_Class_dis,width=0.001,dissolve=TRUE)
Vector_Bramble_Class<-disaggregate(Vector_Bramble_Class_buf)
Vector_Bramble_Class<-as(Vector_Bramble_Class,"SpatialPolygonsDataFrame")
writeOGR(Vector_Bramble_Class,"Products/Shapefiles/","Vector_Bramble_Class_Dis",driver="ESRI Shapefile")
rm(Vector_Bramble_Class_agg,Vector_Bramble_Class_dis,Vector_Bramble_Class_buf)
# Select only large focal areas
Vector_Bramble_Class$Obj_ID<-1:nrow(Vector_Bramble_Class)
Vector_Bramble_Class$Area_sqm<-area(Vector_Bramble_Class)
Vector_Bramble_Class<-st_as_sf(Vector_Bramble_Class)
Vector_Bramble_Class_Filter<-Vector_Bramble_Class %>%
  filter(Area_sqm>=31.25)
st_write(Vector_Bramble_Class_Filter,"Products/Shapefiles/Vector_Bramble_Class_Dis_Filter.shp")

# Generate random points
#https://www.jla-data.net/eng/creating-and-pruning-random-points-and-polygons/
#Vector_Bramble_Class<-st_read("Products/Shapefiles/Vector_Bramble_Class_Dis_Filter.shp")
Vector_Bramble_Class<-st_as_sf(Vector_Bramble_Class_Filter)
set.seed(8534)
Random_Points_Bramble_Class<-st_sample(Vector_Bramble_Class,size=10000,type="random") %>%
  st_sf() %>%
  st_transform(32736)
# Prune to avoid spatial auto-correlation
i<-1
buffer_size<-500
repeat({
  buffer<-st_buffer(Random_Points_Bramble_Class[i,],buffer_size) 
  offending<-Random_Points_Bramble_Class %>%
    st_intersects(buffer,sparse=FALSE)
  offending[i]<-FALSE
  Random_Points_Bramble_Class<-Random_Points_Bramble_Class[!offending,] 
  if(i>=nrow(Random_Points_Bramble_Class)){
    break 
  }else{
    i<-i+1
  }
})
st_write(Random_Points_Bramble_Class,"Products/Shapefiles/Random_Points_Bramble_Class.shp")
rm(Random_Points_Bramble_Class,Vector_Bramble_Class,i,buffer_size,buffer,offending,Vector_Bramble_Class_Filter)

# Clip bramble into two sections
Bramble_Class<-raster("Products/Bramble_Class.tif")
Plantation_Boundaries<-st_read("Shapefiles/Plantation_Boundaries.shp")
Plantation_Boundaries<-st_buffer(Plantation_Boundaries,10)
Production_Areas<-raster::mask(Bramble_Class,Plantation_Boundaries)
writeRaster(Production_Areas,filename="Products/Production_Areas_Bramble.tif",format="GTiff",overwrite=TRUE)
Conservation_Areas<-raster::mask(Bramble_Class,Plantation_Boundaries,inverse=TRUE)
writeRaster(Conservation_Areas,filename="Products/Conservation_Areas_Bramble.tif",format="GTiff",overwrite=TRUE)
rm(Production_Areas,Conservation_Areas)

# Clip points
Random_Points_Bramble_Class<-st_read("Products/Shapefiles/Random_Points_Bramble_Class.shp")
Random_Points_Production_Areas<-Random_Points_Bramble_Class[Plantation_Boundaries,]
st_write(Random_Points_Production_Areas,"Products/Shapefiles/Random_Points_Production_Areas.shp")
Subset<-sapply(st_intersects(Random_Points_Bramble_Class,Plantation_Boundaries),function(x){length(x)==0})
Random_Points_Conservation_Areas<-Random_Points_Bramble_Class[Subset,]
st_write(Random_Points_Conservation_Areas,"Products/Shapefiles/Random_Points_Conservation_Areas.shp")
rm(Random_Points_Bramble_Class,Plantation_Boundaries,Subset)

# Calculate bramble abundance
Production_Areas<-raster("Products/Production_Areas_Bramble.tif")
Conservation_Areas<-raster("Products/Conservation_Areas_Bramble.tif")
Random_Points_Production_Areas_Buf<-st_buffer(Random_Points_Production_Areas,50)
Random_Points_Conservation_Areas_Buf<-st_buffer(Random_Points_Conservation_Areas,50)
Bramble_Production_Abundance<-zonal.stats(Random_Points_Production_Areas_Buf,
                                          Production_Areas,stat=c("sum"))
Bramble_Production_Abundance<-Bramble_Production_Abundance*6.25
Bramble_Conservation_Abundance<-zonal.stats(Random_Points_Conservation_Areas_Buf,
                                            Conservation_Areas,stat=c("sum"))
Bramble_Conservation_Abundance<-Bramble_Conservation_Abundance*6.25
Bramble_Production_Abundance<-data.frame(Bramble_Production_Abundance)
colnames(Bramble_Production_Abundance)<-c("Bram_Abu")
Bramble_Production_Abundance$Section<-c("Production")
Bramble_Conservation_Abundance<-data.frame(Bramble_Conservation_Abundance)
colnames(Bramble_Conservation_Abundance)<-c("Bram_Abu")
Bramble_Conservation_Abundance$Section<-c("Conservation")
Bramble_Abundance<-as.data.frame(rbind(Bramble_Production_Abundance,Bramble_Conservation_Abundance))
write.csv(Bramble_Abundance,"Products/df_Bramble_Response.csv")
rm(Bramble_Abundance,Bramble_Conservation_Abundance,Bramble_Production_Abundance,
   Conservation_Areas,Random_Points_Production_Areas,Random_Points_Production_Areas_Buf,
   Production_Areas,Random_Points_Conservation_Areas,Random_Points_Conservation_Areas_Buf)

#### Create landscape drivers of bramble invasion ####
# Rasterize landscape elements
template_rst<-raster(xmn=200250.2,xmx=261647.7,ymn=6707441,ymx=6766053,resolution=5,
                     crs=projection("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
Roads<-st_read("Shapefiles/Roads_Clipped_Projected.shp")
Roads<-st_zm(Roads)
Roads<-fasterize(Roads,template_rst,field="ID")
Rivers<-st_read("Shapefiles/Rivers_Clipped_Projected.shp")
Rivers<-st_zm(Rivers)
Rivers<-fasterize(Rivers,template_rst,field="ID")

# Calculate distance to elements
Estate1<-readOGR("Shapefiles/Estate1.shp")
Estate2<-readOGR("Shapefiles/Estate2.shp")
Estate3<-readOGR("Shapefiles/Estate3.shp")
Estate4<-readOGR("Shapefiles/Estate4.shp")
Roads_Estate1<-crop(Roads,Estate1)
Roads_Estate1<-mask(Roads_Estate1,Estate1)
Roads_Estate1<-raster::distance(Roads_Estate1,datatype="INT2U")
Roads_Estate2<-crop(Roads,Estate2)
Roads_Estate2<-mask(Roads_Estate2,Estate2)
Roads_Estate2<-raster::distance(Roads_Estate2,datatype="INT2U")
Roads_Estate3<-crop(Roads,Estate3)
Roads_Estate3<-mask(Roads_Estate3,Estate3)
Roads_Estate3<-raster::distance(Roads_Estate3,datatype="INT2U")
Roads_Estate4<-crop(Roads,Estate4)
Roads_Estate4<-mask(Roads_Estate4,Estate4)
Roads_Estate4<-raster::distance(Roads_Estate4,datatype="INT2U")
Roads<-merge(Roads_Estate1,Roads_Estate2,Roads_Estate3,Roads_Estate4)
writeRaster(Roads,filename="Products/Distance_Roads.tif",format="GTiff",overwrite=TRUE)
rm(Roads,Roads_Estate1,Roads_Estate2,Roads_Estate3,Roads_Estate4)
Rivers_Estate1<-crop(Rivers,Estate1)
Rivers_Estate1<-mask(Rivers_Estate1,Estate1)
Rivers_Estate1<-raster::distance(Rivers_Estate1,datatype="INT2U")
Rivers_Estate2<-crop(Rivers,Estate2)
Rivers_Estate2<-mask(Rivers_Estate2,Estate2)
Rivers_Estate2<-raster::distance(Rivers_Estate2,datatype="INT2U")
Rivers_Estate3<-crop(Rivers,Estate3)
Rivers_Estate3<-mask(Rivers_Estate3,Estate3)
Rivers_Estate3<-raster::distance(Rivers_Estate3,datatype="INT2U")
Rivers_Estate4<-crop(Rivers,Estate4)
Rivers_Estate4<-mask(Rivers_Estate4,Estate4)
Rivers_Estate4<-raster::distance(Rivers_Estate4,datatype="INT2U")
Rivers<-merge(Rivers_Estate1,Rivers_Estate2,Rivers_Estate3,Rivers_Estate4)
writeRaster(Rivers,filename="Products/Distance_Rivers.tif",format="GTiff",overwrite=TRUE)
rm(Rivers,Rivers_Estate1,Rivers_Estate2,Rivers_Estate3,Rivers_Estate4)

# Extract distance to woodland areas
Woodland<-raster("Products/Map_ranger_SR_S2_VI.tif")
Woodland<-raster::resample(Woodland,template_rst,method="ngb")
Woodland[Woodland<8]<-NA
Woodland_Estate1<-crop(Woodland,Estate1)
Woodland_Estate1<-mask(Woodland_Estate1,Estate1)
Woodland_Estate1<-raster::distance(Woodland_Estate1,datatype="INT2U")
Woodland_Estate2<-crop(Woodland,Estate2)
Woodland_Estate2<-mask(Woodland_Estate2,Estate2)
Woodland_Estate2<-raster::distance(Woodland_Estate2,datatype="INT2U")
Woodland_Estate3<-crop(Woodland,Estate3)
Woodland_Estate3<-mask(Woodland_Estate3,Estate3)
Woodland_Estate3<-raster::distance(Woodland_Estate3,datatype="INT2U")
Woodland_Estate4<-crop(Woodland,Estate4)
Woodland_Estate4<-mask(Woodland_Estate4,Estate4)
Woodland_Estate4<-raster::distance(Woodland_Estate4,datatype="INT2U")
Woodland<-merge(Woodland_Estate1,Woodland_Estate2,Woodland_Estate3,Woodland_Estate4)
writeRaster(Woodland,filename="Products/Distance_Woodland.tif",format="GTiff",overwrite=TRUE)
rm(Woodland,Woodland_Estate1,Woodland_Estate2,Woodland_Estate3,Woodland_Estate4)

# Calculate distance to plantation boundaries
Plantation_Boundaries<-st_read("Shapefiles/Plantation_Boundaries.shp")
Plantation_Boundaries<-fasterize(Plantation_Boundaries,template_rst)
Plantation_Dist_Estate1<-crop(Plantation_Boundaries,Estate1)
Plantation_Dist_Estate1<-mask(Plantation_Dist_Estate1,Estate1)
Plantation_Dist_Estate1<-raster::distance(Plantation_Dist_Estate1,datatype="INT2U")
Plantation_Dist_Estate2<-crop(Plantation_Boundaries,Estate2)
Plantation_Dist_Estate2<-mask(Plantation_Dist_Estate2,Estate2)
Plantation_Dist_Estate2<-raster::distance(Plantation_Dist_Estate2,datatype="INT2U")
Plantation_Dist_Estate3<-crop(Plantation_Boundaries,Estate3)
Plantation_Dist_Estate3<-mask(Plantation_Dist_Estate3,Estate3)
Plantation_Dist_Estate3<-raster::distance(Plantation_Dist_Estate3,datatype="INT2U")
Plantation_Dist_Estate4<-crop(Plantation_Boundaries,Estate4)
Plantation_Dist_Estate4<-mask(Plantation_Dist_Estate4,Estate4)
Plantation_Dist_Estate4<-raster::distance(Plantation_Dist_Estate4,datatype="INT2U")
Plantation_Dist<-merge(Plantation_Dist_Estate1,Plantation_Dist_Estate2,Plantation_Dist_Estate3,Plantation_Dist_Estate4)
writeRaster(Plantation_Dist,filename="Products/Distance_Plantation.tif",format="GTiff",overwrite=TRUE)
rm(template_rst,Plantation_Dist,Plantation_Boundaries,Plantation_Dist_Estate1,
   Plantation_Dist_Estate2,Plantation_Dist_Estate3,Plantation_Dist_Estate4)

# Calculate distance to plantation harvesting
NBR_2020<-raster("Raster_Data/NBR/Max_NBR_Landsat_2020.tif")
NBR_2019<-raster("Raster_Data/NBR/Max_NBR_Landsat_2019.tif")
NBR_2018<-raster("Raster_Data/NBR/Max_NBR_Landsat_2018.tif")
NBR_2017<-raster("Raster_Data/NBR/Max_NBR_Landsat_2017.tif")
NBR_2016<-raster("Raster_Data/NBR/Max_NBR_Landsat_2016.tif")
Plantation_Boundaries<-st_read("Shapefiles/Plantation_Boundaries.shp")
NBR_2020<-mask(NBR_2020,Plantation_Boundaries)
NBR_2019<-mask(NBR_2019,Plantation_Boundaries)
NBR_2018<-mask(NBR_2018,Plantation_Boundaries)
NBR_2017<-mask(NBR_2017,Plantation_Boundaries)
NBR_2016<-mask(NBR_2016,Plantation_Boundaries)
NBR_2020[NBR_2020>0]<-NA
NBR_2020[NBR_2020<0]<-1
NBR_2019[NBR_2019>0]<-NA
NBR_2019[NBR_2019<0]<-1
NBR_2018[NBR_2018>0]<-NA
NBR_2018[NBR_2018<0]<-1
NBR_2017[NBR_2017>0]<-NA
NBR_2017[NBR_2017<0]<-1
NBR_2016[NBR_2016>0]<-NA
NBR_2016[NBR_2016<0]<-1
NBR_2020<-raster::distance(NBR_2020,datatype="INT2U")
NBR_2019<-raster::distance(NBR_2019,datatype="INT2U")
NBR_2018<-raster::distance(NBR_2018,datatype="INT2U")
NBR_2017<-raster::distance(NBR_2017,datatype="INT2U")
NBR_2016<-raster::distance(NBR_2016,datatype="INT2U")
NBR_2020<-crop(NBR_2020,Plantation_Boundaries)
NBR_2019<-crop(NBR_2019,Plantation_Boundaries)
NBR_2018<-crop(NBR_2018,Plantation_Boundaries)
NBR_2017<-crop(NBR_2017,Plantation_Boundaries)
NBR_2016<-crop(NBR_2016,Plantation_Boundaries)
writeRaster(NBR_2020,filename="Products/Distance_NBR_2020.tif",format="GTiff",overwrite=TRUE)
writeRaster(NBR_2019,filename="Products/Distance_NBR_2019.tif",format="GTiff",overwrite=TRUE)
writeRaster(NBR_2018,filename="Products/Distance_NBR_2018.tif",format="GTiff",overwrite=TRUE)
writeRaster(NBR_2017,filename="Products/Distance_NBR_2017.tif",format="GTiff",overwrite=TRUE)
writeRaster(NBR_2016,filename="Products/Distance_NBR_2016.tif",format="GTiff",overwrite=TRUE)
rm(NBR_2020,NBR_2019,NBR_2018,NBR_2017,NBR_2016,
   Estate1,Estate2,Estate3,Estate4,Plantation_Boundaries)

# Create dataframes
df<-read.csv("Products/df_Bramble_Response.csv")
Bram_Con_df<-subset(df,Section=="Conservation")
Bram_Pro_df<-subset(df,Section=="Production")
rm(df)
Con_Points<-st_read("Products/Shapefiles/Random_Points_Conservation_Areas.shp")
Pro_Points<-st_read("Products/Shapefiles/Random_Points_Production_Areas.shp")
# Load raster
Roads_Dist<-raster("Products/Distance_Roads.tif")
Rivers_Dist<-raster("Products/Distance_Rivers.tif")
Woodland_Dist<-raster("Products/Distance_Woodland.tif")
Plantation_Dist<-raster("Products/Distance_Plantation.tif")
NBR_2020_Dist<-raster("Products/Distance_NBR_2020.tif")
NBR_2019_Dist<-raster("Products/Distance_NBR_2019.tif")
NBR_2018_Dist<-raster("Products/Distance_NBR_2018.tif")
NBR_2017_Dist<-raster("Products/Distance_NBR_2017.tif")
NBR_2016_Dist<-raster("Products/Distance_NBR_2016.tif")
# Extract data over Conservation areas
Con_Roads<-raster::extract(Roads_Dist,Con_Points,method="simple")
Con_Rivers<-raster::extract(Rivers_Dist,Con_Points,method="simple")
Con_Woodland<-raster::extract(Woodland_Dist,Con_Points,method="simple")
Con_Plantation<-raster::extract(Plantation_Dist,Con_Points,method="simple")
Con_NBR_2020<-raster::extract(NBR_2020_Dist,Con_Points,method="simple")
Con_NBR_2019<-raster::extract(NBR_2019_Dist,Con_Points,method="simple")
Con_NBR_2018<-raster::extract(NBR_2018_Dist,Con_Points,method="simple")
Con_NBR_2017<-raster::extract(NBR_2017_Dist,Con_Points,method="simple")
Con_NBR_2016<-raster::extract(NBR_2016_Dist,Con_Points,method="simple")
# Load fire histories from google earth engine
Fire_2020<-read.csv("Excel_Sheets/NBR_2020.csv")
Fire_2020<-Fire_2020[order(Fire_2020$FID,decreasing=FALSE),]
Fire_2020<-Fire_2020$mean
Fire_2019<-read.csv("Excel_Sheets/NBR_2019.csv")
Fire_2019<-Fire_2019[order(Fire_2019$FID,decreasing=FALSE),]
Fire_2019<-Fire_2019$mean
Fire_2018<-read.csv("Excel_Sheets/NBR_2018.csv")
Fire_2018<-Fire_2018[order(Fire_2018$FID,decreasing=FALSE),]
Fire_2018<-Fire_2018$mean
Fire_2017<-read.csv("Excel_Sheets/NBR_2017.csv")
Fire_2017<-Fire_2017[order(Fire_2017$FID,decreasing=FALSE),]
Fire_2017<-Fire_2017$mean
Fire_2016<-read.csv("Excel_Sheets/NBR_2016.csv")
Fire_2016<-Fire_2016[order(Fire_2016$FID,decreasing=FALSE),]
Fire_2016<-Fire_2016$mean

# Construct conservation dataframe
Bram_Con_df<-as.data.frame(cbind(Bram_Con_df,"Con_Roads"=Con_Roads,"Con_Rivers"=Con_Rivers,"Con_Woodland"=Con_Woodland,
                                 "Con_Plantation"=Con_Plantation,"Con_Harvest_2020"=Con_NBR_2020,"Con_Harvest_2019"=Con_NBR_2019,
                                 "Con_Harvest_2018"=Con_NBR_2018,"Con_Harvest_2017"=Con_NBR_2017,"Con_Harvest_2016"=Con_NBR_2016,
                                 "Con_Fire_2020"=Fire_2020,"Con_Fire_2019"=Fire_2019,"Con_Fire_2018"=Fire_2018,
                                 "Con_Fire_2017"=Fire_2017,"Con_Fire_2016"=Fire_2016))
write.csv(Bram_Con_df,"Products/df_Bram_Con.csv")
rm(Bram_Con_df,Con_Roads,Con_Rivers,Con_Woodland,Con_Plantation,Con_NBR_2020,Con_NBR_2019,Con_NBR_2018,Con_NBR_2017,Con_NBR_2016,
   Fire_2020,Fire_2019,Fire_2018,Fire_2017,Fire_2016,Con_Points)
# Extract data over Production areas
Pro_Roads<-raster::extract(Roads_Dist,Pro_Points,method="simple")
Pro_Rivers<-raster::extract(Rivers_Dist,Pro_Points,method="simple")
Pro_Woodland<-raster::extract(Woodland_Dist,Pro_Points,method="simple")
Pro_Plantation<-raster::extract(Plantation_Dist,Pro_Points,method="simple")
Pro_NBR_2020<-raster::extract(NBR_2020_Dist,Pro_Points,method="simple")
Pro_NBR_2019<-raster::extract(NBR_2019_Dist,Pro_Points,method="simple")
Pro_NBR_2018<-raster::extract(NBR_2018_Dist,Pro_Points,method="simple")
Pro_NBR_2017<-raster::extract(NBR_2017_Dist,Pro_Points,method="simple")
Pro_NBR_2016<-raster::extract(NBR_2016_Dist,Pro_Points,method="simple")
# Construct production dataframe
Bram_Pro_df<-as.data.frame(cbind(Bram_Pro_df,"Pro_Roads"=Pro_Roads,"Pro_Rivers"=Pro_Rivers,"Pro_Woodland"=Pro_Woodland,
                                 "Pro_Plantation"=Pro_Plantation,"Pro_Harvest_2020"=Pro_NBR_2020,"Pro_Harvest_2019"=Pro_NBR_2019,
                                 "Pro_Harvest_2018"=Pro_NBR_2018,"Pro_Harvest_2017"=Pro_NBR_2017,"Pro_Harvest_2016"=Pro_NBR_2016))
write.csv(Bram_Pro_df,"Products/df_Bram_Pro.csv")
rm(Bram_Pro_df,Pro_Roads,Pro_Rivers,Pro_Woodland,Pro_Plantation,Pro_NBR_2020,Pro_NBR_2019,Pro_NBR_2018,Pro_NBR_2017,Pro_NBR_2016,
   Pro_Points,NBR_2016_Dist,NBR_2017_Dist,NBR_2018_Dist,NBR_2019_Dist,NBR_2020_Dist,Plantation_Dist,Roads_Dist,Woodland_Dist,Rivers_Dist)

#### Modeling drivers of bramble invasion ####
df_Bram_Con<-read.csv("Products/df_Bram_Con.csv")
df_Bram_Con<-df_Bram_Con[,3:18]
df_Bram_Con$Bram_Abu<-round(df_Bram_Con$Bram_Abu)
df_Bram_Pro<-read.csv("Products/df_Bram_Pro.csv")
df_Bram_Pro<-df_Bram_Pro[,3:13]
df_Bram_Pro$Bram_Abu<-round(df_Bram_Pro$Bram_Abu)

# Normality
hist(df_Bram_Con$Bram_Abu)
qqnorm(df_Bram_Con$Bram_Abu)
qqline(df_Bram_Con$Bram_Abu,lty=2)
shapiro.test(df_Bram_Con$Bram_Abu)
summary(goodfit(df_Bram_Con$Bram_Abu,type="nbinomial",method="ML"))
hist(df_Bram_Pro$Bram_Abu)
qqnorm(df_Bram_Pro$Bram_Abu)
qqline(df_Bram_Pro$Bram_Abu,lty=2)
shapiro.test(df_Bram_Pro$Bram_Abu)
summary(goodfit(df_Bram_Pro$Bram_Abu,type="nbinomial",method="ML"))

#Transforming response
df_Bram_Con$Bram_Abu_Sq<-sqrt(df_Bram_Con$Bram_Abu)
df_Bram_Con$Bram_Abu_Sq<-round(df_Bram_Con$Bram_Abu_Sq)
df_Bram_Pro$Bram_Abu_Sq<-sqrt(df_Bram_Pro$Bram_Abu)
df_Bram_Pro$Bram_Abu_Sq<-round(df_Bram_Pro$Bram_Abu_Sq)
hist(df_Bram_Con$Bram_Abu_Sq)
qqnorm(df_Bram_Con$Bram_Abu_Sq)
qqline(df_Bram_Con$Bram_Abu_Sq,lty=2)
shapiro.test(df_Bram_Con$Bram_Abu_Sq)
summary(goodfit(df_Bram_Con$Bram_Abu_Sq,type="nbinomial",method="ML"))
hist(df_Bram_Pro$Bram_Abu_Sq)
qqnorm(df_Bram_Pro$Bram_Abu_Sq)
qqline(df_Bram_Pro$Bram_Abu_Sq,lty=2)
shapiro.test(df_Bram_Pro$Bram_Abu_Sq)
summary(goodfit(df_Bram_Pro$Bram_Abu_Sq,type="nbinomial",method="ML"))
df_Bram_Con$Section=NULL
df_Bram_Pro$Section=NULL

# Check for correlation between features
features<-df_Bram_Con[,2:15]
features$Con_Harvest_2019=NULL
features$Con_Harvest_2018=NULL
features$Con_Harvest_2017=NULL
features$Con_Harvest_2016=NULL
features$Con_Fire_2019=NULL
features$Con_Fire_2018=NULL
features$Con_Fire_2017=NULL
features$Con_Fire_2016=NULL
ggcorrplot(cor(features,method="spearman"),
           type="lower",lab=TRUE,lab_size=3,tl.cex= 10)
features<-df_Bram_Pro[,2:10]
features$Pro_Harvest_2019=NULL
features$Pro_Harvest_2018=NULL
features$Pro_Harvest_2017=NULL
features$Pro_Harvest_2016=NULL
ggcorrplot(cor(features,method="spearman"),
           type="lower",lab=TRUE,lab_size=3,tl.cex= 10)
rm(features)

# Spatial autocorrelation
Con_Points_GPS<-st_read("Products/Shapefiles/Random_Points_Conservation_Areas.shp")
Pro_Points_GPS<-st_read("Products/Shapefiles/Random_Points_Production_Areas.shp")
Con_Points_GPS<-as.data.frame(st_coordinates(Con_Points_GPS))
Pro_Points_GPS<-as.data.frame(st_coordinates(Pro_Points_GPS))
Data.dist.inv_Con<-1/(as.matrix(dist(cbind(Con_Points_GPS$X,Con_Points_GPS$Y))))
Data.dist.inv_Pro<-1/(as.matrix(dist(cbind(Pro_Points_GPS$X,Pro_Points_GPS$Y))))
Data.dist.inv_Con[is.infinite(Data.dist.inv_Con)]<-0
Data.dist.inv_Pro[is.infinite(Data.dist.inv_Pro)]<-0
Moran.I(df_Bram_Con$Bram_Abu,Data.dist.inv_Con)
Moran.I(df_Bram_Pro$Bram_Abu,Data.dist.inv_Pro)
Moran.I(df_Bram_Con$Bram_Abu_Sq,Data.dist.inv_Con)
Moran.I(df_Bram_Pro$Bram_Abu_Sq,Data.dist.inv_Pro)
rm(Con_Points_GPS,Pro_Points_GPS,Data.dist.inv_Con,Data.dist.inv_Pro)

# Create random variable
Con_Points_GPS<-readOGR("Products/Shapefiles/Random_Points_Conservation_Areas.shp")
Pro_Points_GPS<-readOGR("Products/Shapefiles/Random_Points_Production_Areas.shp")
Plantation<-readOGR("Shapefiles/Plantations_Clipped.shp")
Plantation_Con<-over(Con_Points_GPS,Plantation)
Plantation_Pro<-over(Pro_Points_GPS,Plantation)
Plantation_Con$WPU_NAME=NULL
Plantation_Pro$WPU_NAME=NULL
rm(Con_Points_GPS,Pro_Points_GPS,Plantation)

# Scaling explanatory variables
df_Bram_Con_Bram_Abu<-df_Bram_Con$Bram_Abu
df_Bram_Con_Bram_Abu_Sq<-df_Bram_Con$Bram_Abu_Sq
df_Bram_Pro_Bram_Abu<-df_Bram_Pro$Bram_Abu
df_Bram_Pro_Bram_Abu_Sq<-df_Bram_Pro$Bram_Abu_Sq
df_Bram_Con_scale<-as.data.frame(scale(df_Bram_Con))
df_Bram_Con_scale$Bram_Abu=NULL
df_Bram_Con_scale$Bram_Abu_Sq=NULL
df_Bram_Con_scale<-as.data.frame(cbind("Bram_Abu"=df_Bram_Con_Bram_Abu,"Bram_Abu_Sq"=df_Bram_Con_Bram_Abu_Sq,df_Bram_Con_scale,Plantation_Con))
df_Bram_Pro_scale<-as.data.frame(scale(df_Bram_Pro))
df_Bram_Pro_scale$Bram_Abu=NULL
df_Bram_Pro_scale$Bram_Abu_Sq=NULL
df_Bram_Pro_scale<-as.data.frame(cbind("Bram_Abu"=df_Bram_Pro_Bram_Abu,"Bram_Abu_Sq"=df_Bram_Pro_Bram_Abu_Sq,df_Bram_Pro_scale,Plantation_Pro))
df_Bram_Con<-as.data.frame(cbind(df_Bram_Con,Plantation_Con))
df_Bram_Pro<-as.data.frame(cbind(df_Bram_Pro,Plantation_Pro))
rm(df_Bram_Con_Bram_Abu,df_Bram_Con_Bram_Abu_Sq,df_Bram_Pro_Bram_Abu,df_Bram_Pro_Bram_Abu_Sq,Plantation_Con,Plantation_Pro)

# Add GPS data and remove NAs
Con_Points_GPS<-st_read("Products/Shapefiles/Random_Points_Conservation_Areas.shp")
Pro_Points_GPS<-st_read("Products/Shapefiles/Random_Points_Production_Areas.shp")
Con_Points_GPS<-as.data.frame(st_coordinates(Con_Points_GPS))
Pro_Points_GPS<-as.data.frame(st_coordinates(Pro_Points_GPS))
df_Bram_Con<-as.data.frame(cbind(df_Bram_Con,Con_Points_GPS))
df_Bram_Pro<-as.data.frame(cbind(df_Bram_Pro,Pro_Points_GPS))
df_Bram_Con_scale<-as.data.frame(cbind(df_Bram_Con_scale,Con_Points_GPS))
df_Bram_Pro_scale<-as.data.frame(cbind(df_Bram_Pro_scale,Pro_Points_GPS))
df_Bram_Con<-df_Bram_Con[complete.cases(df_Bram_Con),]
df_Bram_Con_scale<-df_Bram_Con_scale[complete.cases(df_Bram_Con_scale),]
df_Bram_Pro<-df_Bram_Pro[complete.cases(df_Bram_Pro),]
df_Bram_Pro_scale<-df_Bram_Pro_scale[complete.cases(df_Bram_Pro_scale),]
rm(Con_Points_GPS,Pro_Points_GPS)

# Outliers
df_Bram_Con<-df_Bram_Con[-c(32,279),]
df_Bram_Con_scale<-df_Bram_Con_scale[-c(32,279),]

# Model building: Conservation areas
# 2020
Model_Con_2020_nb<-glmmTMB(Bram_Abu_Sq~Con_Rivers+Con_Woodland+Con_Harvest_2020+Con_Fire_2020+
                             Con_Rivers:Con_Fire_2020+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2020_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.26         1.12      0.79
             Con_Woodland 1.14         1.07      0.88
         Con_Harvest_2020 1.07         1.03      0.94
            Con_Fire_2020 1.24         1.12      0.80
 Con_Rivers:Con_Fire_2020 1.22         1.11      0.82
'''
options(na.action="na.fail")
Model_Con_2020_nb_Dredge<-dredge(Model_Con_2020_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2020_nb_Subset<-subset(Model_Con_2020_nb_Dredge,delta<2)
Model_Con_2020_nb_Ave<-model.avg(Model_Con_2020_nb_Subset)
sw(Model_Con_2020_nb_Ave)
'''
                     cond(Con_Fire_2020) cond(Con_Rivers) cond(Con_Harvest_2020) cond(Con_Fire_2020:Con_Rivers)
Sum of weights:      1.00                1.00             0.32                   0.20                          
N containing models:    3                   3                1                      1
'''
confint(Model_Con_2020_nb_Ave)
'''
                                     2.5 %       97.5 %
cond((Int))                     3.45178990  3.581718441
cond(Con_Fire_2020)            -0.15364749 -0.054664679 *
cond(Con_Rivers)               -0.11050916 -0.006687638 *
cond(Con_Harvest_2020)         -0.01923586  0.070295133
cond(Con_Fire_2020:Con_Rivers) -0.06624015  0.038287990
'''
summary(Model_Con_2020_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.516754   0.033004    0.033146 106.100  < 2e-16 ***
cond(Con_Fire_2020)            -0.104156   0.025145    0.025251   4.125 3.71e-05 ***
cond(Con_Rivers)               -0.058598   0.026373    0.026486   2.212   0.0269 *  
cond(Con_Harvest_2020)          0.008199   0.017556    0.017597   0.466   0.6412    
cond(Con_Fire_2020:Con_Rivers) -0.002734   0.012987    0.013032   0.210   0.8338    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.51675    0.03300     0.03315 106.100  < 2e-16 ***
cond(Con_Fire_2020)            -0.10416    0.02514     0.02525   4.125 3.71e-05 ***
cond(Con_Rivers)               -0.05860    0.02637     0.02649   2.212   0.0269 *  
cond(Con_Harvest_2020)          0.02553    0.02274     0.02284   1.118   0.2637    
cond(Con_Fire_2020:Con_Rivers) -0.01398    0.02655     0.02667   0.524   0.6002
'''
# Model assumption
plot(simulateResiduals(Model_Con_2020_nb))
testDispersion(simulateResiduals(Model_Con_2020_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Con_2020_nb),x=df_Bram_Con$X,y=df_Bram_Con$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Con_2020_nb),group=df_Bram_Con$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Con_2020_nb,Model_Con_2020_nb_Dredge,Model_Con_2020_nb_Subset,
   Model_Con_2020_nb_Ave,res_scaled,Plantation,Loc)

#2019
Model_Con_2019_nb<-glmmTMB(Bram_Abu_Sq~Con_Rivers+Con_Woodland+Con_Harvest_2019+Con_Fire_2019+
                             Con_Rivers:Con_Fire_2019+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2019_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.03         1.01      0.97
             Con_Woodland 1.16         1.08      0.86
         Con_Harvest_2019 1.08         1.04      0.92
            Con_Fire_2019 1.21         1.10      0.83
 Con_Rivers:Con_Fire_2019 1.03         1.02      0.97
'''
options(na.action="na.fail")
Model_Con_2019_nb_Dredge<-dredge(Model_Con_2019_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2019_nb_Subset<-subset(Model_Con_2019_nb_Dredge,delta<2)
Model_Con_2019_nb_Ave<-model.avg(Model_Con_2019_nb_Subset)
sw(Model_Con_2019_nb_Ave)
'''
                     cond(Con_Fire_2019) cond(Con_Harvest_2019) cond(Con_Rivers) cond(Con_Fire_2019:Con_Rivers)
Sum of weights:      1.00                1.00                   0.81             0.28                          
N containing models:    3                   3                      2                1  
'''
confint(Model_Con_2019_nb_Ave)
'''
                                     2.5 %       97.5 %
cond((Int))                     3.46608883  3.568112422
cond(Con_Fire_2019)            -0.16321642 -0.067472062 *
cond(Con_Harvest_2019)         -0.09986252 -0.007673383 *
cond(Con_Rivers)               -0.09911270 -0.001196847 *
cond(Con_Fire_2019:Con_Rivers) -0.02749491  0.075936662
'''
summary(Model_Con_2019_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.517101   0.025916    0.026027 135.133  < 2e-16 ***
cond(Con_Fire_2019)            -0.115344   0.024322    0.024425   4.722  2.3e-06 ***
cond(Con_Harvest_2019)         -0.053768   0.023418    0.023518   2.286   0.0222 *  
cond(Con_Rivers)               -0.040398   0.029874    0.029945   1.349   0.1773    
cond(Con_Fire_2019:Con_Rivers)  0.006794   0.017664    0.017711   0.384   0.7013    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.51710    0.02592     0.02603 135.133  < 2e-16 ***
cond(Con_Fire_2019)            -0.11534    0.02432     0.02443   4.722  2.3e-06 ***
cond(Con_Harvest_2019)         -0.05377    0.02342     0.02352   2.286   0.0222 *  
cond(Con_Rivers)               -0.05015    0.02487     0.02498   2.008   0.0447 *  
cond(Con_Fire_2019:Con_Rivers)  0.02422    0.02627     0.02639   0.918   0.3586
'''
# Model assumption
plot(simulateResiduals(Model_Con_2019_nb))
testDispersion(simulateResiduals(Model_Con_2019_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Con_2019_nb),x=df_Bram_Con$X,y=df_Bram_Con$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Con_2019_nb),group=df_Bram_Con$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Con_2019_nb,Model_Con_2019_nb_Dredge,Model_Con_2019_nb_Subset,
   Model_Con_2019_nb_Ave,res_scaled,Plantation,Loc)

#2018
Model_Con_2018_nb<-glmmTMB(Bram_Abu_Sq~Con_Rivers+Con_Woodland+Con_Harvest_2018+Con_Fire_2018+
                             Con_Rivers:Con_Fire_2018+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2018_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.13         1.06      0.88
             Con_Woodland 1.14         1.07      0.88
         Con_Harvest_2018 1.10         1.05      0.91
            Con_Fire_2018 1.15         1.07      0.87
 Con_Rivers:Con_Fire_2018 1.10         1.05      0.91
'''
options(na.action="na.fail")
Model_Con_2018_nb_Dredge<-dredge(Model_Con_2018_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2018_nb_Subset<-subset(Model_Con_2018_nb_Dredge,delta<2)
Model_Con_2018_nb_Ave<-model.avg(Model_Con_2018_nb_Subset)
sw(Model_Con_2018_nb_Ave)
'''
                     cond(Con_Fire_2018) cond(Con_Rivers) cond(Con_Harvest_2018) cond(Con_Fire_2018:Con_Rivers)
Sum of weights:      1.00                1.00             0.30                   0.21                          
N containing models:    3                   3                1                      1  
'''
confint(Model_Con_2018_nb_Ave)
'''
                                     2.5 %        97.5 %
cond((Int))                     3.45414257  3.5823208331
cond(Con_Fire_2018)            -0.17321303 -0.0770990946 *
cond(Con_Rivers)               -0.10096458 -0.0009751924 *
cond(Con_Harvest_2018)         -0.07053335  0.0207921728
cond(Con_Fire_2018:Con_Rivers) -0.07181300  0.0373116865
'''
summary(Model_Con_2018_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.518232   0.032560    0.032699 107.594   <2e-16 ***
cond(Con_Fire_2018)            -0.125156   0.024415    0.024519   5.104    3e-07 ***
cond(Con_Rivers)               -0.050970   0.025401    0.025508   1.998   0.0457 *  
cond(Con_Harvest_2018)         -0.007553   0.017153    0.017194   0.439   0.6605    
cond(Con_Fire_2018:Con_Rivers) -0.003598   0.014471    0.014518   0.248   0.8043    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.51823    0.03256     0.03270 107.594   <2e-16 ***
cond(Con_Fire_2018)            -0.12516    0.02442     0.02452   5.104    3e-07 ***
cond(Con_Rivers)               -0.05097    0.02540     0.02551   1.998   0.0457 *  
cond(Con_Harvest_2018)         -0.02487    0.02320     0.02330   1.068   0.2857    
cond(Con_Fire_2018:Con_Rivers) -0.01725    0.02772     0.02784   0.620   0.5355
'''
# Model assumption
plot(simulateResiduals(Model_Con_2018_nb))
testDispersion(simulateResiduals(Model_Con_2018_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Con_2018_nb),x=df_Bram_Con$X,y=df_Bram_Con$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Con_2018_nb),group=df_Bram_Con$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Con_2018_nb,Model_Con_2018_nb_Dredge,Model_Con_2018_nb_Subset,
   Model_Con_2018_nb_Ave,res_scaled,Plantation,Loc)

#2017
Model_Con_2017_nb<-glmmTMB(Bram_Abu_Sq~Con_Rivers+Con_Woodland+Con_Harvest_2017+Con_Fire_2017+
                             Con_Rivers:Con_Fire_2017+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2017_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.07         1.04      0.93
             Con_Woodland 1.12         1.06      0.90
         Con_Harvest_2017 1.17         1.08      0.85
            Con_Fire_2017 1.12         1.06      0.89
 Con_Rivers:Con_Fire_2017 1.04         1.02      0.96
'''
options(na.action="na.fail")
Model_Con_2017_nb_Dredge<-dredge(Model_Con_2017_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2017_nb_Subset<-subset(Model_Con_2017_nb_Dredge,delta<2)
Model_Con_2017_nb_Ave<-model.avg(Model_Con_2017_nb_Subset)
sw(Model_Con_2017_nb_Ave)
'''
                     cond(Con_Fire_2017) cond(Con_Rivers) cond(Con_Harvest_2017) cond(Con_Fire_2017:Con_Rivers)
Sum of weights:      1.00                0.76             0.17                   0.16                          
N containing models:    4                   3                1                      1   
'''
confint(Model_Con_2017_nb_Ave)
'''
                                     2.5 %       97.5 %
cond((Int))                     3.45339922  3.588448613
cond(Con_Fire_2017)            -0.20648096 -0.106408534 *
cond(Con_Rivers)               -0.09195442  0.004174907
cond(Con_Harvest_2017)         -0.05466921  0.036044872
cond(Con_Fire_2017:Con_Rivers) -0.04531160  0.061692779
'''
summary(Model_Con_2017_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.520924   0.034305    0.034452 102.198   <2e-16 ***
cond(Con_Fire_2017)            -0.156445   0.025421    0.025529   6.128   <2e-16 ***
cond(Con_Rivers)               -0.033562   0.028330    0.028399   1.182    0.237    
cond(Con_Harvest_2017)         -0.001555   0.010035    0.010073   0.154    0.877    
cond(Con_Fire_2017:Con_Rivers)  0.001319   0.011316    0.011361   0.116    0.908    
 
(conditional average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     3.520924   0.034305    0.034452 102.198   <2e-16 ***
cond(Con_Fire_2017)            -0.156445   0.025421    0.025529   6.128   <2e-16 ***
cond(Con_Rivers)               -0.043890   0.024419    0.024523   1.790   0.0735 .  
cond(Con_Harvest_2017)         -0.009312   0.023043    0.023142   0.402   0.6874    
cond(Con_Fire_2017:Con_Rivers)  0.008191   0.027181    0.027298   0.300   0.7641
'''
# Model assumption
plot(simulateResiduals(Model_Con_2017_nb))
testDispersion(simulateResiduals(Model_Con_2017_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Con_2017_nb),x=df_Bram_Con$X,y=df_Bram_Con$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Con_2017_nb),group=df_Bram_Con$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Con_2017_nb,Model_Con_2017_nb_Dredge,Model_Con_2017_nb_Subset,
   Model_Con_2017_nb_Ave,res_scaled,Plantation,Loc)

#2016
Model_Con_2016_nb<-glmmTMB(Bram_Abu_Sq~Con_Rivers+Con_Woodland+Con_Harvest_2016+Con_Fire_2016+
                             Con_Rivers:Con_Fire_2016+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2016_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.10         1.05      0.91
             Con_Woodland 1.11         1.05      0.90
         Con_Harvest_2016 1.13         1.06      0.89
            Con_Fire_2016 1.13         1.06      0.88
 Con_Rivers:Con_Fire_2016 1.04         1.02      0.96
'''
options(na.action="na.fail")
Model_Con_2016_nb_Dredge<-dredge(Model_Con_2016_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2016_nb_Subset<-subset(Model_Con_2016_nb_Dredge,delta<2)
Model_Con_2016_nb_Ave<-model.avg(Model_Con_2016_nb_Subset)
sw(Model_Con_2016_nb_Ave)
'''
                     cond(Con_Fire_2016) cond(Con_Rivers) cond(Con_Woodland)
Sum of weights:      1.00                0.72             0.19              
N containing models:    3                   2                1
'''
confint(Model_Con_2016_nb_Ave)
'''
                          2.5 %       97.5 %
cond((Int))          3.45233143  3.583781556
cond(Con_Fire_2016) -0.20559829 -0.111696650 *
cond(Con_Rivers)    -0.09113861  0.003663894
cond(Con_Woodland)  -0.05320036  0.037810673
'''
summary(Model_Con_2016_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                     Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))          3.518056   0.033391    0.033534 104.911   <2e-16 ***
cond(Con_Fire_2016) -0.158647   0.023853    0.023955   6.623   <2e-16 ***
cond(Con_Rivers)    -0.031292   0.028361    0.028424   1.101    0.271    
cond(Con_Woodland)  -0.001494   0.010631    0.010673   0.140    0.889    
 
(conditional average) 
                     Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))          3.518056   0.033391    0.033534 104.911   <2e-16 ***
cond(Con_Fire_2016) -0.158647   0.023853    0.023955   6.623   <2e-16 ***
cond(Con_Rivers)    -0.043737   0.024082    0.024185   1.808   0.0705 .  
cond(Con_Woodland)  -0.007695   0.023118    0.023218   0.331   0.7403
'''
# Model assumption
plot(simulateResiduals(Model_Con_2016_nb))
testDispersion(simulateResiduals(Model_Con_2016_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Con_2016_nb),x=df_Bram_Con$X,y=df_Bram_Con$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Con_2016_nb),group=df_Bram_Con$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Con_2016_nb,Model_Con_2016_nb_Dredge,Model_Con_2016_nb_Subset,
   Model_Con_2016_nb_Ave,res_scaled,Plantation,Loc)

# Plot significant variables together
a<-ggplot(data=df_Bram_Con,aes(x=Con_Rivers,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nDistance to rivers (m)",tag="a)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
b<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2020,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nFire severity 2020",tag="b)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
c<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2019,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nFire severity 2019",tag="c)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
d<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2018,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nFire severity 2018",tag="d)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
e<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2017,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nFire severity 2017",tag="e)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
f<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2016,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nFire severity 2016",tag="f)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
g<-ggplot(data=df_Bram_Con,aes(x=Con_Harvest_2019,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands (squared)\n",x="\nDistance to harvested plantations 2019",tag="g)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
ggarrange(a,b,c,d,e,f,g,ncol=3,nrow=3) # 14 x 8.5
rm(a,b,c,d,e,f,g,df_Bram_Con,df_Bram_Con_scale)

# Model building: Production areas
# 2020
Model_Pro_2020_nb<-glmmTMB(Bram_Abu_Sq~Pro_Rivers+Pro_Woodland+Pro_Harvest_2020+
                             Pro_Rivers:Pro_Harvest_2020+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2020_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.03         1.02      0.97
                Pro_Woodland 1.16         1.08      0.86
            Pro_Harvest_2020 1.13         1.06      0.88
 Pro_Rivers:Pro_Harvest_2020 1.03         1.01      0.97
'''
options(na.action="na.fail")
Model_Pro_2020_nb_Dredge<-dredge(Model_Pro_2020_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2020_nb_Subset<-subset(Model_Pro_2020_nb_Dredge,delta<2)
Model_Pro_2020_nb_Ave<-model.avg(Model_Pro_2020_nb_Subset)
sw(Model_Pro_2020_nb_Ave)
'''
                     cond(Pro_Harvest_2020) cond(Pro_Rivers) cond(Pro_Woodland)
Sum of weights:      1.00                   0.62             0.20              
N containing models:    3                      2                1 
'''
confint(Model_Pro_2020_nb_Ave)
'''
                             2.5 %      97.5 %
cond((Int))             3.17651364  3.36113497
cond(Pro_Harvest_2020) -0.17670340 -0.04872737 *
cond(Pro_Rivers)       -0.11133548  0.01282449
cond(Pro_Woodland)     -0.04124414  0.09394682
'''
summary(Model_Pro_2020_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                        Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.268824   0.046837    0.047098  69.405  < 2e-16 ***
cond(Pro_Harvest_2020) -0.112715   0.032469    0.032648   3.452 0.000555 ***
cond(Pro_Rivers)       -0.030554   0.034451    0.034550   0.884 0.376517    
cond(Pro_Woodland)      0.005223   0.018534    0.018604   0.281 0.778901    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.26882    0.04684     0.04710  69.405  < 2e-16 ***
cond(Pro_Harvest_2020) -0.11272    0.03247     0.03265   3.452 0.000555 ***
cond(Pro_Rivers)       -0.04926    0.03150     0.03167   1.555 0.119928    
cond(Pro_Woodland)      0.02635    0.03430     0.03449   0.764 0.444825 
'''
# Model assumption
plot(simulateResiduals(Model_Pro_2020_nb))
testDispersion(simulateResiduals(Model_Pro_2020_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Pro_2020_nb),x=df_Bram_Pro$X,y=df_Bram_Pro$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Pro_2020_nb),group=df_Bram_Pro$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Pro_2020_nb,Model_Pro_2020_nb_Dredge,Model_Pro_2020_nb_Subset,
   Model_Pro_2020_nb_Ave,res_scaled,Plantation,Loc)

#2019
Model_Pro_2019_nb<-glmmTMB(Bram_Abu_Sq~Pro_Rivers+Pro_Woodland+Pro_Harvest_2019+
                             Pro_Rivers:Pro_Harvest_2019+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2019_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.04         1.02      0.96
                Pro_Woodland 1.07         1.04      0.93
            Pro_Harvest_2019 1.06         1.03      0.94
 Pro_Rivers:Pro_Harvest_2019 1.03         1.01      0.97
'''
options(na.action="na.fail")
Model_Pro_2019_nb_Dredge<-dredge(Model_Pro_2019_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2019_nb_Subset<-subset(Model_Pro_2019_nb_Dredge,delta<2)
Model_Pro_2019_nb_Ave<-model.avg(Model_Pro_2019_nb_Subset)
sw(Model_Pro_2019_nb_Ave)
'''
                     cond(Pro_Harvest_2019) cond(Pro_Rivers) cond(Pro_Woodland) cond(Pro_Harvest_2019:Pro_Rivers)
Sum of weights:      1.00                   0.67             0.15               0.14                             
N containing models:    4                      3                1                  1 
'''
confint(Model_Pro_2019_nb_Ave)
'''
                                        2.5 %      97.5 %
cond((Int))                        3.22412980  3.33418969
cond(Pro_Harvest_2019)            -0.26441996 -0.15519212 *
cond(Pro_Rivers)                  -0.10323074  0.01126615
cond(Pro_Woodland)                -0.04648224  0.07297030
cond(Pro_Harvest_2019:Pro_Rivers) -0.03863956  0.05536188
'''
summary(Model_Pro_2019_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                   Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        3.279160   0.027921    0.028077 116.792   <2e-16 ***
cond(Pro_Harvest_2019)            -0.209806   0.027710    0.027865   7.529   <2e-16 ***
cond(Pro_Rivers)                  -0.030930   0.032142    0.032240   0.959    0.337    
cond(Pro_Woodland)                 0.001946   0.012527    0.012588   0.155    0.877    
cond(Pro_Harvest_2019:Pro_Rivers)  0.001187   0.009446    0.009494   0.125    0.901    
 
(conditional average) 
                                   Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        3.279160   0.027921    0.028077 116.792   <2e-16 ***
cond(Pro_Harvest_2019)            -0.209806   0.027710    0.027865   7.529   <2e-16 ***
cond(Pro_Rivers)                  -0.045982   0.029047    0.029209   1.574    0.115    
cond(Pro_Woodland)                 0.013244   0.030303    0.030473   0.435    0.664    
cond(Pro_Harvest_2019:Pro_Rivers)  0.008361   0.023847    0.023980   0.349    0.727 
'''
# Model assumption
plot(simulateResiduals(Model_Pro_2019_nb))
testDispersion(simulateResiduals(Model_Pro_2019_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Pro_2019_nb),x=df_Bram_Pro$X,y=df_Bram_Pro$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Pro_2019_nb),group=df_Bram_Pro$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Pro_2019_nb,Model_Pro_2019_nb_Dredge,Model_Pro_2019_nb_Subset,
   Model_Pro_2019_nb_Ave,res_scaled,Plantation,Loc)

#2018
Model_Pro_2018_nb<-glmmTMB(Bram_Abu_Sq~Pro_Rivers+Pro_Woodland+Pro_Harvest_2018+
                             Pro_Rivers:Pro_Harvest_2018+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2018_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.04         1.02      0.96
                Pro_Woodland 1.03         1.01      0.97
            Pro_Harvest_2018 1.01         1.00      0.99
 Pro_Rivers:Pro_Harvest_2018 1.03         1.01      0.97
'''
options(na.action="na.fail")
Model_Pro_2018_nb_Dredge<-dredge(Model_Pro_2018_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2018_nb_Subset<-subset(Model_Pro_2018_nb_Dredge,delta<2)
Model_Pro_2018_nb_Ave<-model.avg(Model_Pro_2018_nb_Subset)
sw(Model_Pro_2018_nb_Ave)
'''
                     cond(Pro_Harvest_2018) cond(Pro_Woodland) cond(Pro_Rivers) cond(Pro_Harvest_2018:Pro_Rivers)
Sum of weights:      1.00                   1.00               0.78             0.21                             
N containing models:    3                      3                  2                1 
'''
confint(Model_Pro_2018_nb_Ave)
'''
                                         2.5 %       97.5 %
cond((Int))                        3.206549038  3.350467476
cond(Pro_Harvest_2018)            -0.221678407 -0.102904563 *
cond(Pro_Rivers)                  -0.121201474 -0.001128266 *
cond(Pro_Woodland)                 0.008298699  0.129432361 *
cond(Pro_Harvest_2018:Pro_Rivers) -0.066678723  0.042768046
'''
summary(Model_Pro_2018_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                   Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        3.278508   0.036510    0.036715  89.297   <2e-16 ***
cond(Pro_Harvest_2018)            -0.162291   0.030131    0.030300   5.356    1e-07 ***
cond(Pro_Rivers)                  -0.047815   0.036928    0.037038   1.291   0.1967    
cond(Pro_Woodland)                 0.068866   0.030731    0.030902   2.229   0.0258 *  
cond(Pro_Harvest_2018:Pro_Rivers) -0.002562   0.013758    0.013825   0.185   0.8530    
 
(conditional average) 
                                  Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        3.27851    0.03651     0.03671  89.297   <2e-16 ***
cond(Pro_Harvest_2018)            -0.16229    0.03013     0.03030   5.356    1e-07 ***
cond(Pro_Rivers)                  -0.06116    0.03046     0.03063   1.997   0.0458 *  
cond(Pro_Woodland)                 0.06887    0.03073     0.03090   2.229   0.0258 *  
cond(Pro_Harvest_2018:Pro_Rivers) -0.01196    0.02776     0.02792   0.428   0.6685
'''
# Model assumption
plot(simulateResiduals(Model_Pro_2018_nb))
testDispersion(simulateResiduals(Model_Pro_2018_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Pro_2018_nb),x=df_Bram_Pro$X,y=df_Bram_Pro$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Pro_2018_nb),group=df_Bram_Pro$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Pro_2018_nb,Model_Pro_2018_nb_Dredge,Model_Pro_2018_nb_Subset,
   Model_Pro_2018_nb_Ave,res_scaled,Plantation,Loc)

#2017
Model_Pro_2017_nb<-glmmTMB(Bram_Abu_Sq~Pro_Rivers+Pro_Woodland+Pro_Harvest_2017+
                             Pro_Rivers:Pro_Harvest_2017+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2017_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.02         1.01      0.98
                Pro_Woodland 1.04         1.02      0.96
            Pro_Harvest_2017 1.06         1.03      0.95
 Pro_Rivers:Pro_Harvest_2017 1.03         1.02      0.97
'''
options(na.action="na.fail")
Model_Pro_2017_nb_Dredge<-dredge(Model_Pro_2017_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2017_nb_Subset<-subset(Model_Pro_2017_nb_Dredge,delta<2)
Model_Pro_2017_nb_Ave<-model.avg(Model_Pro_2017_nb_Subset)
sw(Model_Pro_2017_nb_Ave)
'''
                     cond(Pro_Harvest_2017) cond(Pro_Woodland) cond(Pro_Rivers)
Sum of weights:      1.00                   1.00               0.62            
N containing models:    2                      2                  1
'''
confint(Model_Pro_2017_nb_Ave)
'''
                              2.5 %       97.5 %
cond((Int))             3.236044187  3.355092475
cond(Pro_Harvest_2017) -0.155371237 -0.036000125 *
cond(Pro_Rivers)       -0.119634505  0.006095771
cond(Pro_Woodland)      0.008461412  0.138424586 *
'''
summary(Model_Pro_2017_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.29557    0.03020     0.03037 108.514  < 2e-16 ***
cond(Pro_Harvest_2017) -0.09569    0.03028     0.03045   3.142  0.00168 ** 
cond(Pro_Rivers)       -0.03524    0.03729     0.03738   0.943  0.34582    
cond(Pro_Woodland)      0.07344    0.03297     0.03315   2.215  0.02675 *  
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.29557    0.03020     0.03037 108.514  < 2e-16 ***
cond(Pro_Harvest_2017) -0.09569    0.03028     0.03045   3.142  0.00168 ** 
cond(Pro_Rivers)       -0.05677    0.03190     0.03207   1.770  0.07674 .  
cond(Pro_Woodland)      0.07344    0.03297     0.03315   2.215  0.02675 *
'''
# Model assumption
plot(simulateResiduals(Model_Pro_2017_nb))
testDispersion(simulateResiduals(Model_Pro_2017_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Pro_2017_nb),x=df_Bram_Pro$X,y=df_Bram_Pro$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Pro_2017_nb),group=df_Bram_Pro$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Pro_2017_nb,Model_Pro_2017_nb_Dredge,Model_Pro_2017_nb_Subset,
   Model_Pro_2017_nb_Ave,res_scaled,Plantation,Loc)

#2016
Model_Pro_2016_nb<-glmmTMB(Bram_Abu_Sq~Pro_Rivers+Pro_Woodland+Pro_Harvest_2016+
                             Pro_Rivers:Pro_Harvest_2016+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2016_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.02         1.01      0.98
                Pro_Woodland 1.03         1.01      0.97
            Pro_Harvest_2016 1.01         1.01      0.99
 Pro_Rivers:Pro_Harvest_2016 1.01         1.01      0.99
'''
options(na.action="na.fail")
Model_Pro_2016_nb_Dredge<-dredge(Model_Pro_2016_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2016_nb_Subset<-subset(Model_Pro_2016_nb_Dredge,delta<2)
Model_Pro_2016_nb_Ave<-model.avg(Model_Pro_2016_nb_Subset)
sw(Model_Pro_2016_nb_Ave)
'''
                     cond(Pro_Rivers) cond(Pro_Woodland) cond(Pro_Harvest_2016)
Sum of weights:      0.67             0.66               0.13                  
N containing models:    3                3                  1 
'''
confint(Model_Pro_2016_nb_Ave)
'''
                             2.5 %      97.5 %
cond((Int))             3.20138575 3.373620110
cond(Pro_Rivers)       -0.12255385 0.005704191
cond(Pro_Woodland)     -0.00899035 0.124701466
cond(Pro_Harvest_2016) -0.07362661 0.050729113
'''
summary(Model_Pro_2016_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                        Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.287503   0.043695    0.043938  74.821   <2e-16 ***
cond(Pro_Rivers)       -0.039348   0.038259    0.038362   1.026    0.305    
cond(Pro_Woodland)      0.038163   0.038864    0.038972   0.979    0.327    
cond(Pro_Harvest_2016) -0.001492   0.012024    0.012085   0.123    0.902    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             3.28750    0.04370     0.04394  74.821   <2e-16 ***
cond(Pro_Rivers)       -0.05842    0.03254     0.03272   1.786   0.0742 .  
cond(Pro_Woodland)      0.05786    0.03392     0.03411   1.696   0.0898 .  
cond(Pro_Harvest_2016) -0.01145    0.03155     0.03172   0.361   0.7182 
'''
# Model assumption
plot(simulateResiduals(Model_Pro_2016_nb))
testDispersion(simulateResiduals(Model_Pro_2016_nb))
testSpatialAutocorrelation(simulateResiduals(Model_Pro_2016_nb),x=df_Bram_Pro$X,y=df_Bram_Pro$Y)
# Recalculate residuals grouping per plantation
res_scaled<-recalculateResiduals(simulateResiduals(Model_Pro_2016_nb),group=df_Bram_Pro$PLANTATI_1)
Plantation<-st_read("Shapefiles/Plantations_Clipped.shp")
Loc<-st_centroid(Plantation)
Loc<-as.data.frame(st_coordinates(Loc))
testSpatialAutocorrelation(res_scaled,x=Loc$X,y=Loc$Y)
rm(Model_Pro_2016_nb,Model_Pro_2016_nb_Dredge,Model_Pro_2016_nb_Subset,
   Model_Pro_2016_nb_Ave,res_scaled,Plantation,Loc)

# Plot significant variables together
a<-ggplot(data=df_Bram_Pro,aes(x=Pro_Woodland,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to woodlands (m)",tag="a)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
b<-ggplot(data=df_Bram_Pro,aes(x=Pro_Rivers,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to rivers (m)",tag="b)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
c<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2020,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to harvested trees 2020 (m)",tag="c)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
d<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2019,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to harvested trees 2019 (m)",tag="d)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
e<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2018,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to harvested trees 2018 (m)",tag="e)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
f<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2017,y=Bram_Abu_Sq))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area (squared)\n",x="\nDistance to harvested trees 2017 (m)",tag="f)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
ggarrange(a,b,c,d,e,f,ncol=3,nrow=2) # 15 x 10
rm(a,b,c,d,e,f,df_Bram_Pro,df_Bram_Pro_scale)

###########################################################################
# What is the local impact of bramble on plants and grasshoppers?
###########################################################################

#### Bramble impact on plants ####
Veg<-read.csv("Excel_Sheets/Vegetation_Plot_Data.csv")
Veg$Itteration=NULL
Veg$Time=NULL
Veg$Date=NULL
Veg$Bram_Hei=NULL
Veg$Bram_Rich=NULL
clip<-Veg[,36:39]
Veg<-Veg%>%dplyr::select(!ends_with(c("Cov")))
Veg<-Veg%>%dplyr::select(!ends_with(c("Abu")))
Veg<-as.data.frame(cbind(Veg,clip))
rm(clip)

# Extract vegetation richness
Veg_Rich<-Veg%>%dplyr::select(ends_with(c("Rich")))
Veg_Rich<-as.data.frame(rowSums(Veg_Rich,na.rm=T))
colnames(Veg_Rich)<-c("Veg_Rich")
# Extract bramble cover and abundance
Bram<-Veg%>%dplyr::select(starts_with(c("Bram")))
# Select only growth form richness
Veg_Func_Rich<-Veg%>%dplyr::select(ends_with(c("Rich")))
# Create data frame
Veg<-as.data.frame(cbind(Veg[,1:3],"Veg_Rich"=Veg_Rich,Veg_Func_Rich,Bram))
rm(Veg_Rich,Bram,Veg_Func_Rich)

# Normality
hist(Veg$Veg_Rich)
qqnorm(Veg$Veg_Rich)
qqline(Veg$Veg_Rich,lty=2)
shapiro.test(Veg$Veg_Rich)

# Scaling explanatory variables
Clip<-Veg[,1:13]
features<-Veg[,14:15]
features<-as.data.frame(scale(features))
Veg_Scale<-as.data.frame(cbind(Clip,features))
rm(Clip,features)

# Modelling Bramble abundance
Model_Veg_Rich<-glmmTMB(Veg_Rich~Bram_Abu+(1|Site_ID)+(1|Plantation),
                        data=Veg_Scale,family="poisson")# nbinom2
confint(Model_Veg_Rich)
'''
                                          2.5 %      97.5 %    Estimate
cond.(Intercept)                     1.27495868  1.61768460  1.44632164
cond.Bram_Abu                       -0.11658762 -0.02713555 -0.07186159 *
Site_ID.cond.Std.Dev.(Intercept)     0.23304343  0.37873076  0.29708705
Plantation.cond.Std.Dev.(Intercept)  0.06446514  0.40723193  0.16202550
'''
summary(Model_Veg_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.44632    0.08743  16.542  < 2e-16 ***
Bram_Abu    -0.07186    0.02282  -3.149  0.00164 ** 
'''
r.squaredGLMM(Model_Veg_Rich)
'''
                 R2m       R2c
delta     0.01516688 0.3514905
lognormal 0.01617875 0.3749404
trigamma  0.01408491 0.3264159
'''
# Model testing
plot(residuals(Model_Veg_Rich))
plot(simulateResiduals(Model_Veg_Rich))
rm(Model_Veg_Rich)
# Plot variables
a<-ggplot(data=Veg,aes(x=Bram_Abu,y=Veg_Rich))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Plant species richness\n",x="\nBramble abundance")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
rm(Veg,Veg_Scale)

#### Bramble impact on grasshoppers richness ####
# Calculate averages from plot data
Veg<-read.csv("Excel_Sheets/Vegetation_Plot_Data.csv")
Veg$Bram_Hei=NULL

# Extract bramble abundance and cover
Bramble<-Veg[,40:41]
Bramble<-as.data.frame(cbind(Veg%>%dplyr::select(c("Site_ID")),Bramble))
Bramble<-Bramble %>% 
  group_by(Site_ID) %>% 
  summarise_each(funs(mean(.,na.rm=TRUE)))
Bramble$Site_ID=NULL
rm(Veg)

# Calculate grasshopper variables per conservation priority
hopper<-read.csv("Excel_Sheets/Grasshopper_Assemblage_Data.csv")
hop_Rich<-rowSums(hopper!=0)
hop_exShan<-exp(diversity(hopper,index="shannon"))
groups<-c("Low","Intermediate","Intermediate","Low","Intermediate","Intermediate","Low","Low","High","Low",
          "Low","Intermediate","Intermediate","Low","Low","Intermediate","Low","Low","Low","Low","Intermediate",
          "Low","Low","High","Intermediate","Low","Intermediate","Intermediate","High","Intermediate","High",
          "Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","Low","Intermediate",
          "Intermediate","High","Intermediate","Low","High","Intermediate","Low","High","High","High","Low",
          "Intermediate","High","Intermediate","Low","Intermediate","High","Intermediate","Intermediate",
          "Intermediate","Intermediate","Intermediate","Intermediate","High","High","Intermediate","High",
          "Low","Intermediate","Low","Intermediate","Low")
groups2<-c("Caelifera","Caelifera","Caelifera","Ensifera","Ensifera","Ensifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Ensifera","Ensifera","Caelifera",
           "Caelifera","Ensifera","Caelifera","Caelifera","Ensifera","Caelifera","Caelifera","Ensifera",
           "Ensifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Ensifera","Caelifera",
           "Caelifera","Ensifera","Ensifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera")
hopperT<-t(hopper)
hopperTT<-as.data.frame(cbind(hopperT,"Groups"=groups2))
hopperT<-as.data.frame(cbind(hopperT,"Groups"=groups))
Low<-subset(hopperT,Groups=="Low")
Low$Groups=NULL
Low<-as.data.frame(t(Low))
indx<-sapply(Low,is.factor)
Low[indx]<-lapply(Low[indx],function(x) as.numeric(as.character(x)))
Low_hop_Rich<-rowSums(Low!=0)
Low_exShan<-exp(diversity(Low,index="shannon"))
Intermediate<-subset(hopperT,Groups=="Intermediate")
Intermediate$Groups=NULL
Intermediate<-as.data.frame(t(Intermediate))
indx<-sapply(Intermediate,is.factor)
Intermediate[indx]<-lapply(Intermediate[indx],function(x) as.numeric(as.character(x)))
Intermediate_hop_Rich<-rowSums(Intermediate!=0)
Intermediate_exShan<-exp(diversity(Intermediate,index="shannon"))
High<-subset(hopperT,Groups=="High")
High$Groups=NULL
High<-as.data.frame(t(High))
indx<-sapply(High,is.factor)
High[indx]<-lapply(High[indx],function(x) as.numeric(as.character(x)))
High_hop_Rich<-rowSums(High!=0)
High_exShan<-exp(diversity(High,index="shannon"))
Caelifera<-subset(hopperTT,Groups=="Caelifera")
Caelifera$Groups=NULL
Caelifera<-as.data.frame(t(Caelifera))
indx<-sapply(Caelifera,is.factor)
Caelifera[indx]<-lapply(Caelifera[indx],function(x) as.numeric(as.character(x)))
Caelifera_hop_Rich<-rowSums(Caelifera!=0)
Caelifera_exShan<-exp(diversity(Caelifera,index="shannon"))
Ensifera<-subset(hopperTT,Groups=="Ensifera")
Ensifera$Groups=NULL
Ensifera<-as.data.frame(t(Ensifera))
indx<-sapply(Ensifera,is.factor)
Ensifera[indx]<-lapply(Ensifera[indx],function(x) as.numeric(as.character(x)))
Ensifera_hop_Rich<-rowSums(Ensifera!=0)
Ensifera_exShan<-exp(diversity(Ensifera,index="shannon"))

# Create data frame
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
df<-as.data.frame(cbind(GPS,hop_Rich,Low_hop_Rich,Intermediate_hop_Rich,High_hop_Rich,
                        Caelifera_hop_Rich,Ensifera_hop_Rich,hop_exShan,Low_exShan,Intermediate_exShan,
                        High_exShan,Caelifera_exShan,Ensifera_exShan,Bramble))
rm(hopper,High,Intermediate,Low,indx,hopperT,hopperTT,hop_Rich,groups,groups2,High_hop_Rich,Caelifera,
   Intermediate_hop_Rich,Low_hop_Rich,GPS,Caelifera_hop_Rich,Ensifera_hop_Rich,Ensifera,hop_exShan,
   Low_exShan,Intermediate_exShan,High_exShan,Caelifera_exShan,Ensifera_exShan,Bramble)

# Normality
hist(df$hop_Rich)
qqnorm(df$hop_Rich)
qqline(df$hop_Rich,lty=2)
shapiro.test(df$hop_Rich)
hist(df$hop_exShan)
qqnorm(df$hop_exShan)
qqline(df$hop_exShan,lty=2)
shapiro.test(df$hop_exShan)
hist(df$Low_hop_Rich)
qqnorm(df$Low_hop_Rich)
qqline(df$Low_hop_Rich,lty=2)
shapiro.test(df$Low_hop_Rich)
hist(df$Low_exShan)
qqnorm(df$Low_exShan)
qqline(df$Low_exShan,lty=2)
shapiro.test(df$Low_exShan)
hist(df$Intermediate_hop_Rich)
qqnorm(df$Intermediate_hop_Rich)
qqline(df$Intermediate_hop_Rich,lty=2)
shapiro.test(df$Intermediate_hop_Rich)
hist(df$Intermediate_exShan)
qqnorm(df$Intermediate_exShan)
qqline(df$Intermediate_exShan,lty=2)
shapiro.test(df$Intermediate_exShan)
hist(df$High_hop_Rich)
qqnorm(df$High_hop_Rich)
qqline(df$High_hop_Rich,lty=2)
shapiro.test(df$High_hop_Rich)
hist(df$High_exShan)
qqnorm(df$High_exShan)
qqline(df$High_exShan,lty=2)
shapiro.test(df$High_exShan)
hist(df$Caelifera_hop_Rich)
qqnorm(df$Caelifera_hop_Rich)
qqline(df$Caelifera_hop_Rich,lty=2)
shapiro.test(df$Caelifera_hop_Rich)
hist(df$Caelifera_exShan)
qqnorm(df$Caelifera_exShan)
qqline(df$Caelifera_exShan,lty=2)
shapiro.test(df$Caelifera_exShan)
hist(df$Ensifera_hop_Rich)
qqnorm(df$Ensifera_hop_Rich)
qqline(df$Ensifera_hop_Rich,lty=2)
shapiro.test(df$Ensifera_hop_Rich)
hist(df$Ensifera_exShan)
qqnorm(df$Ensifera_exShan)
qqline(df$Ensifera_exShan,lty=2)
shapiro.test(df$Ensifera_exShan)

# Scaling explanatory variables
Clip<-df[,1:16]
features<-df[,17:18]
features<-as.data.frame(scale(features))
df_Scale<-as.data.frame(cbind(Clip,features))
rm(Clip,features)

# Modelling Grasshopper richness
Model_hop_Rich<-glmmTMB(hop_Rich~Bram_Abu+(1|Plantation),
                        data=df_Scale,family="gaussian")
confint(Model_hop_Rich)
'''
                                         2.5 %     97.5 %   Estimate
cond.(Intercept)                    11.4310300 15.9573333 13.6941817
cond.Bram_Abu                       -1.5743472  0.6580187 -0.4581642
Plantation.cond.Std.Dev.(Intercept)  0.8346435  4.8654412  2.0151697
'''
summary(Model_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  13.6942     1.1547  11.860   <2e-16 ***
Bram_Abu     -0.4582     0.5695  -0.805    0.421 
'''
r.squaredGLMM(Model_hop_Rich)
'''
            R2m       R2c
[1,] 0.01087945 0.2213482
'''
# Model assumption
plot(residuals(Model_hop_Rich))
plot(simulateResiduals(Model_hop_Rich))
testDispersion(simulateResiduals(Model_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_hop_Rich)

# Modelling Low value grasshopper species
Model_Low_hop_Rich<-glmmTMB(Low_hop_Rich~Bram_Abu+(1|Plantation),
                            data=df_Scale,family="gaussian")
confint(Model_Low_hop_Rich)
'''
                                         2.5 %     97.5 %  Estimate
cond.(Intercept)                     6.5374984 10.0052618 8.2713801
cond.Bram_Abu                       -0.6295681  0.9187244 0.1445781
Plantation.cond.Std.Dev.(Intercept)  0.6844756  3.6788211 1.5868407
'''
summary(Model_Low_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   8.2714     0.8846   9.350   <2e-16 ***
Bram_Abu      0.1446     0.3950   0.366    0.714  
'''
r.squaredGLMM(Model_Low_hop_Rich)
'''
             R2m       R2c
[1,] 0.002147609 0.2608596
'''
# Model assumption
plot(residuals(Model_Low_hop_Rich))
plot(simulateResiduals(Model_Low_hop_Rich))
testDispersion(simulateResiduals(Model_Low_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Low_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Low_hop_Rich)

# Modelling Intermediate value grasshopper species
Model_Intermediate_hop_Rich<-glmmTMB(Intermediate_hop_Rich~Bram_Abu+(1|Plantation),
                                     data=df_Scale,family="gaussian")
confint(Model_Intermediate_hop_Rich)
'''
                                           2.5 %     97.5 %   Estimate
cond.(Intercept)                     4.139569398  5.1719119  4.6557407
cond.Bram_Abu                       -1.072983503 -0.1061482 -0.5895659 *
Plantation.cond.Std.Dev.(Intercept)  0.001754767 22.8061478  0.2000487
'''
summary(Model_Intermediate_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   4.6557     0.2634   17.68   <2e-16 ***
Bram_Abu     -0.5896     0.2466   -2.39   0.0168 *
'''
r.squaredGLMM(Model_Intermediate_hop_Rich)
'''
           R2m     R2c
[1,] 0.1061665 0.11839
'''
# Model assumption
plot(residuals(Model_Intermediate_hop_Rich))
plot(simulateResiduals(Model_Intermediate_hop_Rich))
testDispersion(simulateResiduals(Model_Intermediate_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Intermediate_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Intermediate_hop_Rich)
# Plot variables
b<-ggplot(data=df,aes(x=Bram_Abu,y=Intermediate_hop_Rich))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Intermediate value grasshopper species\n",x="\nBramble cover")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))

# Modelling High value grasshopper species
Model_High_hop_Rich<-glmmTMB(High_hop_Rich~Bram_Abu+(1|Plantation),
                             data=df_Scale,family="poisson")
confint(Model_High_hop_Rich)
'''
                                         2.5 %    97.5 %    Estimate
cond.(Intercept)                    -0.6283453 0.1476541 -0.24034562
cond.Bram_Abu                       -0.3852968 0.2814424 -0.05192722
Plantation.cond.Std.Dev.(Intercept)  0.0299690 1.6508780  0.22243014
'''
summary(Model_High_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)
(Intercept) -0.24035    0.19796  -1.214    0.225
Bram_Abu    -0.05193    0.17009  -0.305    0.760
'''
r.squaredGLMM(Model_High_hop_Rich)
'''
                  R2m        R2c
delta     0.002088372 0.04040655
lognormal 0.003141972 0.06079197
trigamma  0.001161981 0.02248242
'''
# Model assumption
plot(residuals(Model_High_hop_Rich))
plot(simulateResiduals(Model_High_hop_Rich))
testDispersion(simulateResiduals(Model_High_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_High_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_High_hop_Rich)

# Modelling Caelifera richness
Model_Caelifera_hop_Rich<-glmmTMB(Caelifera_hop_Rich~Bram_Abu+(1|Plantation),
                                  data=df_Scale,family="gaussian")
confint(Model_Caelifera_hop_Rich)
'''
                                         2.5 %    97.5 %    Estimate
cond.(Intercept)                     9.2823023 13.786692 11.53449736
cond.Bram_Abu                       -1.0635037  0.863758 -0.09987283
Plantation.cond.Std.Dev.(Intercept)  0.8991449  4.820534  2.08191237
'''
summary(Model_Caelifera_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept) 11.53450    1.14910  10.038   <2e-16 ***
Bram_Abu    -0.09987    0.49166  -0.203    0.839 
'''
r.squaredGLMM(Model_Caelifera_hop_Rich)
'''
              R2m       R2c
[1,] 0.0006468186 0.2817156
'''
# Model assumption
plot(residuals(Model_Caelifera_hop_Rich))
plot(simulateResiduals(Model_Caelifera_hop_Rich))
testDispersion(simulateResiduals(Model_Caelifera_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Caelifera_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Caelifera_hop_Rich)

# Modelling Ensifera richness 
Model_Ensifera_hop_Rich<-glmmTMB(Ensifera_hop_Rich~Bram_Abu+(1|Plantation),
                                 data=df_Scale,family="poisson")
confint(Model_Ensifera_hop_Rich)
'''
                                          2.5 %     97.5 %   Estimate
cond.(Intercept)                     0.43903436 1.08694776  0.7629911
cond.Bram_Abu                       -0.38344382 0.06356727 -0.1599383
Plantation.cond.Std.Dev.(Intercept)  0.04167159 1.23946586  0.2272675
'''
summary(Model_Ensifera_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   0.7630     0.1653   4.616 3.91e-06 ***
Bram_Abu     -0.1599     0.1140  -1.403    0.161 
'''
r.squaredGLMM(Model_Ensifera_hop_Rich)
'''
                 R2m       R2c
delta     0.04881181 0.1473705
lognormal 0.05727722 0.1729289
trigamma  0.04008006 0.1210079
'''
# Model assumption
plot(residuals(Model_Ensifera_hop_Rich))
plot(simulateResiduals(Model_Ensifera_hop_Rich))
testDispersion(simulateResiduals(Model_Ensifera_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Ensifera_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Ensifera_hop_Rich)

# Plot
ggarrange(a,b,ncol=2,nrow=1) # 8x4
rm(a,b,Model_Ensifera_hop_Rich,df_Scale,df)

###########################################################################
# How will uncontrolled bramble invasion impact the landscape?
###########################################################################

#### Grasshopper distribution modeling ####
# Load and clip variables
Study_ROI<-readOGR("Shapefiles/Plantations_Clipped.shp")
NDVI<-raster("Products/Cliped_SR_S2_VI.tif",band=11)
NDVI<-reclassify(NDVI,cbind(-Inf,-1,NA),right=FALSE)
NDVI<-reclassify(NDVI,cbind(1,Inf,NA))
Land_Use<-raster("Products/Map_ranger_SR_S2_VI.tif")
Rivers<-raster("Products/Distance_Rivers.tif")
Elv<-raster("Raster_Data/SUDEM/SUDEM_Clip_Pro.tif")
Aspect<-raster("Raster_Data/SUDEM/Aspect_Clip_Pro.tif")
Plantation_Edge<-raster("Products/Distance_Plantation.tif")
Fire_hist<-raster("Raster_Data/NBR/Min_NBR_2017_2020.tif")
NDVI<-crop(NDVI,Study_ROI)
NDVI<-mask(NDVI,Study_ROI)
Land_Use<-crop(Land_Use,Study_ROI)
Land_Use<-mask(Land_Use,Study_ROI)
Rivers<-crop(Rivers,Study_ROI)
Rivers<-mask(Rivers,Study_ROI)
Elv<-crop(Elv,Study_ROI)
Elv<-mask(Elv,Study_ROI)
Aspect<-crop(Aspect,Study_ROI)
Aspect<-mask(Aspect,Study_ROI)
Plantation_Edge<-crop(Plantation_Edge,Study_ROI)
Plantation_Edge<-mask(Plantation_Edge,Study_ROI)
Fire_hist<-crop(Fire_hist,Study_ROI)
Fire_hist<-mask(Fire_hist,Study_ROI)
rm(Study_ROI)

# Rescale variables
template_rst<-raster(xmn=200250.2,xmx=261647.7,ymn=6707441,ymx=6766053,resolution=10,
                     crs=projection("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
NDVI<-raster::resample(NDVI,template_rst,method="bilinear")
Land_Use<-raster::resample(Land_Use,template_rst,method="ngb")
Rivers<-raster::resample(Rivers,template_rst,method="bilinear")
Elv<-raster::resample(Elv,template_rst,method="bilinear")
Aspect<-raster::resample(Aspect,template_rst,method="bilinear")
Plantation_Edge<-raster::resample(Plantation_Edge,template_rst,method="bilinear")
Fire_hist<-raster::resample(Fire_hist,template_rst,method="bilinear")
rm(template_rst)

# Clip rasters
Plantation_Boundaries<-readOGR("Shapefiles/Plantation_Boundaries.shp")
NDVI<-mask(NDVI,Plantation_Boundaries,inverse=TRUE)
Land_Use<-mask(Land_Use,Plantation_Boundaries,inverse=TRUE)
Rivers<-mask(Rivers,Plantation_Boundaries,inverse=TRUE)
Elv<-mask(Elv,Plantation_Boundaries,inverse=TRUE)
Aspect<-mask(Aspect,Plantation_Boundaries,inverse=TRUE)
Plantation_Edge<-mask(Plantation_Edge,Plantation_Boundaries,inverse=TRUE)
Fire_hist<-mask(Fire_hist,Plantation_Boundaries,inverse=TRUE)
writeRaster(NDVI,filename="Products/SDM_NDVI.tif",format="GTiff",overwrite=TRUE)
writeRaster(Land_Use,filename="Products/SDM_Land_Use.tif",format="GTiff",overwrite=TRUE,datatype="INT1U")
writeRaster(Rivers,filename="Products/SDM_Rivers.tif",format="GTiff",overwrite=TRUE,datatype="INT2U")
writeRaster(Elv,filename="Products/SDM_Elv.tif",format="GTiff",overwrite=TRUE,datatype="INT2U")
writeRaster(Aspect,filename="Products/SDM_Aspect.tif",format="GTiff",overwrite=TRUE,datatype="INT2U")
writeRaster(Plantation_Edge,filename="Products/SDM_Plantation_Edge.tif",format="GTiff",overwrite=TRUE,datatype="INT2U")
writeRaster(Fire_hist,filename="Products/SDM_Fire_hist.tif",format="GTiff",overwrite=TRUE)
rm(NDVI,Land_Use,Aspect,Elv,Plantation_Edge,Plantation_Boundaries,Rivers,Fire_hist)

# Filter grasshopper data
hopper<-read.csv("Excel_Sheets/Grasshopper_Assemblage_Data.csv")
groups<-c("Low","Intermediate","Intermediate","Low","Intermediate","Intermediate","Low","Low","High","Low",
          "Low","Intermediate","Intermediate","Low","Low","Intermediate","Low","Low","Low","Low","Intermediate",
          "Low","Low","High","Intermediate","Low","Intermediate","Intermediate","High","Intermediate","High",
          "Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","Low","Intermediate",
          "Intermediate","High","Intermediate","Low","High","Intermediate","Low","High","High","High","Low",
          "Intermediate","High","Intermediate","Low","Intermediate","High","Intermediate","Intermediate",
          "Intermediate","Intermediate","Intermediate","Intermediate","High","High","Intermediate","High",
          "Low","Intermediate","Low","Intermediate","Low")
groups2<-c("Caelifera","Caelifera","Caelifera","Ensifera","Ensifera","Ensifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Ensifera","Ensifera","Caelifera",
           "Caelifera","Ensifera","Caelifera","Caelifera","Ensifera","Caelifera","Caelifera","Ensifera",
           "Ensifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Ensifera","Caelifera",
           "Caelifera","Ensifera","Ensifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera",
           "Caelifera","Caelifera","Caelifera","Caelifera","Caelifera","Caelifera")
hopperT<-t(hopper)
hopperTT<-as.data.frame(cbind(hopperT,"Groups"=groups2))
hopperT<-as.data.frame(cbind(hopperT,"Groups"=groups))
Low<-subset(hopperT,Groups=="Low")
Low$Groups=NULL
Low<-as.data.frame(t(Low))
indx<-sapply(Low,is.factor)
Low[indx]<-lapply(Low[indx],function(x) as.numeric(as.character(x)))
Intermediate<-subset(hopperT,Groups=="Intermediate")
Intermediate$Groups=NULL
Intermediate<-as.data.frame(t(Intermediate))
indx<-sapply(Intermediate,is.factor)
Intermediate[indx]<-lapply(Intermediate[indx],function(x) as.numeric(as.character(x)))
High<-subset(hopperT,Groups=="High")
High$Groups=NULL
High<-as.data.frame(t(High))
indx<-sapply(High,is.factor)
High[indx]<-lapply(High[indx],function(x) as.numeric(as.character(x)))
Caelifera<-subset(hopperTT,Groups=="Caelifera")
Caelifera$Groups=NULL
Caelifera<-as.data.frame(t(Caelifera))
indx<-sapply(Caelifera,is.factor)
Caelifera[indx]<-lapply(Caelifera[indx],function(x) as.numeric(as.character(x)))
Ensifera<-subset(hopperTT,Groups=="Ensifera")
Ensifera$Groups=NULL
Ensifera<-as.data.frame(t(Ensifera))
indx<-sapply(Ensifera,is.factor)
Ensifera[indx]<-lapply(Ensifera[indx],function(x) as.numeric(as.character(x)))
rm(indx,groups,groups2,hopperT,hopperTT)

# Merge GPS data with multivariate abundance data
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Plantation<-GPS$Plantation
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Proj_Cor<-as.data.frame(GPS@coords)

# Removing species with 3 or less occurrences
hopper<-hopper[,colSums(hopper!=0)>3]
Low<-Low[,colSums(Low!=0)>3]
Intermediate<-Intermediate[,colSums(Intermediate!=0)>3]
High<-High[,colSums(High!=0)>3]
Caelifera<-Caelifera[,colSums(Caelifera!=0)>3]
Ensifera<-Ensifera[,colSums(Ensifera!=0)>3]
hopper<-as.data.frame(cbind(Plantation,Proj_Cor,hopper))
Low<-as.data.frame(cbind(Plantation,Proj_Cor,Low))
Intermediate<-as.data.frame(cbind(Plantation,Proj_Cor,Intermediate))
High<-as.data.frame(cbind(Plantation,Proj_Cor,High))
Caelifera<-as.data.frame(cbind(Plantation,Proj_Cor,Caelifera))
Ensifera<-as.data.frame(cbind(Plantation,Proj_Cor,Ensifera))
rm(GPS,Proj_Cor,Plantation)

# Transform multivariate abundance data to occurrence matrix
hopper<-hopper %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Low<-Low %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Intermediate<-Intermediate %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
High<-High %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Caelifera<-Caelifera %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Ensifera<-Ensifera %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
write.csv(hopper,"Products/Occurrence_Data.csv",row.names=FALSE)
write.csv(Low,"Products/Occurrence_Low_Data.csv",row.names=FALSE)
write.csv(Intermediate,"Products/Occurrence_Intermediate_Data.csv",row.names=FALSE)
write.csv(High,"Products/Occurrence_High_Data.csv",row.names=FALSE)
write.csv(Caelifera,"Products/Occurrence_Caelifera_Data.csv",row.names=FALSE)
write.csv(Ensifera,"Products/Occurrence_Ensifera_Data.csv",row.names=FALSE)
rm(Caelifera,Ensifera,High,Intermediate,Low,hopper)

# Raster correlations
NDVI<-raster("Products/SDM_NDVI.tif")
Land_Use<-raster("Products/SDM_Land_Use.tif")
Rivers<-raster("Products/SDM_Rivers.tif")
Elv<-raster("Products/SDM_Elv.tif")
Aspect<-raster("Products/SDM_Aspect.tif")
Plantation_Edge<-raster("Products/SDM_Plantation_Edge.tif")
Fire_hist<-raster("Products/SDM_Fire_hist.tif")
Stack<-stack(NDVI,Land_Use,Rivers,Elv,Aspect,Plantation_Edge,Fire_hist)
jnk<-layerStats(Stack,'pearson',na.rm=T)
corr_matrix<-jnk[["pearson correlation coefficient"]]
corr_matrix
'''
                       SDM_NDVI     SDM_Land_Use   SDM_Rivers     SDM_Elv         SDM_Aspect    SDM_Plantation_Edge   SDM_Fire_hist
SDM_NDVI               0.99999999   0.54273970     -0.18829174    -0.22187075     0.15317677      -0.04467688         0.44902443
SDM_Land_Use          0.54273970    1.00000000     -0.16528368    -0.24076508     0.16390949      -0.05759286         0.52309394
SDM_Rivers            -0.18829174   -0.16528368    0.99999999     0.23764179      -0.03963255     0.16909927          -0.17321863
SDM_Elv               -0.22187075   -0.24076508    0.23764179     0.99999998      -0.05738102     0.28459798          -0.30358735
SDM_Aspect            0.15317677    0.16390949     -0.03963255    -0.05738102     0.99999999      0.04941656          0.15820725
SDM_Plantation_Edge   -0.04467688   -0.05759286    0.16909927     0.28459798      0.04941656      0.99999996          -0.27892126
SDM_Fire_hist         0.4490244     0.52309394     -0.1732186      -0.3035873      0.1582072      -0.2789213           1.0000000
'''
rm(corr_matrix,jnk,Stack,NDVI,Land_Use,Rivers,Elv,Aspect,Plantation_Edge,Fire_hist)

# Load environmental raster data
Env<-load_var(path="~/Desktop/Products/",
              files=c("SDM_NDVI.tif","SDM_Land_Use.tif","SDM_Rivers.tif","SDM_Aspect.tif",
                      "SDM_Plantation_Edge.tif","SDM_Elv.tif","SDM_Fire_hist.tif"),
              format=".tif",Norm=TRUE,categorical="SDM_Land_Use.tif")

# Load occurrence data for ALL SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="All_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM)

# Load occurrence data for LOW SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_Low_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="Low_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM)

# Load occurrence data for Intermediate SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_Intermediate_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="Intermediate_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM)

# Load occurrence data for High SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_High_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="High_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM)

# Load occurrence data for Caelifera SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_Caelifera_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="Caelifera_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM)

# Load occurrence data for Ensifera SPECIES
Occ<-load_occ(path="~/Desktop/Products/",Env=Env,
              file="Occurrence_Ensifera_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",rep=10,tmp=TRUE,
                     save=TRUE,name="Ensifera_Species_SDM",path="~/Desktop/Products/")
rm(Occ,SDM,Env)

#### Identify overlapping areas ####
# Load maps and select highest 95% of suitable pixel values
All_Species_Map<-raster("Products/All_Species_SDM/Stack/Rasters/Diversity.tif")
All_Species_Map<-All_Species_Map[[1]]>=quantile(All_Species_Map,.95)
All_Species_Map[All_Species_Map==0]<-NA
Low_Map<-raster("Products/Low_Species_SDM/Stack/Rasters/Diversity.tif")
Low_Map<-Low_Map[[1]]>=quantile(Low_Map,.95)
Low_Map[Low_Map==0]<-NA
Intermediate_Map<-raster("Products/Intermediate_Species_SDM/Stack/Rasters/Diversity.tif")
Intermediate_Map<-Intermediate_Map[[1]]>=quantile(Intermediate_Map,.95)
Intermediate_Map[Intermediate_Map==0]<-NA
High_Map<-raster("Products/High_Species_SDM/Stack/Rasters/Diversity.tif")
High_Map<-High_Map[[1]]>=quantile(High_Map,.95)
High_Map[High_Map==0]<-NA
Caelifera_Map<-raster("Products/Caelifera_Species_SDM/Stack/Rasters/Diversity.tif")
Caelifera_Map<-Caelifera_Map[[1]]>=quantile(Caelifera_Map,.95)
Caelifera_Map[Caelifera_Map==0]<-NA
Ensifera_Map<-raster("Products/Ensifera_Species_SDM/Stack/Rasters/Diversity.tif")
Ensifera_Map<-Ensifera_Map[[1]]>=quantile(Ensifera_Map,.95)
Ensifera_Map[Ensifera_Map==0]<-NA

# Select bramble from classification
Bramble_Map<-raster("Products/SDM_Land_Use.tif")
Bramble_Map[Bramble_Map>1]<-NA
# Vectorize bramble classification for map
Vector_Focal_Areas<-rasterToPolygons(Bramble_Map)
Vector_Focal_Areas$Obj_ID<-1:nrow(Vector_Focal_Areas)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas,dissolve=TRUE)
Bramble_Map_Class<-disaggregate(Vector_Focal_Areas_agg)
Bramble_Map_Class<-st_as_sf(Bramble_Map_Class)
st_write(Bramble_Map_Class,"Products/Bramble_Class_Vector.shp")
rm(Bramble_Map_Class,Vector_Focal_Areas,Vector_Focal_Areas_agg)

# Identify overlapping areas
Bram_Overall<-Bramble_Map+All_Species_Map
Bram_Low<-Bramble_Map+Low_Map
Bram_Intermediate<-Bramble_Map+Intermediate_Map
Bram_High<-Bramble_Map+High_Map
Bram_Caelifera<-Bramble_Map+Caelifera_Map
Bram_Ensifera<-Bramble_Map+Ensifera_Map

# Save
writeRaster(Bram_Overall,filename="Products/Overlap/Bram_Overall.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Low,filename="Products/Overlap/Bram_Low.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Intermediate,filename="Products/Overlap/Bram_Intermediate.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_High,filename="Products/Overlap/Bram_High.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Caelifera,filename="Products/Overlap/Bram_Caelifera.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Ensifera,filename="Products/Overlap/Bram_Ensifera.tif",format="GTiff",overwrite=TRUE)

# Calculate percentage overlap
All_Species_Map_area<-cellStats(All_Species_Map,stat="sum")
Bram_Overall_area<-cellStats(Bram_Overall,"sum")
Bram_Low_area<-cellStats(Bram_Low,"sum")
Bram_Intermediate_area<-cellStats(Bram_Intermediate,"sum")
Bram_High_area<-cellStats(Bram_High,"sum")
Bram_Caelifera_area<-cellStats(Bram_Caelifera,"sum")
Bram_Ensifera_area<-cellStats(Bram_Ensifera,"sum")
Overall<-((Bram_Overall_area/All_Species_Map_area)*100)
Low<-((Bram_Low_area/All_Species_Map_area)*100)
Intermediate<-((Bram_Intermediate_area/All_Species_Map_area)*100)
High<-((Bram_High_area/All_Species_Map_area)*100)
Caelifera<-((Bram_Caelifera_area/All_Species_Map_area)*100)
Ensifera<-((Bram_Ensifera_area/All_Species_Map_area)*100)

# Export table
Percentages<-as.data.frame(rbind(Overall,Low,Intermediate,High,Caelifera,Ensifera))
Percentages<-add_rownames(Percentages,var="Groups")
colnames(Percentages)<-c("Groups","Percentage")
write.csv(Percentages,"Products/Percentages.csv",row.names=FALSE)

# Calculate future invasion habitat overlap
Bramble_Map<-rasterToPolygons(Bramble_Map)
Bramble_Map_agg<-aggregate(Bramble_Map,dissolve=TRUE)
Bramble_Map_dis<-disaggregate(Bramble_Map_agg)
Bramble_Map_buf<-buffer(Bramble_Map_dis,width=0.001,dissolve=TRUE)
Bramble_Map<-disaggregate(Bramble_Map_buf)
Bramble_Map<-as(Bramble_Map,"SpatialPolygonsDataFrame")
rm(Bramble_Map_agg,Bramble_Map_dis,Bramble_Map_buf)
Bramble_Buf<-buffer(Bramble_Map,width=10)
template_rst<-raster(xmn=200250.2,xmx=261647.7,ymn=6707441,ymx=6766053,resolution=10,
                     crs=projection("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
Bramble_Buf<-rasterize(Bramble_Buf,template_rst)

# Identify overlapping areas
Bram_Overall<-Bramble_Buf+All_Species_Map
Bram_Low<-Bramble_Buf+Low_Map
Bram_Intermediate<-Bramble_Buf+Intermediate_Map
Bram_High<-Bramble_Buf+High_Map
Bram_Caelifera<-Bramble_Buf+Caelifera_Map
Bram_Ensifera<-Bramble_Buf+Ensifera_Map

# Save
writeRaster(Bram_Overall,filename="Products/Overlap/Future_Bram_Overall.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Low,filename="Products/Overlap/Future_Bram_Low.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Intermediate,filename="Products/Overlap/Future_Bram_Intermediate.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_High,filename="Products/Overlap/Future_Bram_High.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Caelifera,filename="Products/Overlap/Future_Bram_Caelifera.tif",format="GTiff",overwrite=TRUE)
writeRaster(Bram_Ensifera,filename="Products/Overlap/Future_Bram_Ensifera.tif",format="GTiff",overwrite=TRUE)

# Calculate percentage overlap
All_Species_Map_area<-cellStats(All_Species_Map,stat="sum")
Bram_Overall_area<-cellStats(Bram_Overall,"sum")
Bram_Low_area<-cellStats(Bram_Low,"sum")
Bram_Intermediate_area<-cellStats(Bram_Intermediate,"sum")
Bram_High_area<-cellStats(Bram_High,"sum")
Bram_Caelifera_area<-cellStats(Bram_Caelifera,"sum")
Bram_Ensifera_area<-cellStats(Bram_Ensifera,"sum")
Overall<-((Bram_Overall_area/All_Species_Map_area)*100)
Low<-((Bram_Low_area/All_Species_Map_area)*100)
Intermediate<-((Bram_Intermediate_area/All_Species_Map_area)*100)
High<-((Bram_High_area/All_Species_Map_area)*100)
Caelifera<-((Bram_Caelifera_area/All_Species_Map_area)*100)
Ensifera<-((Bram_Ensifera_area/All_Species_Map_area)*100)

# Export table
Percentages<-as.data.frame(rbind(Overall,Low,Intermediate,High,Caelifera,Ensifera))
Percentages<-add_rownames(Percentages,var="Groups")
colnames(Percentages)<-c("Groups","Percentage")
write.csv(Percentages,"Products/Percentages_Future.csv",row.names=FALSE)
rm(All_Species_Map_area,Bram_Overall_area,Bram_Low_area,Bram_Intermediate_area,
   Bram_High_area,Bram_Caelifera_area,Bram_Ensifera_area,Overall,Low,Intermediate,
   High,Caelifera,Ensifera,Bramble_Map,All_Species_Map,Low_Map,Intermediate_Map,
   High_Map,Caelifera_Map,Ensifera_Map,Bram_Overall,Bram_Low,Bram_Intermediate,
   Bram_High,Bram_Caelifera,Bram_Ensifera,Percentages,template_rst,Bramble_Buf)

# Filter potential vulnerable areas
All_Species_Map<-raster("Products/Overlap/Future_Bram_Overall.tif")
Vector_Focal_Areas<-rasterToPolygons(All_Species_Map)
Vector_Focal_Areas$Obj_ID<-1:nrow(Vector_Focal_Areas)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas,dissolve=TRUE)
All_Species_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas,Vector_Focal_Areas_agg)
Low_Map<-raster("Products/Overlap/Future_Bram_Low.tif")
Vector_Focal_Areas_Low<-rasterToPolygons(Low_Map)
Vector_Focal_Areas_Low$Obj_ID<-1:nrow(Vector_Focal_Areas_Low)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Low,dissolve=TRUE)
Low_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas_Low,Vector_Focal_Areas_agg)
Intermediate_Map<-raster("Products/Overlap/Future_Bram_Intermediate.tif")
Vector_Focal_Areas_Intermediate<-rasterToPolygons(Intermediate_Map)
Vector_Focal_Areas_Intermediate$Obj_ID<-1:nrow(Vector_Focal_Areas_Intermediate)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Intermediate,dissolve=TRUE)
Intermediate_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas_Intermediate,Vector_Focal_Areas_agg)
High_Map<-raster("Products/Overlap/Future_Bram_High.tif")
Vector_Focal_Areas_High<-rasterToPolygons(High_Map)
Vector_Focal_Areas_High$Obj_ID<-1:nrow(Vector_Focal_Areas_High)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_High,dissolve=TRUE)
High_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas_High,Vector_Focal_Areas_agg)
Caelifera_Map<-raster("Products/Overlap/Future_Bram_Caelifera.tif")
Vector_Focal_Areas_Caelifera<-rasterToPolygons(Caelifera_Map)
Vector_Focal_Areas_Caelifera$Obj_ID<-1:nrow(Vector_Focal_Areas_Caelifera)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Caelifera,dissolve=TRUE)
Caelifera_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas_Caelifera,Vector_Focal_Areas_agg)
Ensifera_Map<-raster("Products/Overlap/Future_Bram_Ensifera.tif")
Vector_Focal_Areas_Ensifera<-rasterToPolygons(Ensifera_Map)
Vector_Focal_Areas_Ensifera$Obj_ID<-1:nrow(Vector_Focal_Areas_Ensifera)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Ensifera,dissolve=TRUE)
Ensifera_Map<-disaggregate(Vector_Focal_Areas_agg)
rm(Vector_Focal_Areas_Ensifera,Vector_Focal_Areas_agg)

# Select only large focal areas
Vector_Focal_Areas_features<-as(All_Species_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
All_Species_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)
Vector_Focal_Areas_features<-as(Low_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Low_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)
Vector_Focal_Areas_features<-as(Intermediate_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Intermediate_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)
Vector_Focal_Areas_features<-as(High_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
High_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)
Vector_Focal_Areas_features<-as(Caelifera_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Caelifera_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)
Vector_Focal_Areas_features<-as(Ensifera_Map,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Ensifera_Map<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=600)

# Save
st_write(All_Species_Map,"Products/Overlap/Filter_Future_Bram_Overall.shp")
st_write(Low_Map,"Products/Overlap/Filter_Future_Bram_Low.shp")
st_write(Intermediate_Map,"Products/Overlap/Filter_Future_Bram_Intermediate.shp")
st_write(High_Map,"Products/Overlap/Filter_Future_Bram_High.shp")
st_write(Caelifera_Map,"Products/Overlap/Filter_Future_Bram_Caelifera.shp")
st_write(Ensifera_Map,"Products/Overlap/Filter_Future_Bram_Ensifera.shp")
rm(All_Species_Map,Low_Map,Intermediate_Map,High_Map,Caelifera_Map,Ensifera_Map,
   Vector_Focal_Areas_features)
