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

# Community
library(mvabund)
library(ade4)
library(vegan)

# Set working directory
setwd("~/")

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
predict_ms210<-stack("Raster_Data/Super_Resolve/Bramble_decjan2019_Results/SuperResolution/predict_ms210.tif")
predict_ms210<-crop(predict_ms210,Study_ROI)
predict_ms210<-mask(predict_ms210,Study_ROI)
predict_ms220<-stack("Raster_Data/Super_Resolve/Bramble_decjan2019_Results/SuperResolution/predict_ms220.tif")
predict_ms220<-crop(predict_ms220,Study_ROI)
predict_ms220<-mask(predict_ms220,Study_ROI)
SR_S2<-stack(predict_ms210,predict_ms220)
SR_S2<-subset(SR_S2,order(c(3,2,1,4,5,6,7,8,9,10)))
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
writeRaster(SR_S2,filename="Products/Cliped_SR_S2.tif",format="GTiff",overwrite=TRUE)
rm(Study_ROI,predict_ms210,predict_ms220,SR_S2)

#### Process training data ####
# Load points
Training_Points<-read.csv("Excel_Sheets/Training_Points.csv")
# Project
Training_Points<-st_as_sf(Training_Points,coords=c("x","y"),crs=4326)
Training_Points<-st_transform(Training_Points,crs=32736)

# Extract training data for Sentinel imagery
S2<-stack("Products/Cliped_S2.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
df_S2<-raster::extract(S2,Training_Points,method="simple",df=TRUE)
df_S2$ID=NULL
Class<-as.matrix(Training_Points$Class)
# Save csv file
df_S2<-as.data.frame(cbind(df_S2,"Class"=Class))
df_S2<-df_S2[complete.cases(df_S2),]
write.csv(df_S2,"Products/df_S2.csv",row.names=FALSE)
rm(Class,S2,Training_Points,df_S2)

# Create buffers for extracting SR_S2 and PS
# Load points
Training_Points<-read.csv("Excel_Sheets/Training_Points.csv")
# Project
Training_Points<-st_as_sf(Training_Points,coords=c("x","y"),crs=4326)
Training_Points<-st_transform(Training_Points,crs=32736)
# Create buffer
Training_Points<-st_buffer(Training_Points,2.5)
st_write(Training_Points,"Products/Shapefiles/Training_Points_Buffer.shp")

# Extract training data for PlanetScope imagery
PS<-stack("Products/Cliped_PS.tif")
names(PS)<-c("Blue","Green","Red","Infrared")
template_rst<-raster(extent(PS[[1]]),resolution=3,crs=projection(PS[[1]]))
Training_Points_rst<-fasterize(Training_Points,template_rst,field="ID")
Training_Points_df<-as.data.frame(rasterToPoints(Training_Points_rst))
colnames(Training_Points_df)<-c("x","y","ID")
points<-SpatialPointsDataFrame(coords=Training_Points_df[,1:2],data=Training_Points_df,proj4string=Training_Points_rst@crs)
writeOGR(points,"Products/Shapefiles/","Raster_to_Points_PS",driver="ESRI Shapefile")
df_PS<-raster::extract(PS,points,method="simple",df=TRUE)
df_PS$ID=NULL
# Assign Class labels
ID<-as.matrix(points@data[["ID"]])
Class<-replace(ID,ID==1,"Bramble")
Class<-replace(Class,Class==2,"Water")
Class<-replace(Class,Class==3,"Plantation")
Class<-replace(Class,Class==4,"Ground")
Class<-replace(Class,Class==5,"Shrubland")
Class<-replace(Class,Class==6,"Grassland")
Class<-replace(Class,Class==7,"Woodland")
Class<-replace(Class,Class==8,"Thicket")
# Save csv file
df_PS<-as.data.frame(cbind(df_PS,"Class"=Class))
df_PS<-df_PS[complete.cases(df_PS),]
write.csv(df_PS,"Products/df_PS.csv",row.names=FALSE)
rm(ID,PS,template_rst,Training_Points_rst,Training_Points_df,points,df_PS,Class)

# Extract training data for super resolution imagery
SR_S2<-stack("Products/Cliped_SR_S2.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
template_rst<-raster(extent(SR_S2[[1]]),resolution=2.5,crs=projection(SR_S2[[1]]))
Training_Points_rst<-fasterize(Training_Points,template_rst,field="ID")
Training_Points_df<-as.data.frame(rasterToPoints(Training_Points_rst))
colnames(Training_Points_df)<-c("x","y","ID")
points<-SpatialPointsDataFrame(coords=Training_Points_df[,1:2],data=Training_Points_df,proj4string=Training_Points_rst@crs)
writeOGR(points,"Products/Shapefiles/","Raster_to_Points_SR_S2",driver="ESRI Shapefile")
df_SR_S2<-raster::extract(SR_S2,points,method="simple",df=TRUE)
df_SR_S2$ID=NULL
# Assign Class labels
ID<-as.matrix(points@data[["ID"]])
Class<-replace(ID,ID==1,"Bramble")
Class<-replace(Class,Class==2,"Water")
Class<-replace(Class,Class==3,"Plantation")
Class<-replace(Class,Class==4,"Ground")
Class<-replace(Class,Class==5,"Shrubland")
Class<-replace(Class,Class==6,"Grassland")
Class<-replace(Class,Class==7,"Woodland")
Class<-replace(Class,Class==8,"Thicket")
# Save csv file
df_SR_S2<-as.data.frame(cbind(df_SR_S2,"Class"=Class))
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]
write.csv(df_SR_S2,"Products/df_SR_S2.csv",row.names=FALSE)
rm(ID,SR_S2,template_rst,Training_Points_rst,Training_Points_df,Training_Points,df_SR_S2,points,Class)

#### Classification ####
# Load and classify Sentinel 2
df_S2<-read.csv("Products/df_S2.csv")
df_S2<-df_S2[complete.cases(df_S2),]
S2<-stack("Products/Cliped_S2.tif")
names(S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
             "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
# Inspect class balance
plot(df_S2$Class,main="Class Frequencies S2")
# Split dataset
set.seed(321)
df_Split<-createDataPartition(df_S2$Class,p=0.7,list=FALSE)
df_S2_train<-df_S2[df_Split,]
df_S2_test<-df_S2[-df_Split,]
rm(df_Split,df_S2)

# Classification: Random Forest (ranger)
set.seed(333)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2,
                           method="ranger",importance="impurity",
                           data=df_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_S2.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_S2_test)
confusionMatrix(data=predict_ranger,df_S2_test$Class)#,mode="prec_recall"
# Variable importance
plot(varImp(model_ranger),main="Variable Importance Ranger S2")
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_S2.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map)

# Load and classify PlanetScope
df_PS<-read.csv("Products/df_PS.csv")
df_PS<-df_PS[complete.cases(df_PS),]
PS<-stack("Products/Cliped_PS.tif")
names(PS)<-c("Blue","Green","Red","Infrared")
# Inspect class balance
plot(df_PS$Class,main="Class Frequencies PS")
# Split dataset
set.seed(321)
df_Split<-createDataPartition(df_PS$Class,p=0.7,list=FALSE)
df_PS_train<-df_PS[df_Split,]
df_PS_test<-df_PS[-df_Split,]
rm(df_Split,df_PS)

# Classification: Random Forest (ranger)
set.seed(333)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared,
                           method="ranger",importance="impurity",
                           data=df_PS_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_PS.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_PS_test)
confusionMatrix(data=predict_ranger,df_PS_test$Class)#,mode="prec_recall"
# Variable importance
plot(varImp(model_ranger),main="Variable Importance Ranger PS")
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(PS,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_PS.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map)

# Load and classify Super-resolve
df_SR_S2<-read.csv("Products/df_SR_S2_SR_S2.csv")
df_SR_S2<-df_SR_S2[complete.cases(df_SR_S2),]
SR_S2<-stack("Products/Cliped_SR_S2.tif")
names(SR_S2)<-c("Blue","Green","Red","Infrared","RedEdge1",
                "RedEdge2","RedEdge3","RedEdge4","SWIR1","SWIR2")
# Inspect class balance
plot(df_SR_S2$Class,main="Class Frequencies SR S2")
# Split dataset
set.seed(321)
df_SR_S2_Split<-createDataPartition(df_SR_S2$Class,p=0.7,list=FALSE)
df_SR_S2_train<-df_SR_S2[df_SR_S2_Split,]
df_SR_S2_test<-df_SR_S2[-df_SR_S2_Split,]
rm(df_SR_S2_Split,df_SR_S2)

# Classification: Random Forest (ranger)
set.seed(333)
model_ranger<-caret::train(Class~Blue+Green+Red+Infrared+RedEdge1+RedEdge2+RedEdge3+RedEdge4+SWIR1+SWIR2,
                           method="ranger",importance="impurity",
                           data=df_SR_S2_train,tuneLength=10,metric="Kappa",
                           trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE))
# Saving model
saveRDS(model_ranger,file="Products/model_ranger_SR_S2.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,df_SR_S2_test)
confusionMatrix(data=predict_ranger,df_SR_S2_test$Class)#,mode="prec_recall"
# Variable importance
plot(varImp(model_ranger),main="Variable Importance Ranger SR S2")
# Predict and save
beginCluster(7,type='SOCK')
ranger_map<-clusterR(SR_S2,predict,args=list(model=model_ranger,type="raw"))
endCluster()
writeRaster(ranger_map,filename="Products/Map_ranger_SR_S2.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,model_ranger,ranger_map,df_SR_S2_test,df_SR_S2_train,SR_S2)

#### Create bramble response variable ####
# Extract bramble class
Bramble_Class<-raster("Products/Map_ranger_SR_S2.tif")
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
  filter(Area_sqm>=31.25) #2.5 pixel size = 6.25m squared: 
st_write(Vector_Bramble_Class_Filter,"Products/Shapefiles/Vector_Bramble_Class_Dis_Filter.shp")

# Generate random points
#https://www.jla-data.net/eng/creating-and-pruning-random-points-and-polygons/
Vector_Bramble_Class<-st_as_sf(Vector_Bramble_Class)
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
rm(Random_Points_Bramble_Class,Vector_Bramble_Class,i,buffer_size,buffer,offending)
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
Woodland<-raster("Products/Map_ranger_SR_S2.tif")
Woodland<-raster::resample(Woodland,template_rst,method="ngb")
Woodland[Woodland<7]<-NA
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
NBR_2020<-raster("Raster_Data/NBR/Min_NBR_Landsat_2020.tif")
NBR_2019<-raster("Raster_Data/NBR/Min_NBR_Landsat_2019.tif")
NBR_2018<-raster("Raster_Data/NBR/Min_NBR_Landsat_2018.tif")
NBR_2017<-raster("Raster_Data/NBR/Min_NBR_Landsat_2017.tif")
NBR_2016<-raster("Raster_Data/NBR/Min_NBR_Landsat_2016.tif")
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

# Outliers Conservation dataframe
df_Bram_Con<-df_Bram_Con[-c(31,55,256,266,270),]
dotchart(df_Bram_Con$Bram_Abu,ylab="Bramble Abundance")
dotchart(df_Bram_Con$Con_Roads,ylab="Distance to roads")
dotchart(df_Bram_Con$Con_Rivers,ylab="Distance to rivers")
dotchart(df_Bram_Con$Con_Woodland,ylab="Distance to Woodlands")
dotchart(df_Bram_Con$Con_Plantation,ylab="Distance to Plantation boundary")
dotchart(df_Bram_Con$Con_Harvest_2020,ylab="Distance to 2020 harvesting")
dotchart(df_Bram_Con$Con_Harvest_2019,ylab="Distance to 2019 harvesting")
dotchart(df_Bram_Con$Con_Harvest_2018,ylab="Distance to 2018 harvesting")
dotchart(df_Bram_Con$Con_Harvest_2017,ylab="Distance to 2017 harvesting")
dotchart(df_Bram_Con$Con_Harvest_2016,ylab="Distance to 2016 harvesting")
dotchart(df_Bram_Con$Con_Fire_2020,ylab="Fire 2020")
dotchart(df_Bram_Con$Con_Fire_2019,ylab="Fire 2019")
dotchart(df_Bram_Con$Con_Fire_2018,ylab="Fire 2018")
dotchart(df_Bram_Con$Con_Fire_2017,ylab="Fire 2017")
dotchart(df_Bram_Con$Con_Fire_2016,ylab="Fire 2016")
df_Bram_Pro<-df_Bram_Pro[-c(95,134,217),]
dotchart(df_Bram_Pro$Bram_Abu,ylab="Bramble Abundance")
dotchart(df_Bram_Pro$Pro_Roads,ylab="Distance to roads")
dotchart(df_Bram_Pro$Pro_Rivers,ylab="Distance to rivers")
dotchart(df_Bram_Pro$Pro_Woodland,ylab="Distance to Woodlands")
dotchart(df_Bram_Pro$Pro_Plantation,ylab="Distance to Plantation boundary")
dotchart(df_Bram_Pro$Pro_Harvest_2020,ylab="Distance to 2020 harvesting")
dotchart(df_Bram_Pro$Pro_Harvest_2019,ylab="Distance to 2019 harvesting")
dotchart(df_Bram_Pro$Pro_Harvest_2018,ylab="Distance to 2018 harvesting")
dotchart(df_Bram_Pro$Pro_Harvest_2017,ylab="Distance to 2017 harvesting")
dotchart(df_Bram_Pro$Pro_Harvest_2016,ylab="Distance to 2016 harvesting")

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
Con_Points_GPS<-Con_Points_GPS[-c(31,55,256,266,270),]
Pro_Points_GPS<-Pro_Points_GPS[-c(95,134,217),]
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
Plantation_Con<-Plantation_Con[-c(31,55,256,266,270),]
Plantation_Pro<-Plantation_Pro[-c(95,134,217),]
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
Con_Points_GPS<-Con_Points_GPS[-c(31,55,256,266,270),]
Pro_Points_GPS<-Pro_Points_GPS[-c(95,134,217),]
df_Bram_Con<-as.data.frame(cbind(df_Bram_Con,Con_Points_GPS))
df_Bram_Pro<-as.data.frame(cbind(df_Bram_Pro,Pro_Points_GPS))
df_Bram_Con_scale<-as.data.frame(cbind(df_Bram_Con_scale,Con_Points_GPS))
df_Bram_Pro_scale<-as.data.frame(cbind(df_Bram_Pro_scale,Pro_Points_GPS))
df_Bram_Con<-df_Bram_Con[complete.cases(df_Bram_Con),]
df_Bram_Con_scale<-df_Bram_Con_scale[complete.cases(df_Bram_Con_scale),]
df_Bram_Pro<-df_Bram_Pro[complete.cases(df_Bram_Pro),]
df_Bram_Pro_scale<-df_Bram_Pro_scale[complete.cases(df_Bram_Pro_scale),]
rm(Con_Points_GPS,Pro_Points_GPS)

# Model building: Conservation areas
# 2020
Model_Con_2020_nb<-glmmTMB(Bram_Abu~Con_Rivers+Con_Woodland+Con_Harvest_2020+Con_Fire_2020+
                             Con_Rivers:Con_Fire_2020+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2020_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.20         1.09      0.84
             Con_Woodland 1.32         1.15      0.76
         Con_Harvest_2020 1.17         1.08      0.86
            Con_Fire_2020 1.37         1.17      0.73
 Con_Rivers:Con_Fire_2020 1.16         1.08      0.86
'''
options(na.action="na.fail")
Model_Con_2020_nb_Dredge<-dredge(Model_Con_2020_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")Model_Con_2020_nb_Ave
Model_Con_2020_nb_Subset<-subset(Model_Con_2020_nb_Dredge,delta<2)
Model_Con_2020_nb_Ave<-model.avg(Model_Con_2020_nb_Subset)
importance(Model_Con_2020_nb_Ave)
'''
                     cond(Con_Fire_2020) cond(Con_Rivers) cond(Con_Fire_2020:Con_Rivers) cond(Con_Harvest_2020)
Sum of weights:      1.00                0.83             0.24                           0.17                  
N containing models:    4                   3                1                              1
'''
confint(Model_Con_2020_nb_Ave)
'''
                                     2.5 %       97.5 %
cond((Int))                     6.91040010  7.176241899
cond(Con_Fire_2020)            -0.29163077 -0.096622927 *
cond(Con_Rivers)               -0.19911693 -0.005521042 *
cond(Con_Fire_2020:Con_Rivers) -0.14297391  0.046366343
cond(Con_Harvest_2020)         -0.06116333  0.112394848
'''
summary(Model_Con_2020_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.038880   0.066318    0.066626 105.647  < 2e-16 ***
cond(Con_Fire_2020)            -0.199850   0.048914    0.049142   4.067 4.77e-05 ***
cond(Con_Rivers)               -0.108396   0.049156    0.049377   2.195   0.0281 *  
cond(Con_Fire_2020:Con_Rivers) -0.015702   0.035502    0.035594   0.441   0.6591    
cond(Con_Harvest_2020)          0.006016   0.023399    0.023481   0.256   0.7978    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.03888    0.06632     0.06663 105.647  < 2e-16 ***
cond(Con_Fire_2020)            -0.19985    0.04891     0.04914   4.067 4.77e-05 ***
cond(Con_Rivers)               -0.10840    0.04916     0.04938   2.195   0.0281 *  
cond(Con_Fire_2020:Con_Rivers) -0.05171    0.04784     0.04807   1.076   0.2820    
cond(Con_Harvest_2020)          0.02855    0.04422     0.04442   0.643   0.5204
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
Model_Con_2019_nb<-glmmTMB(Bram_Abu~Con_Rivers+Con_Woodland+Con_Harvest_2019+Con_Fire_2019+
                             Con_Rivers:Con_Fire_2019+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2019_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.19         1.09      0.84
             Con_Woodland 1.41         1.19      0.71
         Con_Harvest_2019 1.20         1.09      0.84
            Con_Fire_2019 1.35         1.16      0.74
 Con_Rivers:Con_Fire_2019 1.10         1.05      0.91
'''
options(na.action="na.fail")
Model_Con_2019_nb_Dredge<-dredge(Model_Con_2019_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2019_nb_Subset<-subset(Model_Con_2019_nb_Dredge,delta<2)
Model_Con_2019_nb_Ave<-model.avg(Model_Con_2019_nb_Subset)
importance(Model_Con_2019_nb_Ave)
'''
                     cond(Con_Fire_2019) cond(Con_Rivers) cond(Con_Fire_2019:Con_Rivers)
Sum of weights:      1.00                0.67             0.19                          
N containing models:    3                   2                1   
'''
confint(Model_Con_2019_nb_Ave)
'''
                                    2.5 %      97.5 %
cond((Int))                     6.9242216  7.17674058
cond(Con_Fire_2019)            -0.2409542 -0.04639960 *
cond(Con_Rivers)               -0.1735418  0.01376803
cond(Con_Fire_2019:Con_Rivers) -0.0627881  0.10472244
'''
summary(Model_Con_2019_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.050481   0.064120    0.064419 109.447  < 2e-16 ***
cond(Con_Fire_2019)            -0.143677   0.049407    0.049632   2.895  0.00379 ** 
cond(Con_Rivers)               -0.053455   0.054098    0.054229   0.986  0.32427    
cond(Con_Fire_2019:Con_Rivers)  0.003983   0.020280    0.020360   0.196  0.84491    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.05048    0.06412     0.06442 109.447  < 2e-16 ***
cond(Con_Fire_2019)            -0.14368    0.04941     0.04963   2.895  0.00379 ** 
cond(Con_Rivers)               -0.07989    0.04756     0.04778   1.672  0.09456 .  
cond(Con_Fire_2019:Con_Rivers)  0.02097    0.04253     0.04273   0.491  0.62367
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
Model_Con_2018_nb<-glmmTMB(Bram_Abu~Con_Rivers+Con_Woodland+Con_Harvest_2018+Con_Fire_2018+
                             Con_Rivers:Con_Fire_2018+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2018_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.27         1.13      0.79
             Con_Woodland 1.30         1.14      0.77
         Con_Harvest_2018 1.26         1.12      0.80
            Con_Fire_2018 1.29         1.14      0.77
 Con_Rivers:Con_Fire_2018 1.09         1.04      0.92
'''
options(na.action="na.fail")
Model_Con_2018_nb_Dredge<-dredge(Model_Con_2018_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2018_nb_Subset<-subset(Model_Con_2018_nb_Dredge,delta<2)
Model_Con_2018_nb_Ave<-model.avg(Model_Con_2018_nb_Subset)
importance(Model_Con_2018_nb_Ave)
'''
                     cond(Con_Fire_2018) cond(Con_Rivers) cond(Con_Harvest_2018) cond(Con_Woodland)
Sum of weights:      1.00                1.00             0.67                   0.23              
N containing models:    3                   3                2                      1    
'''
confint(Model_Con_2018_nb_Ave)
'''
                             2.5 %      97.5 %
cond((Int))             6.91934692  7.14854670
cond(Con_Fire_2018)    -0.37559854 -0.17799548 *
cond(Con_Harvest_2018) -0.01417652  0.16717924
cond(Con_Rivers)       -0.22466945 -0.03248726 *
cond(Con_Woodland)     -0.15138493  0.05596822
'''
summary(Model_Con_2018_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             7.03395    0.05820     0.05847 120.299  < 2e-16 ***
cond(Con_Fire_2018)    -0.27680    0.05019     0.05041   5.491  < 2e-16 ***
cond(Con_Harvest_2018)  0.05119    0.05210     0.05223   0.980  0.32700    
cond(Con_Rivers)       -0.12858    0.04881     0.04903   2.623  0.00873 ** 
cond(Con_Woodland)     -0.01091    0.03218     0.03227   0.338  0.73529    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             7.03395    0.05820     0.05847 120.299  < 2e-16 ***
cond(Con_Fire_2018)    -0.27680    0.05019     0.05041   5.491  < 2e-16 ***
cond(Con_Harvest_2018)  0.07650    0.04605     0.04627   1.654  0.09822 .  
cond(Con_Rivers)       -0.12858    0.04881     0.04903   2.623  0.00873 ** 
cond(Con_Woodland)     -0.04771    0.05265     0.05290   0.902  0.36711
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
Model_Con_2017_nb<-glmmTMB(Bram_Abu~Con_Rivers+Con_Woodland+Con_Harvest_2017+Con_Fire_2017+
                             Con_Rivers:Con_Fire_2017+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2017_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.25         1.12      0.80
             Con_Woodland 1.39         1.18      0.72
         Con_Harvest_2017 1.29         1.14      0.77
            Con_Fire_2017 1.28         1.13      0.78
 Con_Rivers:Con_Fire_2017 1.10         1.05      0.91
'''
options(na.action="na.fail")
Model_Con_2017_nb_Dredge<-dredge(Model_Con_2017_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2017_nb_Subset<-subset(Model_Con_2017_nb_Dredge,delta<2)
Model_Con_2017_nb_Ave<-model.avg(Model_Con_2017_nb_Subset)
importance(Model_Con_2017_nb_Ave)
'''
                     cond(Con_Fire_2017) cond(Con_Rivers) cond(Con_Harvest_2017) cond(Con_Fire_2017:Con_Rivers)
Sum of weights:      1.00                0.86             0.67                   0.15                          
N containing models:    5                   4                3                      1                          
                     cond(Con_Woodland)
Sum of weights:      0.14              
N containing models:    1    
'''
confint(Model_Con_2017_nb_Ave)
'''
                                      2.5 %        97.5 %
cond((Int))                     6.915156241  7.1752150697
cond(Con_Fire_2017)            -0.252173355 -0.0572167901 *
cond(Con_Harvest_2017)         -0.006316717  0.1764073166
cond(Con_Rivers)               -0.202390120 -0.0006724493 *
cond(Con_Fire_2017:Con_Rivers) -0.074088706  0.1327425613
cond(Con_Woodland)             -0.134921000  0.0906466866
'''
summary(Model_Con_2017_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.045186   0.066034    0.066343 106.194  < 2e-16 ***
cond(Con_Fire_2017)            -0.154695   0.049508    0.049735   3.110  0.00187 ** 
cond(Con_Harvest_2017)          0.057380   0.055135    0.055258   1.038  0.29908    
cond(Con_Rivers)               -0.087138   0.059220    0.059387   1.467  0.14230    
cond(Con_Fire_2017:Con_Rivers)  0.004517   0.023171    0.023258   0.194  0.84599    
cond(Con_Woodland)             -0.003143   0.022921    0.023017   0.137  0.89140    
 
(conditional average) 
                               Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                     7.04519    0.06603     0.06634 106.194  < 2e-16 ***
cond(Con_Fire_2017)            -0.15470    0.04951     0.04973   3.110  0.00187 ** 
cond(Con_Harvest_2017)          0.08505    0.04640     0.04661   1.824  0.06808 .  
cond(Con_Rivers)               -0.10153    0.05123     0.05146   1.973  0.04849 *  
cond(Con_Fire_2017:Con_Rivers)  0.02933    0.05252     0.05276   0.556  0.57834    
cond(Con_Woodland)             -0.02214    0.05727     0.05754   0.385  0.70046
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
Model_Con_2016_nb<-glmmTMB(Bram_Abu~Con_Rivers+Con_Woodland+Con_Harvest_2016+Con_Fire_2016+
                             Con_Rivers:Con_Fire_2016+(1|PLANTATI_1),
                           data=df_Bram_Con_scale,family="nbinom2")
check_collinearity(Model_Con_2016_nb)
'''
Low Correlation

                     Term  VIF Increased SE Tolerance
               Con_Rivers 1.46         1.21      0.69
             Con_Woodland 1.36         1.16      0.74
         Con_Harvest_2016 1.40         1.19      0.71
            Con_Fire_2016 1.28         1.13      0.78
 Con_Rivers:Con_Fire_2016 1.26         1.12      0.80
'''
options(na.action="na.fail")
Model_Con_2016_nb_Dredge<-dredge(Model_Con_2016_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Con_2016_nb_Subset<-subset(Model_Con_2016_nb_Dredge,delta<2)
Model_Con_2016_nb_Ave<-model.avg(Model_Con_2016_nb_Subset)
importance(Model_Con_2016_nb_Ave)
'''
                     cond(Con_Fire_2016) cond(Con_Rivers) cond(Con_Harvest_2016) cond(Con_Woodland)
Sum of weights:      1.00                1.00             0.74                   0.22              
N containing models:    3                   3                2                      1
'''
confint(Model_Con_2016_nb_Ave)
'''
                              2.5 %      97.5 %
cond((Int))             6.931254371  7.15755460
cond(Con_Fire_2016)    -0.331590434 -0.13640715 *
cond(Con_Harvest_2016) -0.005826654  0.17680511 
cond(Con_Rivers)       -0.231396663 -0.02486902 *
cond(Con_Woodland)     -0.143006471  0.07452583
'''
summary(Model_Con_2016_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                        Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             7.044404   0.057461    0.057731 122.022  < 2e-16 ***
cond(Con_Fire_2016)    -0.233999   0.049566    0.049793   4.699  2.6e-06 ***
cond(Con_Harvest_2016)  0.063318   0.054743    0.054878   1.154    0.249    
cond(Con_Rivers)       -0.128133   0.052467    0.052687   2.432    0.015 *  
cond(Con_Woodland)     -0.007511   0.029495    0.029602   0.254    0.800    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             7.04440    0.05746     0.05773 122.022  < 2e-16 ***
cond(Con_Fire_2016)    -0.23400    0.04957     0.04979   4.699  2.6e-06 ***
cond(Con_Harvest_2016)  0.08549    0.04638     0.04659   1.835   0.0665 .  
cond(Con_Rivers)       -0.12813    0.05247     0.05269   2.432   0.0150 *  
cond(Con_Woodland)     -0.03424    0.05523     0.05549   0.617   0.5372 
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
a<-ggplot(data=df_Bram_Con,aes(x=Con_Rivers,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nDistance to rivers (m)",tag="a)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
b<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2020,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nFire severity 2020",tag="b)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
c<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2019,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nFire severity 2019",tag="c)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
d<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2018,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nFire severity 2018",tag="d)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
e<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2017,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nFire severity 2017",tag="e)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
f<-ggplot(data=df_Bram_Con,aes(x=Con_Fire_2016,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in grasslands\n",x="\nFire severity 2016",tag="f)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
ggarrange(a,b,c,d,e,f,ncol=3,nrow=2) # 14 x 8.5
rm(a,b,c,d,e,f,df_Bram_Con,df_Bram_Con_scale)

# Model building: Production areas
# 2020
Model_Pro_2020_nb<-glmmTMB(Bram_Abu~Pro_Rivers+Pro_Woodland+Pro_Harvest_2020+
                             Pro_Rivers:Pro_Harvest_2020+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2020_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.00         1.00      1.00
                Pro_Woodland 1.12         1.06      0.90
            Pro_Harvest_2020 1.13         1.06      0.89
 Pro_Rivers:Pro_Harvest_2020 1.01         1.01      0.99
'''
options(na.action="na.fail")
Model_Pro_2020_nb_Dredge<-dredge(Model_Pro_2020_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2020_nb_Subset<-subset(Model_Pro_2020_nb_Dredge,delta<2)
Model_Pro_2020_nb_Ave<-model.avg(Model_Pro_2020_nb_Subset)
importance(Model_Pro_2020_nb_Ave)
'''
                     cond(Pro_Harvest_2020) cond(Pro_Woodland) cond(Pro_Rivers)
Sum of weights:      1.00                   1.00               0.34            
N containing models:    2                      2                  1 
'''
confint(Model_Pro_2020_nb_Ave)
'''
                             2.5 %     97.5 %
cond((Int))             6.39960762  6.8976532
cond(Pro_Harvest_2020) -0.45146052 -0.2119503 *
cond(Pro_Woodland)      0.06920352  0.3487356 *
cond(Pro_Rivers)       -0.06318513  0.1685239
'''
summary(Model_Pro_2020_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.64863    0.12640     0.12705  52.329  < 2e-16 ***
cond(Pro_Harvest_2020) -0.33171    0.06078     0.06110   5.429    1e-07 ***
cond(Pro_Woodland)      0.20897    0.07094     0.07131   2.930  0.00339 ** 
cond(Pro_Rivers)        0.01814    0.04263     0.04277   0.424  0.67155    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.64863    0.12640     0.12705  52.329  < 2e-16 ***
cond(Pro_Harvest_2020) -0.33171    0.06078     0.06110   5.429    1e-07 ***
cond(Pro_Woodland)      0.20897    0.07094     0.07131   2.930  0.00339 ** 
cond(Pro_Rivers)        0.05267    0.05880     0.05911   0.891  0.37291 
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
Model_Pro_2019_nb<-glmmTMB(Bram_Abu~Pro_Rivers+Pro_Woodland+Pro_Harvest_2019+
                             Pro_Rivers:Pro_Harvest_2019+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2019_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.02         1.01      0.98
                Pro_Woodland 1.07         1.04      0.93
            Pro_Harvest_2019 1.07         1.03      0.93
 Pro_Rivers:Pro_Harvest_2019 1.02         1.01      0.98
'''
options(na.action="na.fail")
Model_Pro_2019_nb_Dredge<-dredge(Model_Pro_2019_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2019_nb_Subset<-subset(Model_Pro_2019_nb_Dredge,delta<2)
Model_Pro_2019_nb_Ave<-model.avg(Model_Pro_2019_nb_Subset)
importance(Model_Pro_2019_nb_Ave)
'''
                     cond(Pro_Harvest_2019) cond(Pro_Woodland) cond(Pro_Rivers) cond(Pro_Harvest_2019:Pro_Rivers)
Sum of weights:      1.00                   1.00               0.47             0.26                             
N containing models:    3                      3                  2                1   
'''
confint(Model_Pro_2019_nb_Ave)
'''
                                        2.5 %      97.5 %
cond((Int))                        6.54207569  6.86274279
cond(Pro_Harvest_2019)            -0.48650042 -0.29347308 *
cond(Pro_Woodland)                 0.08009658  0.34074733 *
cond(Pro_Rivers)                  -0.07805070  0.14378005
cond(Pro_Harvest_2019:Pro_Rivers) -0.18249657  0.01721941
'''
summary(Model_Pro_2019_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                                  Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        6.70241    0.08138     0.08180  81.932  < 2e-16 ***
cond(Pro_Harvest_2019)            -0.38999    0.04899     0.04924   7.920  < 2e-16 ***
cond(Pro_Woodland)                 0.21042    0.06615     0.06649   3.165  0.00155 ** 
cond(Pro_Rivers)                   0.01547    0.04197     0.04215   0.367  0.71356    
cond(Pro_Harvest_2019:Pro_Rivers) -0.02147    0.04450     0.04458   0.482  0.63016    
 
(conditional average) 
                                  Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))                        6.70241    0.08138     0.08180  81.932  < 2e-16 ***
cond(Pro_Harvest_2019)            -0.38999    0.04899     0.04924   7.920  < 2e-16 ***
cond(Pro_Woodland)                 0.21042    0.06615     0.06649   3.165  0.00155 ** 
cond(Pro_Rivers)                   0.03286    0.05630     0.05659   0.581  0.56141    
cond(Pro_Harvest_2019:Pro_Rivers) -0.08264    0.05068     0.05095   1.622  0.10481
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
Model_Pro_2018_nb<-glmmTMB(Bram_Abu~Pro_Rivers+Pro_Woodland+Pro_Harvest_2018+
                             Pro_Rivers:Pro_Harvest_2018+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2018_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.01         1.00      0.99
                Pro_Woodland 1.02         1.01      0.98
            Pro_Harvest_2018 1.02         1.01      0.98
 Pro_Rivers:Pro_Harvest_2018 1.00         1.00      1.00
'''
options(na.action="na.fail")
Model_Pro_2018_nb_Dredge<-dredge(Model_Pro_2018_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2018_nb_Subset<-subset(Model_Pro_2018_nb_Dredge,delta<2)
Model_Pro_2018_nb_Ave<-model.avg(Model_Pro_2018_nb_Subset)
importance(Model_Pro_2018_nb_Ave)
'''
                     cond(Pro_Harvest_2018) cond(Pro_Woodland) cond(Pro_Rivers)
Sum of weights:      1.00                   1.00               0.27            
N containing models:    2                      2                  1
'''
confint(Model_Pro_2018_nb_Ave)
'''
                             2.5 %     97.5 %
cond((Int))             6.50474559  6.8901112
cond(Pro_Harvest_2018) -0.40627716 -0.1977936 *
cond(Pro_Woodland)      0.14792003  0.4174112 *
cond(Pro_Rivers)       -0.09537794  0.1361539
'''
summary(Model_Pro_2018_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                        Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.697428   0.097800    0.098309  68.126  < 2e-16 ***
cond(Pro_Harvest_2018) -0.302035   0.052910    0.053186   5.679  < 2e-16 ***
cond(Pro_Woodland)      0.282666   0.068393    0.068749   4.112 3.93e-05 ***
cond(Pro_Rivers)        0.005517   0.031879    0.032032   0.172    0.863    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.69743    0.09780     0.09831  68.126  < 2e-16 ***
cond(Pro_Harvest_2018) -0.30204    0.05291     0.05319   5.679  < 2e-16 ***
cond(Pro_Woodland)      0.28267    0.06839     0.06875   4.112 3.93e-05 ***
cond(Pro_Rivers)        0.02039    0.05876     0.05907   0.345     0.73 
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
Model_Pro_2017_nb<-glmmTMB(Bram_Abu~Pro_Rivers+Pro_Woodland+Pro_Harvest_2017+
                             Pro_Rivers:Pro_Harvest_2017+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2017_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.00         1.00      1.00
                Pro_Woodland 1.02         1.01      0.98
            Pro_Harvest_2017 1.02         1.01      0.98
 Pro_Rivers:Pro_Harvest_2017 1.02         1.01      0.98
'''
options(na.action="na.fail")
Model_Pro_2017_nb_Dredge<-dredge(Model_Pro_2017_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2017_nb_Subset<-subset(Model_Pro_2017_nb_Dredge,delta<2)
Model_Pro_2017_nb_Ave<-model.avg(Model_Pro_2017_nb_Subset)
importance(Model_Pro_2017_nb_Ave)
'''
                     cond(Pro_Woodland) cond(Pro_Harvest_2017) cond(Pro_Rivers)
Sum of weights:      1.00               0.72                   0.25            
N containing models:    3                  2                      1 
'''
confint(Model_Pro_2017_nb_Ave)
'''
                             2.5 %      97.5 %
cond((Int))             6.50005914 6.939657545
cond(Pro_Harvest_2017) -0.22274726 0.007878911
cond(Pro_Woodland)      0.17029017 0.462603003 *
cond(Pro_Rivers)       -0.06534086 0.177372742
'''
summary(Model_Pro_2017_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.71986    0.11157     0.11214  59.921  < 2e-16 ***
cond(Pro_Harvest_2017) -0.07736    0.06923     0.06942   1.114    0.265    
cond(Pro_Woodland)      0.31645    0.07419     0.07457   4.244  2.2e-05 ***
cond(Pro_Rivers)        0.01400    0.03920     0.03933   0.356    0.722    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.71986    0.11157     0.11214  59.921  < 2e-16 ***
cond(Pro_Harvest_2017) -0.10743    0.05853     0.05883   1.826   0.0678 .  
cond(Pro_Woodland)      0.31645    0.07419     0.07457   4.244  2.2e-05 ***
cond(Pro_Rivers)        0.05602    0.06160     0.06192   0.905   0.3656 
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
Model_Pro_2016_nb<-glmmTMB(Bram_Abu~Pro_Rivers+Pro_Woodland+Pro_Harvest_2016+
                             Pro_Rivers:Pro_Harvest_2016+(1|PLANTATI_1),
                           data=df_Bram_Pro_scale,family="nbinom2")
check_collinearity(Model_Pro_2016_nb)
'''
Low Correlation

                        Term  VIF Increased SE Tolerance
                  Pro_Rivers 1.02         1.01      0.98
                Pro_Woodland 1.01         1.01      0.99
            Pro_Harvest_2016 1.01         1.00      0.99
 Pro_Rivers:Pro_Harvest_2016 1.03         1.01      0.97
'''
options(na.action="na.fail")
Model_Pro_2016_nb_Dredge<-dredge(Model_Pro_2016_nb,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Pro_2016_nb_Subset<-subset(Model_Pro_2016_nb_Dredge,delta<2)
Model_Pro_2016_nb_Ave<-model.avg(Model_Pro_2016_nb_Subset)
importance(Model_Pro_2016_nb_Ave)
'''
                     cond(Pro_Woodland) cond(Pro_Rivers) cond(Pro_Harvest_2016)
Sum of weights:      1.00               0.25             0.24                  
N containing models:    3                  1                1
'''
confint(Model_Pro_2016_nb_Ave)
'''
                             2.5 %    97.5 %
cond((Int))             6.50840071 6.9444387
cond(Pro_Woodland)      0.18053424 0.4730813 *
cond(Pro_Rivers)       -0.07190741 0.1721730
cond(Pro_Harvest_2016) -0.07676976 0.1697576
'''
summary(Model_Pro_2016_nb_Ave)
'''
Model-averaged coefficients:  
(full average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.72642    0.11066     0.11124  60.470  < 2e-16 ***
cond(Pro_Woodland)      0.32681    0.07425     0.07463   4.379 1.19e-05 ***
cond(Pro_Rivers)        0.01258    0.03788     0.03801   0.331    0.741    
cond(Pro_Harvest_2016)  0.01109    0.03642     0.03655   0.303    0.762    
 
(conditional average) 
                       Estimate Std. Error Adjusted SE z value Pr(>|z|)    
cond((Int))             6.72642    0.11066     0.11124  60.470  < 2e-16 ***
cond(Pro_Woodland)      0.32681    0.07425     0.07463   4.379 1.19e-05 ***
cond(Pro_Rivers)        0.05013    0.06194     0.06227   0.805    0.421    
cond(Pro_Harvest_2016)  0.04649    0.06257     0.06289   0.739    0.460 
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
a<-ggplot(data=df_Bram_Pro,aes(x=Pro_Woodland,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area\n",x="\nDistance to woodlands (m)",tag="a)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
b<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2020,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area\n",x="\nDistance to harvested trees 2020 (m)",tag="b)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
c<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2019,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area\n",x="\nDistance to harvested trees 2019 (m)",tag="c)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
d<-ggplot(data=df_Bram_Pro,aes(x=Pro_Harvest_2018,y=Bram_Abu))+
  geom_point(colour="black",shape=21,size=2,fill="light blue")+
  geom_smooth(method="lm",colour="black",se=FALSE)+
  labs(y="Bramble abundance in production area\n",x="\nDistance to harvested trees 2018 (m)",tag="d)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.tag=element_text(hjust=0.5))
ggarrange(a,b,c,d,ncol=2,nrow=2) # 14 x 8.5
rm(a,b,c,d,df_Bram_Pro,df_Bram_Pro_scale)

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

# Outliers
dotchart(Veg$Tree_Rich,ylab="Tree_Rich")
dotchart(Veg$Tree_Hei,ylab="Tree_Hei")
dotchart(Veg$Shrub_Rich,ylab="Shrub_Rich")
dotchart(Veg$Shrub_Hei,ylab="Shrub_Hei")
dotchart(Veg$Suc_Rich,ylab="Suc_Rich")
dotchart(Veg$Suc_Hei,ylab="Suc_Hei")
dotchart(Veg$Fern_Rich,ylab="Fern_Rich")
dotchart(Veg$Fern_Hei,ylab="Fern_Hei")
dotchart(Veg$Sed_Rich,ylab="Sed_Rich")
dotchart(Veg$Sed_Hei,ylab="Sed_Hei")
dotchart(Veg$Forb_Rich,ylab="Forb_Rich")
dotchart(Veg$Forb_Hei,ylab="Forb_Hei")
dotchart(Veg$Bulb_Rich,ylab="Bulb_Rich")
dotchart(Veg$Bulb_Hei,ylab="Bulb_Hei")
dotchart(Veg$Grass_Rich,ylab="Grass_Rich")
dotchart(Veg$Grass_Hei,ylab="Grass_Hei")
dotchart(Veg$Bram_Abu,ylab="Bram_Abu")
dotchart(Veg$Bram_Cov,ylab="Bram_Cov")
dotchart(Veg$Grou_Cov,ylab="Grou_Cov")
dotchart(Veg$Rock_Cov,ylab="Rock_Cov")
dotchart(Veg$Climb_Rich,ylab="Climb_Rich")
dotchart(Veg$Climb_Hei,ylab="Climb_Hei")

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
                        data=Veg_Scale,family="poisson")
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
ggplot(data=Veg,aes(x=Bram_Abu,y=Veg_Rich))+
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

# Modelling Grasshopper richness Bram_Cov
Model_hop_Rich<-glmmTMB(hop_exShan~Bram_Abu+(1|Plantation),
                        data=df_Scale,family="gaussian")
confint(Model_hop_Rich)
'''
                                         2.5 %    97.5 %   Estimate
cond.(Intercept)                     8.1758899 10.651599  9.4137446
cond.Bram_Abu                       -0.8537046  0.643342 -0.1051813
Plantation.cond.Std.Dev.(Intercept)  0.3569254  2.872029  1.0124723
'''
summary(Model_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   9.4137     0.6316  14.905   <2e-16 ***
Bram_Abu     -0.1052     0.3819  -0.275    0.783 
'''
r.squaredGLMM(Model_hop_Rich)
'''
             R2m       R2c
[1,] 0.001413845 0.1324199
'''
# Model assumption
plot(residuals(Model_hop_Rich))
plot(simulateResiduals(Model_hop_Rich))
testDispersion(simulateResiduals(Model_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_hop_Rich)

# Modelling Low value grasshopper species Bram_Abu
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

# Modelling Low value grasshopper species Bram_Cov
Model_Low_hop_Rich<-glmmTMB(Low_exShan~Bram_Abu+(1|Plantation),
                            data=df_Scale,family="gaussian")
confint(Model_Low_hop_Rich)
'''
                                         2.5 %    97.5 %  Estimate
cond.(Intercept)                     5.0409256 6.6815205 5.8612230
cond.Bram_Abu                       -0.4398466 0.6422458 0.1011996
Plantation.cond.Std.Dev.(Intercept)  0.1968465 2.0269204 0.6316583
'''
summary(Model_Low_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   5.8612     0.4185  14.004   <2e-16 ***
Bram_Abu      0.1012     0.2760   0.367    0.714 
'''
r.squaredGLMM(Model_Low_hop_Rich)
'''
             R2m       R2c
[1,] 0.002559794 0.1022865
'''
# Model assumption
plot(residuals(Model_Low_hop_Rich))
plot(simulateResiduals(Model_Low_hop_Rich))
testDispersion(simulateResiduals(Model_Low_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Low_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Low_hop_Rich)

# Modelling Intermediate value grasshopper species Bram_Abu
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
ggplot(data=df,aes(x=Bram_Abu,y=Intermediate_hop_Rich))+
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

# Modelling Intermediate value grasshopper species Bram_Cov
Model_Intermediate_hop_Rich<-glmmTMB(Intermediate_exShan~Bram_Abu+(1|Plantation),
                                     data=df_Scale,family="gaussian")
confint(Model_Intermediate_hop_Rich)
'''
                                            2.5 %       97.5 %    Estimate
cond.(Intercept)                     3.267233e+00 3.963644e+00  3.61543847
cond.Bram_Abu                       -6.654051e-01 4.049755e-02 -0.31245379
Plantation.cond.Std.Dev.(Intercept)  1.041788e-07 6.272423e+04  0.08083647
'''
summary(Model_Intermediate_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   3.6154     0.1777  20.350   <2e-16 ***
Bram_Abu     -0.3125     0.1801  -1.735   0.0827 .
'''
r.squaredGLMM(Model_Intermediate_hop_Rich)
'''
            R2m        R2c
[1,] 0.06119043 0.06528611
'''
# Model assumption
plot(residuals(Model_Intermediate_hop_Rich))
plot(simulateResiduals(Model_Intermediate_hop_Rich))
testDispersion(simulateResiduals(Model_Intermediate_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Intermediate_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Intermediate_hop_Rich)

# Modelling High value grasshopper species Bram_Abu
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

# Modelling High value grasshopper species
Model_High_hop_Rich<-glmmTMB(High_exShan~Bram_Abu+(1|Plantation),
                             data=df_Scale,family="poisson")
confint(Model_High_hop_Rich)
'''
                                          2.5 %    97.5 %     Estimate
cond.(Intercept)                    -0.06136097 0.4382118 1.884254e-01
cond.Bram_Abu                       -0.23374444 0.2625332 1.439437e-02
Plantation.cond.Std.Dev.(Intercept)  0.00000000       Inf 1.570107e-05
'''
summary(Model_High_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  0.18843    0.12744   1.478    0.139
Bram_Abu     0.01439    0.12660   0.114    0.909
'''
r.squaredGLMM(Model_High_hop_Rich)
'''
                   R2m          R2c
delta     0.0002501228 0.0002501231
lognormal 0.0003433119 0.0003433123
trigamma  0.0001648835 0.0001648837
'''
# Model assumption
plot(residuals(Model_High_hop_Rich))
plot(simulateResiduals(Model_High_hop_Rich))
testDispersion(simulateResiduals(Model_High_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_High_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_High_hop_Rich)

# Modelling Grasshopper richness Bram_Abu
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

# Modelling Grasshopper richness Bram_Cov
Model_Caelifera_hop_Rich<-glmmTMB(Caelifera_exShan~Bram_Abu+(1|Plantation),
                                  data=df_Scale,family="gaussian")
confint(Model_Caelifera_hop_Rich)
'''
                                         2.5 %    97.5 %   Estimate
cond.(Intercept)                     6.8297513 8.9502431 7.88999721
cond.Bram_Abu                       -0.6314177 0.7032967 0.03593949
Plantation.cond.Std.Dev.(Intercept)  0.2819053 2.5806181 0.85293014
'''
summary(Model_Caelifera_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  7.89000    0.54095  14.585   <2e-16 ***
Bram_Abu     0.03594    0.34049   0.106    0.916
'''
r.squaredGLMM(Model_Caelifera_hop_Rich)
'''
              R2m       R2c
[1,] 0.0002134731 0.1204472
'''
# Model assumption
plot(residuals(Model_Caelifera_hop_Rich))
plot(simulateResiduals(Model_Caelifera_hop_Rich))
testDispersion(simulateResiduals(Model_Caelifera_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Caelifera_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)
rm(Model_Caelifera_hop_Rich)

# Modelling Grasshopper richness Bram_Abu
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

# Modelling Grasshopper richness Bram_Cov
Model_Ensifera_hop_Rich<-glmmTMB(Ensifera_exShan~Bram_Abu+(1|Plantation),
                                 data=df_Scale,family="poisson")
confint(Model_Ensifera_hop_Rich)
'''
                                         2.5 %     97.5 %      Estimate
cond.(Intercept)                     0.5873385 0.96134763  7.743431e-01
cond.Bram_Abu                       -0.3295511 0.09834826 -1.156014e-01
Plantation.cond.Std.Dev.(Intercept)  0.0000000        Inf  1.624334e-05
'''
summary(Model_Ensifera_hop_Rich)
'''
Conditional model:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  0.77434    0.09541   8.116 4.83e-16 ***
Bram_Abu    -0.11560    0.10916  -1.059     0.29
'''
r.squaredGLMM(Model_Ensifera_hop_Rich)
'''
                 R2m        R2c
delta     0.02833648 0.02833648
lognormal 0.03421355 0.03421355
trigamma  0.02257209 0.02257209
'''
# Model assumption
plot(residuals(Model_Ensifera_hop_Rich))
plot(simulateResiduals(Model_Ensifera_hop_Rich))
testDispersion(simulateResiduals(Model_Ensifera_hop_Rich))
testSpatialAutocorrelation(simulateResiduals(Model_Ensifera_hop_Rich),
                           x=df_Scale$Long_X,y=df_Scale$Lat_Y)

# Plot
ggarrange(a,b,ncol=2,nrow=1) # 8x4
rm(a,b,Model_Ensifera_hop_Rich,df_Scale,df)


###########################################################################
# How will uncontrolled bramble invasion impact the landscape?
###########################################################################

#### Grasshopper distribution modeling ####
# Load and clip variables
Study_ROI<-readOGR("Shapefiles/Plantations_Clipped.shp")
NDVI<-stack("Products/Cliped_SR_S2.tif",bands=c(3,4))
names(NDVI)<-c("Red","Infrared")
NDVI<-(NDVI[["Infrared"]]-NDVI[["Red"]])/(NDVI[["Infrared"]]+NDVI[["Red"]])
NDVI[NDVI>1]<-NA
NDVI[NDVI<(-1)]<-NA
Land_Use<-raster("Products/Map_ranger_SR_S2.tif")
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
NDVI<-resample(NDVI,template_rst)
Land_Use<-resample(Land_Use,template_rst,method="ngb")
Rivers<-resample(Rivers,template_rst)
Elv<-resample(Elv,template_rst)
Aspect<-resample(Aspect,template_rst)
Plantation_Edge<-resample(Plantation_Edge,template_rst)
Fire_hist<-resample(Fire_hist,template_rst)
rm(template_rst)

# Clip and save rasters
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

# Load environmental raster data
Env<-load_var(path="~/Products/",
              files=c("SDM_NDVI.tif","SDM_Land_Use.tif","SDM_Rivers.tif","SDM_Aspect.tif",
                      "SDM_Plantation_Edge.tif","SDM_Elv.tif","SDM_Fire_hist.tif"),
              format=".tif",Norm=TRUE,categorical="SDM_Land_Use.tif")
# Load occurrence data for ALL SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="All_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/All_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/All_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/All_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM)#5x4

# Load occurrence data for LOW SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_Low_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="Low_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/Low_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/Low_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/Low_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM)#5x4

# Load occurrence data for Intermediate SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_Intermediate_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="Intermediate_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/Intermediate_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/Intermediate_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/Intermediate_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM)#5x4

# Load occurrence data for High SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_High_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="High_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/High_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/High_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/High_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM)#5x4

# Load occurrence data for Caelifera SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_Caelifera_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="Caelifera_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/Caelifera_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/Caelifera_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/Caelifera_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM)#5x4

# Load occurrence data for Ensifera SPECIES
Occ<-load_occ(path="~/Products/",Env=Env,
              file="Occurrence_Ensifera_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")
# Perform stacked species distribution models 
SDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                     Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",metric="SES",rep=10,tmp=TRUE,
                     save=TRUE,name="Ensifera_Species_SDM",path="~/Products/",
                     cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))
# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/Ensifera_Species_SDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Dist Riv","Aspect",
                                  "Dist Plan","Elv","Fire Hist")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(variable.importance) #7x5
# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/Ensifera_Species_SDM/Stack/Tables/StackEval.csv",row.names=1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(evaluation) #7x5
# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/Ensifera_Species_SDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(correlation,Occ,SDM,Env)#5x4

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