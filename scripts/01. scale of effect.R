# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############## Landscape analysis for Palheta et al. (2023) ###############
# # # # # # # # # # # # # # By: Lucas Colares # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # 20-08-2023 # # # # # # # ## # # # # # # # #
library(raster)
library(sf)
raster("datasets/rasters/brasil_coverage_2019.tif")->brazil2019
raster("datasets/rasters/brasil_coverage_2018.tif")->brazil2018
st_read("datasets/shapefiles/PA_Municipios_2022.shp")->PAshp
read.csv("datasets/coordinates.csv",sep=";",row.names = 1)->coords
coords[match(rownames(HII),rownames(coords)),]->coords
coords$IIH<-HII$IIH_total
write.csv(coords,"datasets/coordinates.csv",fileEncoding = "latin1",row.names = F)

raster::crop(brazil2019, raster::extent(PAshp))->PAcov2019
raster::crop(brazil2018, raster::extent(PAshp))->PAcov2018
rm(brazil2019,brazil2018)

spdf <- SpatialPointsDataFrame(coords = coords, data = coords, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

Buffers<-{}
for(x in c(seq(30,90,10),seq(100,1000,100))){
  st_buffer(st_as_sf(spdf),dist = x)->Buffers[[length(Buffers)+1]]
}
names(Buffers)=c(seq(30,90,10),seq(100,1000,100))

raster::crop(PAcov2019, extent(Buffers[[length(Buffers)]]))->FinalCov2019
raster::crop(PAcov2018, extent(Buffers[[length(Buffers)]]))->FinalCov2018
rm(PAcov2019,PAcov2018)
gc()

AllValues2019<-{}
for(x in 1:length(Buffers)){
  for(y in 1:nrow(Buffers[[x]])){
    mask2019 = raster::mask(FinalCov2019, Buffers[[x]][y,])
    Values2019<-raster::getValues(mask2019)
    gsub(paste0("^",c(19,39,20,40,62,41,36,46,47,48),"$",collapse = "|"),18,Values2019)->Values2019
    data.frame(table(na.omit(Values2019)))->Freq
    colnames(Freq)[1]="Cat"
    data.frame(Freq,Site=rownames(Buffers[[x]][y,]),Buffer=names(Buffers)[x])->AllValues2019[[length(AllValues2019)+1]]
  }
}
do.call("rbind",AllValues2019)->AllValues2019

AllValues2019[!AllValues2019$Cat==33,]->AllValues2019
data.frame(matrix(NA,nrow=length(unique(AllValues2018$Site)),ncol=length(unique(c(unique(AllValues2019$Cat),unique(AllValues2018$Cat))))))->DataFrameBase
colnames(DataFrameBase)=unique(c(unique(AllValues2019$Cat),unique(AllValues2018$Cat)))
rownames(DataFrameBase)=unique(AllValues2019$Site)

Final2019<-{}
for(x in names(Buffers)){
  AllValues2019[AllValues2019$Buffer==x,]->SelBuf
  NewData=DataFrameBase
  for(u in rownames(NewData)){
    SelBuf[SelBuf$Site==u,]->SelSite
    SelSite$Freq/sum(SelSite$Freq)->SelSite$Freq
    NewData[rownames(NewData)==u,][,match(SelSite$Cat,colnames(NewData))]<-SelSite$Freq
  }
  Final2019[[length(Final2019)+1]]<-NewData
}
names(Final2019)=names(Buffers)

for(x in 1:length(Final2019)){
  Final2019[[x]][is.na(Final2019[[x]])]<-0
}

### Select scale of effect
Adjs<-{}
for(x in 1:length(Final2019)){
  lm(HII$IIH_total~.,data=Final2019[[x]])->LMIHH
  step(LMIHH,direction = "both")->LMIHH
  summary(LMIHH)
  data.frame(Buffer=names(Final2019)[x],RAdj=RsquareAdj(LMIHH)[[2]])->Adjs[[length(Adjs)+1]]
}
do.call("rbind",Adjs)->Adjs
Adjs[order(Adjs$RAdj,decreasing = T),]

rowSums(decostand(comm[,1:ncol(comm)-1],"pa"))->S
rowSums(comm[,1:ncol(comm)-1])->Dens

Adjs_Rich<-{}
for(x in 1:length(Final2019)){
  glm(S~.,data=Final2019[[x]],family = "poisson")->LMIHH
  step(LMIHH,direction = "both")->LMIHH
  summary(LMIHH)
  data.frame(Buffer=names(Final2019)[x],RAdj=with(summary(LMIHH), 1 - deviance/null.deviance)
)->Adjs_Rich[[length(Adjs_Rich)+1]]
}
do.call("rbind",Adjs_Rich)->Adjs_Rich
Adjs_Rich[order(Adjs_Rich$RAdj,decreasing = T),]

Adjs_Dens<-{}
for(x in 1:length(Final2019)){
  lm(log10(Dens)~.,data=Final2019[[x]])->LMIHH
  step(LMIHH,direction = "forward")->LMIHH
  summary(LMIHH)
  data.frame(Buffer=names(Final2019)[x],RAdj=RsquareAdj(LMIHH)[[2]])->Adjs_Dens[[length(Adjs_Dens)+1]]
}
do.call("rbind",Adjs_Dens)->Adjs_Dens
Adjs_Dens[order(Adjs_Dens$RAdj,decreasing = T),]
