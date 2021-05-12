#make map of residential address on a national map
library(maptools)
library(sp)
library(rgdal)
library(maps)
library(sp)
library(ggplot2)
library(gstat) #use the idw() function in gstat, mask spatstat
library(maptools)
library(raster)
library(mapproj)
library(data.table)
mapUSm<- readShapePoly('/Users/cindyhu/Documents/Research/EWG/moveRC/state_shape/US_48states')
proj4string(mapUSm)<-CRS("+proj=longlat +datum=WGS84")

xy<-data.frame(ID=df$id,X=df$gdtlong1990, Y=df$gdtlat1990)
coordinates(xy)<-c('X','Y')
proj4string(xy)<-CRS('+proj=longlat +datum=WGS84')

#the 5 locations in 2016
idList<-c('108067','134102','126013','112072','131456')

postscript("NHS_225_map_030718.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 8, width = 10)
ggplot(mapUSm) +
  geom_polygon(aes(x = long, y = lat, group = group),colour='white', fill='lightblue',size=0.2)+
  coord_equal()+coord_map("lambert", parameters = c(20, 50))+
  geom_point(data=df[is.na(df$batch_pfas),],aes(x=gdtlong1990, y=gdtlat1990),colour='black',pch=21,size=2)+
  geom_point(data=df[!is.na(df$batch_pfas),],aes(x=gdtlong1990, y=gdtlat1990),colour='darkorange',pch=19,size=2)+
  geom_point(data=df[df$id%in%idList,],aes(x=gdtlong1990, y=gdtlat1990),colour='forestgreen',pch=20,size=1.5)+
  theme(panel.background = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())
dev.off()
