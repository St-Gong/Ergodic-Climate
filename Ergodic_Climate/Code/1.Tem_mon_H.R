# 1.Clean up workspace and load or install necessary packages---------
  rm(list=ls())
  want <- c("rasterVis","rgdal","sf","raster","RColorBrewer","sp","maps","maptools","rgeos",
            "utils","ncdf4","tmaptools","rworldmap","Matrix","ggplot2","ggpubr")
  need <- want[!(want %in% installed.packages()[,"Package"])]
  if (length(need)) install.packages(need)
  lapply(want, function(i) require(i, character.only=TRUE))
  rm(want, need)

# 2.Working directories--------
  dir <- list()
  dir$root <- dirname(getwd())
  dir$Temperature_monthly <- paste(dir$root,"/Data/Raw_Data/Temperature_monthly",sep="")
  dir$map <- paste(dir$root,"/Data/Raw_Data/Chinese_Map",sep="")
  dir$output <- paste(dir$root,"/Output",sep="")
  dir$elevation <- paste(dir$root,"/Data/Raw_Data/Elevation",sep="")

# 3.Download Monthly data from "http://hydrology.princeton.edu/data/pgf/v2/0.5deg/monthly/" -----

  #Monthly data
  url <- "http://hydrology.princeton.edu/data/pgf/v2/0.5deg/monthly/tas_monthly_1901-2012.nc"
  destfile <- paste(dir$Temperature_monthly,"/tas_monthly_1901-2012.nc",sep="")
  download.file(url, destfile=destfile, method="curl", quiet = FALSE, mode = "w")

  #Elevation data
  url <- "http://hydrology.princeton.edu/data/pgf/elevation_0.5deg.nc"
  destfile <- paste(dir$elevation,"/elevation_0.5deg.nc",sep="")
  download.file(url, destfile=destfile, method="curl", quiet = FALSE, mode = "w")
  
# 4.Read Data and Manipulation--------
  flist <- list.files(dir$Temperature_monthly, full.names=T)
  Temp_mon <- brick(flist[1])

  cn <- readOGR(paste(dir$map,"/Province.shp",sep=""))
  
  elevation <- raster(paste(dir$elevation,"/elevation_0.5deg.nc",sep=""))
  land <- elevation
  land[] <- ifelse(is.na(elevation[]), 0, 1)  
  id <- land
  id[] <- 1:length(id)
  id[land[]==1] <- NA
  location <- id[!is.na(id)]
  Temp_mon[][location,] <- NA
  
  Temp_mon[] <- Temp_mon[]-273.15
  e <-c(65,137,10,60)
  Temp_mon_cn <- crop(Temp_mon,e)
  cn <- crop(cn,e)
  
  name <- c()
  k <- 1
  for (i in 1901:2012){
    for (j in 1:12){
      name[k] <- paste(i, "_", j, sep="") 
      k <- k+1
    }
  }
  names(Temp_mon_cn) <- name


# 5.Function we need----
  # This part includes two function:
  # 1) A simple function used to shift a column of dataframe upward (To create Xt and Xt-1 when using regression)
  # 2) Function to calculate H
  ## 5.1.Function to shift a column upward-------
  shift <- function(x, n){
    c(x[-(seq(n))], rep(NA, n))
  }

  ## 5.2.Function to calculate H--------
  Calculate_H <- function(tmp){
    data <- cbind(tmp,tmp)
    colnames(data) = c("tmp1", "tmp2") 
    data <- as.data.frame(data)
    data$tmp2 <- shift(data$tmp1, 1) 
    reg<-lm(data$tmp2~data$tmp1,data=data)  
    coe <- summary(reg)$coefficients
    beta <- coe[2,1]
    alpha <- coe[1,1]
    a <- 1-beta
    sigma <- 0
    for (i in 1:112){
      t <- (1-a)^(2*i) 
      sigma <- sigma + t
    }
    H <- 1/2*(log(1+sigma)/log(67))
    return(H)
  }


# 6.Main Code-------
  H_list <- c()
  tmp_list <- c()
  month <- c()
  for (j in 1:12){
    for (i in 1:112){
      month[i] <- paste("X", i+1900, "_", j, sep="") 
  }
    Temp_month<- Temp_mon_cn[][,month]
    for (i in 1:dim(Temp_month[])[1]){
      tmp <- Temp_month[][i,]
      if (sum(!is.na(tmp))>=2){
        tmp_list[i] <- Calculate_H(tmp)
      }
      else {
        tmp_list[i] <- NA
      }
    }
    H_list <- cbind(H_list,tmp_list)
  }

  Temp_mon_cn[][,seq(1,12)] <- H_list

# 8.Draw The map ----- 
  png(paste(dir$output,"/Tem_H.png",sep=""), width=1000*4, height = 1000*3, res=80*4)
  
  for (i in 1:12){
    num <- paste("p",i,sep="")
    plot_data <- as.data.frame(Temp_mon_cn[[i]],xy=TRUE)
    colnames(plot_data) <- c("x","y","value")
    p<-  ggplot()+
        geom_raster(data=plot_data,aes(x=x,y=y,fill=value))+
        scale_fill_gradientn(colours=c(brewer.pal(9,'OrRd')[1],brewer.pal(9,'OrRd')[6]), na.value="white") +
        theme_classic()+
        theme(axis.line = element_blank(),
        legend.key.width=unit(3,"cm"),
        legend.title=element_blank())+
        xlab(NULL)+
        ylab(NULL)+
        geom_polygon(data=cn,aes(long, lat, group = group, fill = hole),colour = "black", fill = NA, size=0.1)
   
  assign(num,p)
    }


  ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, nrow = 3, ncol=4, common.legend = TRUE,labels = month.abb[1:12],legend = "bottom",vjust=0.3, font.label = list(size = 10))+
    theme(plot.margin = margin(1,0.5,1,0.5, "cm")) 

  dev.off()
 
