#Code for structuring and plotting XRD data
#accompanying code to Wichern et al. paper.

#load required packages or install using install.packages()
library(tidyverse)
library(rxylib)
library(powdR)
library(baseline)

#set working directory
setwd("insert own wd")

#-------------------Formatting functions to structure xrdml data----------------------------
formatfunc <- function(x){
  data.frame(t(sapply(x$dataset,c)))
}

formatfunc2 <- function(x){
  unlist(x$data_block)
}

formatfunc3 <- function(x){
  x[2351:4700]
}




#-------------------Structuring and formatting each sample------------------------
#W2c_19:

#extract 2theta column, same for all runs of this sample
W2c_19_C1<-read_xyData(file="W2c-19_C01.xrdml", metaData=TRUE)
W2c_19_C1_unlist<-data.frame(t(sapply(W2c_19_C1$dataset,c)))
W2c_19_C1_xy<-W2c_19_C1_unlist$data_block
W2c_19_C1_xy_fin<-as.data.frame(W2c_19_C1_xy)
x2th_W2c_19<-W2c_19_C1_xy_fin$X2Theta

#load .xrdml data with rxylib
data_frame_names_2c19 <- list.files(pattern = "W2c-19_C*") 
data_frame_names_2c19 
data_frame_list_2c19 <- lapply(data_frame_names_2c19, read_xyData)  # Read all data frames in working directory
data_frame_list_2c19   

#Apply formatting functions
res_2c19 <- lapply(data_frame_list_2c19, formatfunc)
resxy_2c19 <-lapply(res_2c19, formatfunc2)
resxy2_2c19 <-lapply(resxy_2c19, formatfunc3)

resxy3_2c19<-Reduce("+", resxy2_2c19)

W2c_19_full<-cbind(x2th_W2c_19,resxy3_2c19)
colnames(W2c_19_full)<-c("x2th", "counts")
W2c_19_full<-as.data.frame(W2c_19_full)

#Remove background:
W2c_19_full_bg<-bkg(xrd=W2c_19_full)
W2c_19_full_bg<-cbind(W2c_19_full_bg$tth, W2c_19_full_bg$counts-W2c_19_full_bg$background) #convoluted, want it in a dataframe
W2c_19_full_bg<-as.data.frame(W2c_19_full_bg)
colnames(W2c_19_full_bg)<-c("x2th", "counts")
#Final dataframe with x2theta and counts columns




#W5-3:

#extract 2theta column, same for all runs of this sample
W5_3_C1<-read_xyData(file="Ws3_C01.xrdml", metaData=TRUE)
W5_3_C1_unlist<-data.frame(t(sapply(W5_3_C1$dataset,c)))
W5_3_C1_xy<-W5_3_C1_unlist$data_block
W5_3_C1_xy_fin<-as.data.frame(W5_3_C1_xy)
x2th_W5_3<-W5_3_C1_xy_fin$X2Theta

#load .xrdml data with rxylib
data_frame_names_5_3 <- list.files(pattern = "Ws3_C*") 
data_frame_names_5_3
data_frame_list_5_3 <- lapply(data_frame_names_5_3, read_xyData)  # Read all data frames within working directory
data_frame_list_5_3   

#Apply formatting functions
res_5_3 <- lapply(data_frame_list_5_3, formatfunc)
resxy_5_3 <-lapply(res_5_3, formatfunc2)
resxy2_5_3 <-lapply(resxy_5_3, formatfunc3)

resxy3_5_3<-Reduce("+", resxy2_5_3)

W5_3_full<-cbind(x2th_W5_3,resxy3_5_3)
colnames(W5_3_full)<-c("x2th", "counts")
W5_3_full<-as.data.frame(W5_3_full)

#Remove background
W5_3_full_bg<-bkg(xrd=W5_3_full)
W5_3_full_bg<-cbind(W5_3_full_bg$tth, W5_3_full_bg$counts-W5_3_full_bg$background) #convoluted, want it in a dataframe
W5_3_full_bg<-as.data.frame(W5_3_full_bg)
colnames(W5_3_full_bg)<-c("x2th", "counts")
#Final dataframe with x2theta and counts columns




#W26d-0:

#extract 2theta column, same for all runs of this sample
W26d_0_C1<-read_xyData(file="Ws26d-o_C01.xrdml", metaData=TRUE)
W26d_0_C1_unlist<-data.frame(t(sapply(W26d_0_C1$dataset,c)))
W26d_0_C1_xy<-W26d_0_C1_unlist$data_block
W26d_0_C1_xy_fin<-as.data.frame(W26d_0_C1_xy)
x2th_W26d_0<-W26d_0_C1_xy_fin$X2Theta

#load .xrdml data with rxylib
data_frame_names_26d_0 <- list.files(pattern = "Ws26d-o_C*") 
data_frame_names_26d_0
data_frame_list_26d_0 <- lapply(data_frame_names_26d_0, read_xyData)  # Read all data frames within working directory
data_frame_list_26d_0   

#Apply formatting functions
res_26d_0 <- lapply(data_frame_list_26d_0, formatfunc)
resxy_26d_0 <-lapply(res_26d_0, formatfunc2)
resxy2_26d_0 <-lapply(resxy_26d_0, formatfunc3)

resxy3_26d_0<-Reduce("+", resxy2_26d_0)

W26d_0_full<-cbind(x2th_W26d_0,resxy3_26d_0)
colnames(W26d_0_full)<-c("x2th", "counts")
W26d_0_full<-as.data.frame(W26d_0_full)

#Remove background
W26d_0_full_bg<-bkg(xrd=W26d_0_full)
W26d_0_full_bg<-cbind(W26d_0_full_bg$tth, W26d_0_full_bg$counts-W26d_0_full_bg$background) #convoluted, want it in a dataframe
W26d_0_full_bg<-as.data.frame(W26d_0_full_bg)
colnames(W26d_0_full_bg)<-c("x2th", "counts")
#Final dataframe with x2theta and counts columns





#W35ac-2:

#extract 2theta column, same for all runs of this sample
W35ac_2_C1<-read_xyData(file="Ws35ac-2_C01.xrdml", metaData=TRUE)
W35ac_2_C1_unlist<-data.frame(t(sapply(W35ac_2_C1$dataset,c)))
W35ac_2_C1_xy<-W35ac_2_C1_unlist$data_block
W35ac_2_C1_xy_fin<-as.data.frame(W35ac_2_C1_xy)
x2th_W35ac_2<-W35ac_2_C1_xy_fin$X2Theta

#load .xrdml data with rxylib
data_frame_names_35ac_2 <- list.files(pattern = "Ws35ac-2_C*") 
data_frame_names_35ac_2
data_frame_list_35ac_2 <- lapply(data_frame_names_35ac_2, read_xyData)  # Read all data frames within working directory
data_frame_list_35ac_2  

#Apply formatting functions
res_35ac_2 <- lapply(data_frame_list_35ac_2, formatfunc)
resxy_35ac_2 <-lapply(res_35ac_2, formatfunc2)
resxy2_35ac_2 <-lapply(resxy_35ac_2, formatfunc3)

resxy3_35ac_2<-Reduce("+", resxy2_35ac_2)

W35ac_2_full<-cbind(x2th_W35ac_2,resxy3_35ac_2)
colnames(W35ac_2_full)<-c("x2th", "counts")
W35ac_2_full<-as.data.frame(W35ac_2_full)

#Remove background
W35ac_2_full_bg<-bkg(xrd=W35ac_2_full)
W35ac_2_full_bg<-cbind(W35ac_2_full_bg$tth, W35ac_2_full_bg$counts-W35ac_2_full_bg$background) #convoluted, want it in a dataframe
W35ac_2_full_bg<-as.data.frame(W35ac_2_full_bg)
colnames(W35ac_2_full_bg)<-c("x2th", "counts")
#Final dataframe with x2theta and counts columns




#--------------------Plotting data - background------------------------------------------
#Compare data with and without background

#W2c_19:
ggplot()+
  geom_line(data=W2c_19_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W2c_19_full, aes(x=x2th, y=counts), col="red")+
  theme_classic()

#W5_3:
ggplot()+
  geom_line(data=W5_3_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W5_3_full, aes(x=x2th, y=counts), col="red")+
  theme_classic()

#W26d_0:
ggplot()+
  geom_line(data=W26d_0_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W26d_0_full, aes(x=x2th, y=counts), col="red")+
  theme_classic()

#W35ac_2:
ggplot()+
  geom_line(data=W35ac_2_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W35ac_2_full, aes(x=x2th, y=counts), col="red")+
  theme_classic()

  
#--------------------Plotting data - all samples-----------------------------------
#Plot all samples together
p_W2c_19<-ggplot()+
  geom_line(data=W2c_19_full_bg, aes(x=x2th, y=counts), col="black")+
  ggtitle("Wc2-19")+
  theme_classic()

p_W5_3<-ggplot()+
  geom_line(data=W5_3_full_bg, aes(x=x2th, y=counts), col="black")+
  ggtitle("W5-3")+
  theme_classic()

p_W26d_0<-ggplot()+
  geom_line(data=W26d_0_full_bg, aes(x=x2th, y=counts), col="black")+
  ggtitle("W26d-0")+
  theme_classic()

p_W35ac_2<-ggplot()+
  geom_line(data=W35ac_2_full_bg, aes(x=x2th, y=counts), col="black")+
  ggtitle("W35ac-2")+
  theme_classic()

cowplot::plot_grid(p_W2c_19, p_W5_3, p_W26d_0, p_W35ac_2, nrow=4)



#Within one plot:
ggplot()+
  geom_line(data=W2c_19_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W5_3_full_bg, aes(x=x2th, y=counts), col="purple")+
  geom_line(data=W26d_0_full_bg, aes(x=x2th, y=counts), col="blue")+
  geom_line(data=W35ac_2_full_bg, aes(x=x2th, y=counts), col="green")+
  geom_segment(aes(x=3, y=10000, xend=6, yend=10000), col="black")+
  geom_segment(aes(x=3, y=9000, xend=6, yend=9000), col="purple")+
  geom_segment(aes(x=3, y=8000, xend=6, yend=8000), col="blue")+
  geom_segment(aes(x=3, y=7000, xend=6, yend=7000), col="green")+
  annotate("text", x = 8, y = 10000, label = "W2c-19")+
  annotate("text", x = 8, y = 9000, label = "W5-3")+
  annotate("text", x = 8, y = 8000, label = "W26d-0")+
  annotate("text", x = 8, y = 7000, label = "W35ac-2")+
  theme_classic()


#Zooming in to relevant interval for clay minerals
ggplot()+
  geom_line(data=W2c_19_full_bg, aes(x=x2th, y=counts), col="black")+
  geom_line(data=W5_3_full_bg, aes(x=x2th, y=counts), col="purple")+
  geom_line(data=W26d_0_full_bg, aes(x=x2th, y=counts), col="blue")+
  geom_line(data=W35ac_2_full_bg, aes(x=x2th, y=counts), col="green")+
  geom_segment(aes(x=1, y=1000, xend=4, yend=1000), col="black")+
  geom_segment(aes(x=1, y=900, xend=4, yend=900), col="purple")+
  geom_segment(aes(x=1, y=800, xend=4, yend=800), col="blue")+
  geom_segment(aes(x=1, y=700, xend=4, yend=700), col="green")+
  annotate("text", x = 5, y = 1000, label = "W2c-19")+
  annotate("text", x = 5, y = 900, label = "W5-3")+
  annotate("text", x = 5, y = 800, label = "W26d-0")+
  annotate("text", x = 5, y = 700, label = "W35ac-2")+
  xlim(0,15)+
  ylim(0,1000)+
  theme_classic()
