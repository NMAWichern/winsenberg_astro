#NMDS-based analysis of the Winsenberg portable XRF dataset
#Measured on a Bruker Titan S1 800 in fall 2022
#40 kV/20 mA/75 s/no filters

#----------------------DATA PREPARATION------------------------------------------
#Load packages
library(tidyverse)     #packages for plotting
library(RColorBrewer)
library(ggpubr)

library(MVN)           #packages for ordination & statistics
library(MASS)
library(vegan)
library(dml)
library(Hmisc)
library(corrplot)

#Set working directory
setwd("C:/Users/nwichern/OneDrive - Universität Münster/Late Devonian black shales/Submitted manuscripts/Winsenberg_Clim_Past/Supplements Winsenberg/Data/compiled for code testing")


#Load XRF data. We carried out analysis on the counts dataset, as it is more complete than the calibrated dataset.
WXRF_dt<-read.csv("Winsenberg_pXRF_counts.csv", header=T)
WXRF_dt[WXRF_dt<0]<-0 #replace any erronous negative values

#Remove elements that are of no interest:
#with stratigraphic height:
WXRFdt<-WXRF_dt[c("Stratigraphic.level..cm.", "Al.K12..counts.", "Ca.K12..counts.", "Cr.K12..counts.", "Cu.K12..counts.", "Fe.K12..counts.", 
                "K.K12..counts.", "Mn.K12..counts.", "P.K12..counts.", "Rb.K12..counts.", "S.K12..counts.", 
                "Si.K12..counts.", "Sr.K12..counts.", "Ti.K12..counts.", "V.K12..counts.","Zr.K12..counts.")] 
colnames(WXRFdt)<-c("dt", "Al.K12", "Ca.K12", "Cr.K12", "Cu.K12", "Fe.K12", "K.K12", "Mn.K12", "P.K12", 
                    "Rb.K12", "S.K12", "Si.K12", "Sr.K12", "Ti.K12", "V.K12","Zr.K12")
#without:
WXRF<-WXRF_dt[c("Al.K12..counts.", "Ca.K12..counts.", "Cr.K12..counts.", "Cu.K12..counts.", "Fe.K12..counts.", 
                  "K.K12..counts.", "Mn.K12..counts.", "P.K12..counts.", "Rb.K12..counts.", "S.K12..counts.", 
                  "Si.K12..counts.", "Sr.K12..counts.", "Ti.K12..counts.", "V.K12..counts.","Zr.K12..counts.")] 
colnames(WXRF)<-c("Al.K12", "Ca.K12", "Cr.K12", "Cu.K12", "Fe.K12", "K.K12", "Mn.K12", "P.K12", 
                    "Rb.K12", "S.K12", "Si.K12", "Sr.K12", "Ti.K12", "V.K12","Zr.K12")


#----------------------DATA ASSESSMENT--------------------------------------------
hist.data.frame(WXRF) #plot histograms of each element
mvn(WXRF) #statistical description of dataset

#The data are not normally distributed. Non-Metric multiDimensional Scaling (NMDS)
#does not require normally distributed data and is therefore suitable for analysis
#(compared to e.g., Principle Component Analysis).


#-------------------------DATA SUBSETS----------------------------------------
#Analyse Kellwasser and non-Kellwasser intervals separately. 
#The following subsets are made:
#LKW = Lower Kellwasser black shale
#iKW = inter-Kellwasser, in between lower and upper Kellwasser black shales
#UKW = Upper Kellwasser black shale
#pKW = post-Kellwasser interval


#Divisions based on log
WXRF_LKW_l<-WXRF[1:67,]  
WXRF_iKW_l<-WXRF[68:396,]  
WXRF_UKW_l<-WXRF[397:443,]
WXRF_pKW_l<-WXRF[444:556,]

#Added column identifying the interval:
LKW_l<-cbind("LKW", WXRF_LKW_l[,1:15])  #construct column
colnames(LKW_l)[1]<-"int"
iKW_l<-cbind("iKW", WXRF_iKW_l[,1:15])
colnames(iKW_l)[1]<-"int"
UKW_l<-cbind("UKW", WXRF_UKW_l[,1:15])
colnames(UKW_l)[1]<-"int"
pKW_l<-cbind("pKW", WXRF_pKW_l[,1:15])
colnames(pKW_l)[1]<-"int"
WXRF_int_l<-rbind(LKW_l, iKW_l, UKW_l, pKW_l)  

#make separate dataframe with interval names and depth
int_l<-cbind(WXRF_int_l$int, WXRFdt$dt) 
int_l<-as.data.frame(int_l)
colnames(int_l)<-c("interval", "dt")
int_l$interval<-as.factor(int_l$interval)


#plot correlations by interval:
dev.off()
pairs(WXRF, cex=0.5, pch=20, 
      col = brewer.pal(4, name="Set1")[int_l$interval],
      lower.panel=NULL)
par(xpd = TRUE)
legend("bottomleft", fill = brewer.pal(4, name="Set1"), 
       legend = c( levels(int_l$interval)))

#To look at any correlation in more detail:
plot(x=WXRF$Mn.K12, y=WXRF$Ca.K12, cex=1.2, pch=20, 
     col = brewer.pal(4, name="Set1")[int_l$interval],
     lower.panel=NULL)
par(xpd = TRUE)
legend("topleft", fill = brewer.pal(4, name="Set1"), 
       legend = c( levels(int_l$interval)))
#Note that these plots are not included in the accompanying paper.



#-----------------------------------NMDS (Fig. C4)-------------------------------------
#NMDS analysis using the vegan package and the code from Bialik et al. 2022, DOI:10.1002/dep2.161
#NMDS can be applied to a wide variety of datasets. 
#The data contains no ordinal or binary values, so "gower" is a good method for distance calculations
#Note that it is an iterative method, so individual outcomes may vary slightly.

set.seed(1)
WXRF_nmds<- metaMDS(WXRF, distance="gower", autotransform = F, wascores = F) #run iterative stresses
#stress has to be >0.2 to be acceptable-->check.

stressplot(WXRF_nmds) #fit is OK

#Format data
data.scores <- as.data.frame(scores(WXRF_nmds))  
data.scores$site <- rownames(data.scores)  
data.scores$int <- int_l$interval
species_object <- envfit(WXRF_nmds, WXRF)
species.scores <- as.data.frame(species_object$vectors$arrows)
species.scores$species <- rownames(species.scores)

#Plot, including variables (=elements)
p_nmds_whole_1v2_litho<-ggplot() + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),size=3) +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=int,colour=int),size=2.5) + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0, alpha=0.5) +
  scale_colour_brewer(palette="Set1") +
  coord_equal()+
  theme_bw()


#Plot, samples only, grouped by interval
grp.LKW <- data.scores[data.scores$int == "LKW", ][chull(data.scores[data.scores$int == 
                                                                       "LKW", c("NMDS1", "NMDS2")]), ]  
grp.iKW <- data.scores[data.scores$int == "iKW", ][chull(data.scores[data.scores$int == 
                                                                       "iKW", c("NMDS1", "NMDS2")]), ]  
grp.UKW <- data.scores[data.scores$int == "UKW", ][chull(data.scores[data.scores$int == 
                                                                       "UKW", c("NMDS1", "NMDS2")]), ]  
grp.pKW <- data.scores[data.scores$int == "pKW", ][chull(data.scores[data.scores$int == 
                                                                       "pKW", c("NMDS1", "NMDS2")]), ] 

hull.data <- rbind(grp.LKW, grp.iKW, grp.UKW, grp.pKW)  

p_nmds_whole_1v2_hulls_litho<-ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=int,group=int),alpha=0.50)+
  scale_fill_manual(values = c("#FFB74D","#67A9CF","#DCEDC8","#E1BEE7"))+
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=int),size=2.5) + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0, alpha=0.5) +
  scale_colour_brewer(palette="Set1") +
  coord_equal()+
  theme_bw()

cowplot::plot_grid(p_nmds_whole_1v2_litho, p_nmds_whole_1v2_hulls_litho, ncol=2)
