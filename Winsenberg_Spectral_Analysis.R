#Spectral analysis carried out on pXRF data from the Winsenberg Road Section, Rhenish Massif, Germany
#Accompanying code to Wichern et al. 
#This script also includes code to plot (components of) additional figures.

#pXRF instrument: Bruker Titan S1 800
#Settings: 40 kV/20 mA/75 s/no filters
#measured in spring-summer 2022

#load required packages
#(if not installed: run install.packages("package name"))
library(tidyverse)
library(RColorBrewer)
library(astrochron)
library(biwavelet)
library(plotly)
library(dml)
library(Hmisc)
library(viridis)
library(palinsol)

setwd("insert own working directory")



#--------------------Insolation curve (part of Fig. 1)--------------------------------------
#plot insolation curve for an interval that contains a 2.4 Myr minimum.
#There is no solution for the Devonian, so we select a random recent interval to illustrate
#the eccentricity minimum hypothesis.

data(LA04)
orbit=c(eps = la04past$eps, ecc = la04past$ecc, varpi = la04past$varpi)

#calculate insolation time series at december solstice, 20S latitude.
insol = c()
for (i in 1:51001){
  insol[i] = Insol(c(eps = la04past$eps[i], ecc = la04past$ecc[i], varpi = la04past$varpi[i]), long = 3*pi/2, lat = -20*pi/180)
}

insol_data<-cbind(la04past$time, la04past$ecc, la04past$eps, la04past$varpi, insol)
insol_data<-as.data.frame(insol_data)

#Plot 2.4 Myr minimum around 4.5 Myrs ago, an the preceding and superceding normal 405 kyr cycles
ggplot(data=insol_data)+
  #geom_line(aes(x=la04past$time, y=la04past$ecc))+
  geom_line(aes(x=la04past$time/1000, y=insol))+
  xlim(-5.2, -4.0)+
  xlab("time [Myr]")+
  ylab("insolation [W/m^2]")+
  coord_flip()+
  theme_classic()


#--------------------------------p-XRF--------------------------------------------
#---------------------------Upload data------------------------------------------

#load calibrated dataset
W_XRF<-read.csv("Winsenberg_pXRF_cal.csv", header=T)
colnames(W_XRF)[3]<-"ID" #rename for ease
colnames(W_XRF)[4]<-"dt" #rename for ease

#cut down to relevant elements:
W_XRF<-W_XRF[c("dt","ID","Al2O3..wt..", "CaO..wt..","K2O..wt..", "SiO2..wt..", "TiO2..wt..")] 
colnames(W_XRF)<-c("dt", "ID", "Al2O3", "CaO", "K2O", "SiO2", "TiO2") #rename for ease



#---------------------Ti vs Al vs carbonate content (Fig. C6)---------------------

#Plot all datapoints with gradual colour scheme, clearly different slopes visible:
a<-ggplot(data=W_XRF)+
  geom_point(aes(x=Al2O3, y=TiO2, col=CaO*1.78))+ #mulitply with 1.78 to obtain approximate CaCO3 content
  scale_colour_continuous(type = "viridis")+
  ylim(-0.1, 2.7)+
  xlim(0,50)+
  labs(title="Ti2O vs Al2O3 based on carbonate content",
       x ="Al2O3 [wt%]", y = "TiO2 [wt%]") +
  theme_classic()

#datapoints plotted with discrete colours to better visualise differences:
d1<-subset(W_XRF, W_XRF$CaO*1.78>65)
d2<-subset(W_XRF, W_XRF$CaO*1.78>=15 & W_XRF$CaO*1.78<=65)
d3<-subset(W_XRF, W_XRF$CaO*1.78<15)
#boundaries: what best fits the observed changes while still keeping an existing classification scheme
#in mind (Correns 1939: <15 is clay-slightly marly clay, 15-65 is marly clay-marl, >65 is calcareous marl - pure limestone)

ggplot()+
  geom_point(data=d1, aes(x=Al2O3, y=TiO2), fill="#fde725", col="#5ec962", shape=23)+
  geom_point(data=d2, aes(x=Al2O3, y=TiO2), fill="#20A387FF", col="#3b528b", shape=21)+
  geom_point(data=d3, aes(x=Al2O3, y=TiO2), fill="#440154FF", col="black", shape=22)+
  labs(title="Ti2O vs Al2O3 based on carbonate content",
       x ="Al2O3 [wt%]", y = "TiO2 [wt%]")+
  theme_classic()


#groups with regression lines plotted to quantify differences in slope:
lm1<-lm(data=d1, TiO2~Al2O3)
lm2<-lm(data=d2, TiO2~Al2O3)
lm3<-lm(data=d3, TiO2~Al2O3)

lm1$coefficients 
lm2$coefficients
lm3$coefficients

  
b<-ggplot()+
  geom_point(data=d1, aes(x=Al2O3, y=TiO2, fill='>65 wt% (Calcareous marl to pure limestone)', col=">65 wt% (Calcareous marl to pure limestone)"), shape=23)+
  geom_point(data=d2, aes(x=Al2O3, y=TiO2, fill='15-65 wt% (Marly clay to marl)', col="15-65 wt% (Marly clay to marl)"), shape=21)+
  geom_point(data=d3, aes(x=Al2O3, y=TiO2, fill='<15 wt% (Clay to slightly marly clay)', col="<15 wt% (Clay to slightly marly clay)"), shape=22)+
  geom_smooth(data=d1,method='lm', formula=y~x, aes(x=Al2O3, y=TiO2), se=F, col="#5ec962")+ #through the origin: do y~x+0 instead of y~x.
  geom_smooth(data=d2,method='lm', formula=y~x,aes(x=Al2O3, y=TiO2), se=F, col="#3b528b")+
  geom_smooth(data=d3,method='lm', formula=y~x,aes(x=Al2O3, y=TiO2), se=F, col="black")+
  annotate("text", x=12, y=0, label="y = 0.016x - 0.007", col="#5ec962", fontface = "bold")+
  annotate("text", x=28, y=0.4, label="y = 0.036x - 0.22", col="#20A387FF", fontface = "bold")+
  annotate("text", x=28, y=1.3, label="y = 0.20x + 0.45", col="#440154FF", fontface = "bold")+
  annotate("text", x=20, y=2.25, label="Classification after Correns (1939)", col="black", fontface = "italic")+
  labs(title="Ti2O vs Al2O3 based on carbonate content - annotated",
     x ="Al2O3 [wt%]", y = "TiO2 [wt%]")+
  scale_fill_manual(name='Carbonate content (CaO*1.78)',
                       breaks=c('>65 wt% (Calcareous marl to pure limestone)', '15-65 wt% (Marly clay to marl)', '<15 wt% (Clay to slightly marly clay)'),
                       values=c('>65 wt% (Calcareous marl to pure limestone)'='#fde725', '15-65 wt% (Marly clay to marl)'='#20A387FF', '<15 wt% (Clay to slightly marly clay)'='#440154FF'))+
  scale_colour_manual(name='Carbonate content (CaO*1.78)',
                      breaks=c('>65 wt% (Calcareous marl to pure limestone)', '15-65 wt% (Marly clay to marl)', '<15 wt% (Clay to slightly marly clay)'),
                      values=c('>65 wt% (Calcareous marl to pure limestone)'='#fde725', '15-65 wt% (Marly clay to marl)'='#20A387FF', '<15 wt% (Clay to slightly marly clay)'='#440154FF'))+
  guides(colour = guide_legend(override.aes = list(shape = c(23,21,22), size=c(4,4,4))))+
  theme_classic()+
  ylim(-0.1, 2.7)+
  xlim(0,50)+
  theme(legend.position=c(0.55, 0.9),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=11))
  
#plot together:
cowplot::plot_grid(a,b, ncol=2)






#--------------------------Create proxy records-----------------------------------
#Create elemental ratio records in a format that can be analysed by astrochron's spectral analysis functions
#All records are linearly interpolated, which is a requirement for MTM spectral analysis
#Done using the linterp() function in astrochron

#Proxy records: SiO2/CaO, TiO2/Al2O3, K2O/Al2O3.

#Named as single elements instead of oxides to shorten object names

#Single elements:
Si<-cbind(W_XRF$dt, W_XRF$SiO2)
colnames(Si)<-c("dt", "Si")
Si<-as.data.frame(Si)

Ca<-cbind(W_XRF$dt, W_XRF$CaO)
colnames(Ca)<-c("dt", "Ca")
Ca<-as.data.frame(Ca)

Al<-cbind(W_XRF$dt, W_XRF$Al2O3)
colnames(Al)<-c("dt", "Al")
Al<-as.data.frame(Al)

Ti<-cbind(W_XRF$dt, W_XRF$TiO2)
Ti<-trim(Ti)
colnames(Ti)<-c("dt", "Ti")
Ti<-as.data.frame(Ti)

K<-cbind(W_XRF$dt, W_XRF$K2O)
colnames(K)<-c("dt", "K")
K<-as.data.frame(K)



#Ratios
#Si/Ca is transformed to a log10 scale in order to deal with intervals with very low Ca,
#which lead to exceptionally high Si/Ca values.
SiCa<-cbind(W_XRF$dt, W_XRF$SiO2/W_XRF$CaO)
colnames(SiCa)<-c("dt", "SiCa")
SiCa<-as.data.frame(SiCa)
SiCa<-cbind(SiCa$dt, log10(SiCa$SiCa))
SiCa<-as.data.frame(SiCa)
SiCa_int<-linterp(SiCa)

#Outlier values are removed from Ti/Al, Zr/Al, and K/Al using the trim() function in astrochron
#These outliers will otherwise obscure the MTM spectra
TiAl<-cbind(W_XRF$dt, W_XRF$TiO2/W_XRF$Al2O3)
colnames(TiAl)<-c("dt", "TiAl")
TiAl<-as.data.frame(TiAl)
TiAl_t<-trim(TiAl)
TiAl_int<-linterp(TiAl_t)

KAl<-cbind(W_XRF$dt, W_XRF$K2O/W_XRF$Al2O3)
colnames(KAl)<-c("dt", "KAl")
KAl<-as.data.frame(KAl)
KAl_t<-trim(KAl)
KAl_int<-linterp(KAl_t)



#-------------------Plot depth domain elemental records (Fig. C3)----------------
#Creat plots for single elements. This is just for inspection, we will analyse the ratios hereafter.
#Includes relative errors based on duplicate measurements of 5 samples (repackaged after each measurement)
#Si: 1.7
#Ca: 0.5
#Al: 1.5
#Ti: 1.5
#K: 0.6


br<-seq(0, 1300, 100)

#Create plots for each elemental ratio record
Si_p<-ggplot()+
  geom_ribbon(data=W_XRF, aes(x=dt, y=SiO2, ymin = SiO2-SiO2*0.017, ymax = SiO2+SiO2*0.017),alpha=0.7, fill="dodgerblue")+
  geom_line(data=Si, aes(x=dt, y=Si), col="black")+
  labs(title="SiO2",
       x ="depth [cm]", y = "content [wt%]")+
  scale_x_continuous(breaks=br, limits=c(0, 1300))+
  coord_flip()+
  theme_classic()

Ca_p<-ggplot()+
  geom_ribbon(data=W_XRF, aes(x=dt, y=CaO, ymin = CaO-CaO*0.005, ymax = CaO+CaO*0.005),alpha=0.7, fill="dodgerblue")+
  geom_line(data=Ca, aes(x=dt, y=Ca), col="black")+
  labs(title="CaO",
       x ="depth [cm]", y = "content [wt%]")+
  scale_x_continuous(breaks=br, limits=c(0, 1300))+
  coord_flip()+
  theme_classic()

Al_p<-ggplot()+
  geom_ribbon(data=W_XRF, aes(x=dt, y=Al2O3, ymin = Al2O3-Al2O3*0.015, ymax = Al2O3+Al2O3*0.015),alpha=0.7, fill="dodgerblue")+
  geom_line(data=Al, aes(x=dt, y=Al), col="black")+
  labs(title="Al2O3",
       x ="depth [cm]", y = "content [wt%]")+
  scale_x_continuous(breaks=br, limits=c(0, 1300))+
  coord_flip()+
  theme_classic()

Ti_p<-ggplot()+
  geom_ribbon(data=W_XRF, aes(x=dt, y=TiO2, ymin = TiO2-TiO2*0.015, ymax = TiO2+TiO2*0.015),alpha=0.7, fill="dodgerblue")+
  geom_line(data=Ti, aes(x=dt, y=Ti), col="black")+
  labs(title="TiO2",
       x ="depth [cm]", y = "content [wt%]")+
  scale_x_continuous(breaks=br, limits=c(0, 1300))+
  coord_flip()+
  theme_classic()


K_p<-ggplot()+
  geom_ribbon(data=W_XRF, aes(x=dt, y=K2O, ymin = K2O-K2O*0.006, ymax = K2O+K2O*0.006),alpha=0.7, fill="dodgerblue")+
  geom_line(data=K, aes(x=dt, y=K), col="black")+
  labs(title="K2O",
       x ="depth [cm]", y = "content [wt%]")+
  scale_x_continuous(breaks=br, limits=c(0, 1300))+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(Si_p, Al_p, Ti_p, K_p, nrow=1)




#---------------------Plot depth domain ratios (Fig. C5)----------------------------
#Create plots for each elemental ratio record
SiCa_p<-ggplot()+
  geom_line(data=SiCa, aes(x=V1, y=V2), col="dodgerblue4")+
 # geom_vline(xintercept=xrd, col="red")+
  labs(title="log10 SiO2/CaO",
       x ="depth [cm]", y = "ratio of wt%'s")+
  xlim(0, 1250)+
  coord_flip()+
  theme_classic()

TiAl_p<-ggplot()+
  geom_line(data=TiAl_t, aes(x=X1, y=X2), col="dodgerblue4")+
 # geom_vline(xintercept=xrd, col="red")+
  labs(title="TiO2/Al2O3",
       x ="depth [cm]", y = "ratio of wt%'s")+
  xlim(0, 1250)+
  coord_flip()+
  theme_classic()

KAl_p<-ggplot()+
  geom_line(data=KAl_t, aes(x=X1, y=X2), col="dodgerblue4")+
 # geom_vline(xintercept=xrd, col="red")+
  labs(title="K2O/Al2O3",
       x ="depth [cm]", y = "ratio of wt%'s")+
  xlim(0, 1250)+
  coord_flip()+
  theme_classic()


#plot these together:
cowplot::plot_grid(SiCa_p, TiAl_p, KAl_p, nrow=1)


#--------------------SiO2/CaO heat map (part of Fig. 5b)--------------------------

#v1=depth v2=value
ggplot(data=SiCa, aes(x=1, y=V1))+
  geom_tile(aes(fill=V2), height=7)+
  scale_fill_viridis_c(option = "cividis", direction = -1)+
  labs(title="log10 SiO2/CaO colour grading",
       y ="depth [cm]", x = "ratio of wt%'s")+
  theme_classic()





#---------------------------SPECTRAL ANALYSIS--------------------------------------

#-------------------initial depth domain assessment (Fig. D1)----------------------
#First analyse all periodograms; we will then focus on the SiO2/CaO record
#This is done in order to be able to test the tuning with the other records

#Periodograms
periodogram(SiCa_int, xmax=0.09)
periodogram(TiAl_int, xmax=0.09)
periodogram(KAl_int, xmax=0.09)

#MTM, significant at 90% CL:
mtm(SiCa_int, tbw=2, detrend=T) 
#Very noisy spectrum. Most identified peaks are high-frequency, and there is no clear Milankovitch signal.




#-----------------------Evolutive assessment (Fig. 5c-d)--------------------------

#As there are lithological changes, changes in sedimentation are to be expected. 
#Therefore, we plot the EHA spectrum to see if we can identify shifts:
eha(SiCa_int, tbw=2, win=100, step=5, detrend=T,genplot=3)

#As well as wavelet spectra:
#Save these to your local drive
wave<-wt(SiCa_int, do.sig=F, dt=5) #same dt as EHA. 

setwd("insert own wd")

pdf("SiCa_int_wavelet.pdf", height=3, width=6)

plot.biwavelet(wave, fill.cols=viridis(n=10, option="D"), plot.cb=F, type="power.corr.norm")
#you can also just run the above line

dev.off()

#There are clear shifts in the dominant spectra throughout the record.
#These changes coincide with lithological changes, suggesting these shifts are related to sedimentation rate changes.





#----------------------------TUNING (Fig. 5c-e)---------------------------------------------
#We tune to the cycle that is easiest to identify and trace within the EHA spectrum
#see paper for more details.

#create dataset with eha data (probability for each frequency)
SiCa_eha<-eha(SiCa_int, tbw=2, win=100, step=5, detrend=T,genplot=3, output=3)

#trace the identified frequency; we chose this one as it is easiest to trace
SiCa_tf<-traceFreq(SiCa_eha, pl=2)

#adjust the record 
SiCa_f2sr<-freq2sedrate(SiCa_tf, 19) #approximate average precession duration of 19 kyr.
#write.csv(SiCa_f2sr, "SiCa_f2sr.csv") <-- Only run this if you want to use your own tracing.

#note that as this is done manually, it will be slightly different each time you try it.
#Our used traced frequency can be loaded from the supplementary data. 


#traced solution used:
SiCa_f2sr<-read.csv("SiCa_f2sr.csv")
SiCa_f2sr[2]<-SiCa_f2sr[2]*100 #because our data is in cm instead of m
colnames(SiCa_f2sr)<-c("X1", "X2") #for compatibility in astrochron

#convert sedimentation rate to timescale
SiCa_sr2t<-sedrate2time(SiCa_f2sr)
#write.csv(SiCa_sr2t, "SiCa_sr2t.csv") see also supplements

#"tune" the record to this timescale
SiCa_tune<-tune(SiCa, controlPts=SiCa_sr2t, extrapolate=T)
#write.csv(SiCa_tune, "SiCa_tune.csv") see also supplements

#interpolate the record
SiCa_tune_int<-linterp(SiCa_tune)
SiCa_tune_int_0<-cbind(SiCa_tune_int$X1+83.7235, SiCa_tune_int$X2) #start record at t=0
#write.csv(SiCa_tune_int, "SiCa_tune_int.csv") see also supplements



#------------------Spectral analysis in the time domain (Fig. 6)---------------------

#check eha to determine whether the sed rate adjustment was carried out correctly
eha(SiCa_tune_int_0, tbw=2, win=150, step=5, detrend=T,genplot=3, fmax=0.15)


#apply spectral analysis
SiCa_tune_mtm<-mtm(SiCa_tune_int, tbw=2, xmax=0.1, sigID=F, output=1)
periodogram(SiCa_tune_int, background=1, xmax=0.1)



#------------------Bandpasses in the time domain (Fig. 7)----------------------------

#bandpass potential Milankovitch cycles
bandpass(SiCa_tune_int, flow=0.08, fhigh=0.10, xmax=0.2) #p
p<-bandpass(SiCa_tune_int, flow=0.045, fhigh=0.065, xmax=0.1) #p
e1<-bandpass(SiCa_tune_int, flow=0.0019, fhigh=0.003, xmax=0.1) #e1
e2<-bandpass(SiCa_tune_int, flow=0.007, fhigh=0.012, xmax=0.1) #e2
bandpass(SiCa_tune_int, flow=0.025, fhigh=0.04, xmax=0.1) #o
bandpass(SiCa_tune_int, flow=0.022, fhigh=0.07, xmax=0.1) #o+p
#there is no clear obliquity signal, but there is a decent amplitude
#modulation fit between p, e1, and e2.

#Plot the amplitude modulation observed from bandpassed frequencies
pam<-hilbert(p) #hilbert transformation
e2am<-hilbert(e2)

ggplot()+
  geom_line(data=p, aes(x=X1, y=X2), col="blue")+
  geom_line(data=pam, aes(x=V1, y=envelope+0.3), col="orange")+
  geom_line(data=e2, aes(x=X1, y=(X2*-1)-0.5), col="blue")+ #plotted inverted to show fit with p
  geom_line(data=e2am, aes(x=V1, y=(envelope)-0.7), col="orange")+
  geom_line(data=e1, aes(x=X1, y=(X2*-1)+1.1), col="blue")+
  labs(title="Bandpass Amplitude Modulation", x="duration [kyr]", y="amplitude")+
  theme_classic()

br<-seq(-100, 1300, 100)

ggplot()+
  geom_line(data=SiCa_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  geom_line(data=p, aes(x=X1+83.7235, y=X2), col="orange")+
  #geom_line(data=pam, aes(x=V1, y=envelope+0.3), col="orange")+
  geom_line(data=e2, aes(x=X1+83.7235, y=(X2*-1)+0.5), col="red")+ #plotted inverted to show fit with p
 # geom_line(data=e2am, aes(x=V1, y=(envelope)-0.7), col="orange")+
  geom_line(data=e1, aes(x=X1+83.7235, y=(X2*-1)+1.1), col="red")+
  labs(title="Bandpass Amplitude Modulation", x="duration [kyr]", y="amplitude")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()



#----------------------------TimeOpt (Fig. 8)------------------------------------

#Test the amplitude modulation fit with TimeOpt:
targete<-c(405, 130.7, 123.8, 94.9) #eccentricity
targetp<-c(21.3, 20.2, 17.37, 17.50) #precession, taken from Waltham (2015)
#https://davidwaltham.com/wp-content/uploads/2014/01/Milankovitch.html
#(last accessed: 06/07/2023)

#convert cm to m, requirement of TimeOpt
SiCa_tune_m<-cbind(SiCa_tune_int$X1/100, SiCa_tune_int$X2) 

#apply timeopt. 
#TimeOpt is meant for depth records, but here it is used solely to test the AM
#Take the time axis as depth axis here, so from ca. 1100 kyr to ca. 1100 cm (11 m)
#To better compare records: look at distance between d13C peaks. Here: 850 m
#Then, the maximum tested duration is ca. 1600 kyr (Ma et al. 2022)--> 850/1600 = 0.53 cm/kyr
#The minimum tested duration is ca. 600 kyr (De Vleeschouwer et al. 2017)-->850/600 = 1.41 cm/kyr
#So, sedmin and sedmax are set to 0.4 and 1.4 cm/kyr, respectively, to cover the range of suspected durations.
SiCa_tune_to_ep<-timeOpt(dat=SiCa_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                      targetP=targetp, fit=1, flow=0.045, fhigh=0.065, output=2)
#Supports visual interpretation shown earlier.

SiCa_tune_to_ee<-timeOpt(dat=SiCa_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                      targetP=targetp, fit=2, flow=0.007, fhigh=0.012, output=2)
#poor fit in the center of the record, but this is expected as this is a suspected eccentricity minimum
  

#Compare 2 timeOpt records and original tuned record:
x<-ggplot()+
  geom_line(data=SiCa_tune_to_ep, aes(x=time, y=value), col="blue")+
  geom_line(data=SiCa_tune_to_ee, aes(x=time, y=value), col="orange")+ 
  geom_line(data=SiCa_tune_int, aes(x=X1, y=X2), col="red")+
  theme_classic()
ggplotly(x)


#formatting
SiCa_to_ep<-SiCa_tune_to_ep[,1:2]
SiCa_to_ee<-SiCa_tune_to_ee[,1:2]

#mtm
mtm(SiCa_to_ep, tbw=2)
mtm(SiCa_to_ee, tbw=2)

#bandpass potential Milankovitch cycles
p<-bandpass(SiCa_to_ep, flow=0.052, fhigh=0.075, xmax=0.2) #p,shifted
e1<-bandpass(SiCa_to_ep, flow=0.0019, fhigh=0.003, xmax=0.02) #e1
e2<-bandpass(SiCa_to_ep, flow=0.007, fhigh=0.013, xmax=0.2) #e2
bandpass(SiCa_to_ep, flow=0.03, fhigh=0.04, xmax=0.2) #o
bandpass(SiCa_to_ep, flow=0.022, fhigh=0.07, xmax=0.2) #o+p
#there is no clear obliquity signal, but there is a decent amplitude
#modulation fit between p, e1, and e2.

#Plot the amplitude modulation observed from bandpassed frequencies
pam<-hilbert(p) #hilbert transformation
e2am<-hilbert(e2)

ggplot()+
  geom_line(data=p, aes(x=time, y=value), col="blue")+
  geom_line(data=pam, aes(x=V1, y=envelope+0.3), col="orange")+
  geom_line(data=e2, aes(x=time, y=(value*-1)-0.5), col="blue")+ #plotted inverted to show fit with p
  geom_line(data=e2am, aes(x=V1, y=(envelope)-0.7), col="orange")+
  geom_line(data=e1, aes(x=time, y=(value*-1)+1.1), col="blue")+
  labs(title="Bandpass Amplitude Modulation", x="duration [kyr]", y="amplitude")+
  theme_classic()




#bandpass potential Milankovitch cycles
p<-bandpass(SiCa_to_ee, flow=0.045, fhigh=0.065, xmax=0.2) #p, not shifted
e1<-bandpass(SiCa_to_ee, flow=0.0019, fhigh=0.003, xmax=0.02) #e1
e2<-bandpass(SiCa_to_ee, flow=0.007, fhigh=0.013, xmax=0.2) #e2
bandpass(SiCa_to_ee, flow=0.03, fhigh=0.04, xmax=0.2) #o
bandpass(SiCa_to_ee, flow=0.027, fhigh=0.065, xmax=0.2) #o+p
#there is no clear obliquity signal, but there is a decent amplitude
#modulation fit between p, e1, and e2.

#Plot the amplitude modulation observed from bandpassed frequencies
pam<-hilbert(p) #hilbert transformation
e2am<-hilbert(e2)

ggplot()+
  geom_line(data=p, aes(x=time, y=value), col="blue")+
  geom_line(data=pam, aes(x=V1, y=envelope+0.3), col="orange")+
  geom_line(data=e2, aes(x=time, y=(value*-1)-0.5), col="blue")+ #plotted inverted to show fit with p
  geom_line(data=e2am, aes(x=V1, y=(envelope)-0.7), col="orange")+
  geom_line(data=e1, aes(x=time, y=(value*-1)+1.1), col="blue")+
  labs(title="Bandpass Amplitude Modulation", x="duration [kyr]", y="amplitude")+
  theme_classic()




#----------------------Tune other proxy records-------------------------------------

#Transfer the Ti/Al and K/Al records to the constructed floating timescale
#and check with spectral analysis
TiAl_tune<-tune(TiAl_t, controlPts=SiCa_sr2t, extrapolate=T) 
KAl_tune<-tune(KAl_t, controlPts=SiCa_sr2t, extrapolate=T)

write.csv(TiAl_tune, "TiAl_tune.csv") #see supplements; note that these have been formatted 
write.csv(KAl_tune, "KAl_tune.csv")   #in order to comply with PANGAEA standards
write.csv(SiCa_tune, "SiCa_tune.csv")

TiAl_tune_int<-linterp(TiAl_tune)
KAl_tune_int<-linterp(KAl_tune)

write.csv(TiAl_tune_int, "TiAl_tune_int.csv")
write.csv(KAl_tune_int, "KAl_tune_int.csv")

mtm(TiAl_tune_int, tbw=2)
mtm(KAl_tune_int, tbw=2)

periodogram(KAl_tune_int)
periodogram(TiAl_tune_int)

eha(TiAl_tune_int, tbw=2, win=100, step=5, detrend=T,genplot=3)
eha(KAl_tune_int, tbw=2, win=100, step=5, detrend=T,genplot=3, fmax=0.15)




#-----------------Proxies in time domain (Fig. 9)----------------------------
br<-seq(-100, 1300, 100)

#Create plots for each elemental ratio record
SiCa_tune_p<-ggplot()+
  geom_line(data=SiCa_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+ #shift base to 0
 labs(title="log10 SiO2/CaO",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic() 

TiAl_tune_p<-ggplot()+
  geom_line(data=TiAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  labs(title="TiO2/Al2O3",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

KAl_tune_p<-ggplot()+
  geom_line(data=KAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  labs(title="K2O/Al2O3",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()


#plot these together:
cowplot::plot_grid(SiCa_tune_p, TiAl_tune_p, KAl_tune_p, nrow=1)






br<-seq(-100, 1300, 100)

TiAl_tune_lp<-lowpass(TiAl_tune_int, fcut=0.003)
KAl_tune_lp<-lowpass(KAl_tune_int, fcut=0.003)
SiCa_tune_lp<-lowpass(SiCa_tune_int, fcut=0.003)

#Now with high/low value vs mean and lowpass filters
SiCa_tune_p<-ggplot()+
  geom_line(data=SiCa_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+ #shift base to 0
 # stat_difference(data=SiCa_tune, aes(x=X1+83.7235, ymin = mean(X2), ymax = X2), alpha = 0.3)+
  geom_line(data=SiCa_tune_lp, aes(x=X1+83.7235, y=X2), col="black")+
  labs(title="log10 SiO2/CaO",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic() 

TiAl_tune_p<-ggplot()+
  geom_line(data=TiAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
 # stat_difference(data=TiAl_tune, aes(x=X1+83.7235, ymin = mean(X2), ymax = X2), alpha = 0.3)+
  geom_line(data=TiAl_tune_lp, aes(x=X1+83.7235, y=X2), col="black")+
  labs(title="TiO2/Al2O3",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

KAl_tune_p<-ggplot()+
  geom_line(data=KAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  #stat_difference(data=KAl_tune, aes(x=X1+83.7235, ymin = mean(X2), ymax = X2), alpha = 0.3)+
  geom_line(data=KAl_tune_lp, aes(x=X1+83.7235, y=X2), col="black")+
  labs(title="K2O/Al2O3",
       x ="duration [kyr]", y = "ratio of wt%'s")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  scale_y_reverse()+ #to improve readability of interpretation as low K/Al = high chemical weathering
  theme_classic()


#plot these together:
cowplot::plot_grid(SiCa_tune_p, TiAl_tune_p, KAl_tune_p, nrow=1)







#Test the amplitude modulation fit with TimeOpt (not in paper):
targete<-c(405, 130.7, 123.8, 94.9) #eccentricity
targetp<-c(21.3, 20.2, 17.37, 17.50) #precession, taken from Waltham (2015)

#convert cm to m
TiAl_tune_m<-cbind(TiAl_tune_int$X1/100, TiAl_tune_int$X2)
KAl_tune_m<-cbind(KAl_tune_int$X1/100, KAl_tune_int$X2) 

#apply timeopt. 
#TimeOpt is meant for depth records, but here it is used solely to test the AM
#Take the time axis as depth axis here, so from ca. 1100 kyr to ca. 1100 cm (11 m)
#To better compare records: look at distance between d13C peaks. Here: 850 m
#Then, the maximum tested duration is ca. 1600 kyr (Ma et al. 2022)--> 850/1600 = 0.53 cm/kyr
#The minimum tested duration is ca. 600 kyr (De Vleeschouwer et al. 2017)-->850/600 = 1.41 cm/kyr
#So, sedmin and sedmax are set to 0.4 and 1.4 cm/kyr, respectively, to cover the range of suspected durations.
TiAl_tune_to_ep<-timeOpt(dat=TiAl_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                         targetP=targetp, fit=1, flow=0.045, fhigh=0.065, output=2)
#Supports visual interpretation shown earlier.

TiAl_tune_to_ee<-timeOpt(dat=TiAl_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                         targetP=targetp, fit=2, flow=0.007, fhigh=0.012, output=2)
#poor fit in the center of the record, but this is expected as this is a suspected eccentricity minimum


KAl_tune_to_ep<-timeOpt(dat=KAl_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                         targetP=targetp, fit=1, flow=0.045, fhigh=0.065, output=2)
#Supports visual interpretation shown earlier.

KAl_tune_to_ee<-timeOpt(dat=KAl_tune_m, sedmin=0.4, sedmax=1.7, numsed=200, targetE=targete, 
                         targetP=targetp, fit=2, flow=0.007, fhigh=0.012, output=2)
#poor fit in the center of the record, but this is expected as this is a suspected eccentricity minimum






#-----------------------Proxies bandpassed (Fig. D3)-----------------------------

#With bandpassed filters:
pSiCa<-bandpass(SiCa_tune_int, flow=0.045, fhigh=0.065, xmax=0.1) #p
e1SiCa<-bandpass(SiCa_tune_int, flow=0.0019, fhigh=0.003, xmax=0.1) #e1
e2SiCa<-bandpass(SiCa_tune_int, flow=0.007, fhigh=0.012, xmax=0.1) #e2
bandpass(SiCa_tune_int, flow=0.027, fhigh=0.035, xmax=0.1) #o
bandpass(SiCa_tune_int, flow=0.022, fhigh=0.07, xmax=0.1) #o+p

pTiAl<-bandpass(TiAl_tune_int, flow=0.045, fhigh=0.065, xmax=0.1) #p
e1TiAl<-bandpass(TiAl_tune_int, flow=0.0019, fhigh=0.003, xmax=0.1) #e1
e2TiAl<-bandpass(TiAl_tune_int, flow=0.007, fhigh=0.012, xmax=0.1) #e2
bandpass(TiAl_tune_int, flow=0.027, fhigh=0.035, xmax=0.1) #o
bandpass(TiAl_tune_int, flow=0.022, fhigh=0.07, xmax=0.1) #o+p

pKAl<-bandpass(KAl_tune_int, flow=0.045, fhigh=0.065, xmax=0.1) #p
e1KAl<-bandpass(KAl_tune_int, flow=0.0019, fhigh=0.003, xmax=0.1) #e1
e2KAl<-bandpass(KAl_tune_int, flow=0.007, fhigh=0.012, xmax=0.1) #e2
bandpass(KAl_tune_int, flow=0.059, fhigh=0.072, xmax=0.1) #P
bandpass(KAl_tune_int, flow=0.023, fhigh=0.035, xmax=0.1) #o
bandpass(KAl_tune_int, flow=0.022, fhigh=0.07, xmax=0.1) #o+p

br<-seq(-100, 1300, 100)

SiCa_p_tune<-ggplot()+
  geom_line(data=SiCa_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  geom_line(data=pSiCa, aes(x=X1+83.7235, y=X2), col="orange")+
  geom_line(data=e2SiCa, aes(x=X1+83.7235, y=(X2)), col="red")+
  geom_line(data=e1SiCa, aes(x=X1+83.7235, y=(X2)+0.6), col="red")+
  labs(title="log10 SiO2/CaO", x="duration [kyr]", y="amplitude")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

TiAl_p_tune<-ggplot()+
  geom_line(data=TiAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  geom_line(data=pTiAl, aes(x=X1+83.7235, y=X2), col="orange")+
  geom_line(data=e2TiAl, aes(x=X1+83.7235, y=(X2)), col="red")+ 
  geom_line(data=e1TiAl, aes(x=X1+83.7235, y=(X2)+0.006), col="red")+
  labs(title="TiO2/Al2O3", x="duration [kyr]", y="amplitude")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

KAl_p_tune<-ggplot()+
  geom_line(data=KAl_tune, aes(x=X1+83.7235, y=X2), col="dodgerblue4")+
  geom_line(data=pKAl, aes(x=X1+83.7235, y=X2), col="orange")+
  geom_line(data=e2KAl, aes(x=X1+83.7235, y=(X2)), col="red")+ 
  geom_line(data=e1KAl, aes(x=X1+83.7235, y=(X2)+0.015), col="red")+
  labs(title="K2O/Al2O3", x="duration [kyr]", y="amplitude")+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()



#plot these together:
cowplot::plot_grid(SiCa_p_tune, TiAl_p_tune, KAl_p_tune, nrow=1)



#---------------------Tune carbon isotopes (Fig. 9)-----------------------------------

#Load Carbon isotope records
W_d13Corg<-read.csv("d13Corg_winsenberg.csv", header=T, fileEncoding='latin1')
W_d13Corg<-W_d13Corg[,3:4]
colnames(W_d13Corg)<-c("dt", "d13Corg")

W_d13Ccarb<-read.csv("d13Ccarb_Winsenberg.csv", header=T, fileEncoding='latin1')
W_d13Ccarb<-W_d13Ccarb[,c(3,6)]
colnames(W_d13Ccarb)<-c("dt", "d13Ccarb")
W_d13Ccarb<-W_d13Ccarb[2:144,] #first value is out of tuning bounds, remove


#tune using SiO2/CaO record
d13Corg_tune<-tune(W_d13Corg, controlPts=SiCa_sr2t, extrapolate=T) 
d13Ccarb_tune<-tune(W_d13Ccarb, controlPts=SiCa_sr2t, extrapolate=T) 

#Plot:
d13Corg_tune_p<-ggplot()+
  geom_line(data=d13Corg_tune, aes(x=X1+83.7235, y=X2))+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

d13Ccarb_tune_p<-ggplot()+
  geom_line(data=d13Ccarb_tune, aes(x=X1+83.7235, y=X2))+
  scale_x_continuous(breaks=br, limits=c(-10, 1300))+
  coord_flip()+
  theme_classic()

cowplot::plot_grid(d13Ccarb_tune_p, d13Corg_tune_p)
