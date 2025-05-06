#clear that brain

rm(list = ls())

## set directory

setwd("")

#install.packages('glmmTMB')

#Statistical analyses for PAM at T1


library(glmmTMB)
library(lmerTest)
library(emmeans)
library(sjPlot)
library(drc)
library(ggplot2)
library(Rmisc)
library(readxl)

library(dplyr)
library(car)
library(reshape2)
library(grid)
library(gridExtra)
library(corrgram)
library (tidyr)
library (ggpubr) 
library (tidyverse)
library(lubridate)
library(scales)
library (writexl)
library(extrafont)
library(patchwork)
library(grid)


#Load file

AF <- read_excel("AF_DRC.xlsx")

AF$Geno<- as.factor(AF$Geno)
str(AF)

#Subset data
Control_check1<-subset(AF,  Temp=="31")
Control_check<-subset(Control_check1,Timepoint=="Hold")

AF_no_field<-subset(AF, Temp=="31" | Temp=="35" | Temp=="37" | Temp=="40")
AF_no_field$Site_Source<- as.factor(AF_no_field$Site_Source)
AF_no_field$Source<- as.factor(AF_no_field$Source)

Hold_all<-subset(AF_no_field, Timepoint=="Hold")


###### Compute and plot T1 Fv/Fm ED50 for acute (CBASS) 

Hold_all$Temp<-as.character(Hold_all$Temp)
Hold_all$Temp<-as.numeric(Hold_all$Temp)
str(Hold_all)

#### Compare curves and ED50s ####

#compromise model

AF_DRC <- drm(PAM ~ Temp, data = Hold_all, curveid = Site_Source, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl= c(100, 0.8, NA),lowerl = c(10, 0.3, 30))
mselect(AF_DRC, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(AF_DRC)
summary(AF_DRC)

#extract ED50

WWN_coeff<-AF_DRC$coefficients[9]
WWD_coeff<-AF_DRC$coefficients[10]
LWN_coeff<-AF_DRC$coefficients[11]
LWD_coeff<-AF_DRC$coefficients[12]

#### Run individually for plotting ####

### Windward nursery Acropora florida

#Run model
WWN<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Nursery",],fct = LL.3())
summary(WWN)
plot(WWN)

### WW donor AF

WWD<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Donor",],fct = LL.3())
summary(WWD)
plot(WWD)

###leeward nursery AF

LWN<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Nursery",],fct = LL.3())
summary(LWN)
plot(LWN)


### leeward donor

#Run model
LWD<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Donor",],fct = LL.3())
summary(LWD)
plot(LWD)

#### Combine ED50 data plus predict curves from models for plotting ####

AF_coeffs<-data.frame(WWN_coeff, WWD_coeff, LWN_coeff, LWD_coeff)

WWN_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWN_pred = as.data.frame(predict(WWN, newdata = WWN_preddata, interval = 'confidence'))
WWN_preddata = data.frame(WWN_preddata, fvfm = WWN_pred$Prediction, Lower = WWN_pred$Lower, Upper = WWN_pred$Upper)


WWD_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWD_pred = as.data.frame(predict(WWD, newdata = WWD_preddata, interval = 'confidence'))
WWD_preddata = data.frame(WWD_preddata, fvfm = WWD_pred$Prediction, Lower = WWD_pred$Lower, Upper = WWD_pred$Upper)


LWN_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWN_pred = as.data.frame(predict(LWN, newdata = LWN_preddata, interval = 'confidence'))
LWN_preddata = data.frame(LWN_preddata, fvfm = LWN_pred$Prediction, Lower = LWN_pred$Lower, Upper = LWN_pred$Upper)

LWD_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWD_pred = as.data.frame(predict(LWD, newdata = LWD_preddata, interval = 'confidence'))
LWD_preddata = data.frame(LWD_preddata, fvfm = LWD_pred$Prediction, Lower = LWD_pred$Lower, Upper = LWD_pred$Upper)

#### PLOT  ####
levels(Hold_all$Site)

AF_DRC_plot<- ggplot() +
  geom_jitter(data = Hold_all, aes(x = Temp, y = PAM, color = Site_Source), size = 1.5, width = 0.4) +
  scale_x_continuous(limits=c(30,42), breaks=c(31,34,37,40)) +
  scale_y_continuous(limits=c(0.0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = WWN_preddata, aes(x = temp, y = fvfm), color = '#66ccfe', show.legend = FALSE) +
  geom_ribbon(data = WWN_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#66ccfe', linetype=2, alpha = 0.2) +
  geom_vline(data = AF_coeffs, aes(xintercept = WWN_coeff), color = '#66ccfe', show.legend = FALSE) +
  geom_text(data = AF_coeffs, aes(label=round(WWN_coeff, digits=2), size = 30.0), x = 41.0, y = 0.66, show.legend = FALSE, color = '#66ccfe') +
  
  geom_line(data = WWD_preddata, aes(x = temp, y = fvfm), color = 'darkblue', show.legend = FALSE) +
  geom_ribbon(data = WWD_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'darkblue', linetype=2, alpha = 0.2) +
  geom_vline(data = AF_coeffs, aes(xintercept = WWD_coeff), color = 'darkblue', show.legend = FALSE) +
  geom_text(data = AF_coeffs, aes(label=round(WWD_coeff, digits=2), size = 30.0), x = 41.0, y = 0.69, show.legend = FALSE, color = 'darkblue') +
  
  
  geom_line(data = LWN_preddata, aes(x = temp, y = fvfm), color = 'gold', show.legend = FALSE) +
  geom_ribbon(data = LWN_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'gold', linetype=2, alpha = 0.2) +
  geom_vline(data = AF_coeffs, aes(xintercept = LWN_coeff), color = 'gold', show.legend = FALSE) +
  geom_text(data = AF_coeffs, aes(label=round(LWN_coeff, digits=2), size = 30.0), x = 41.0, y = 0.72, show.legend = FALSE, color = 'gold') +  
  
  geom_line(data = LWD_preddata, aes(x = temp, y = fvfm), color = '#e92000', show.legend = FALSE) +
  geom_ribbon(data = LWD_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#e92000', linetype=2, alpha = 0.2) +
  geom_vline(data = AF_coeffs, aes(xintercept = LWD_coeff), color = '#e92000', show.legend = FALSE) +
  geom_text(data = AF_coeffs, aes(label=round(LWD_coeff, digits=2), size = 30.0), x = 41.0, y = 0.75, show.legend = FALSE, color = '#e92000') +
  
  scale_color_manual(values=c('#e92000', 'gold', 'darkblue', '#66ccfe')) +
  ylab("Fv/Fm") +
  xlab("") +
  theme_light()+
  theme(legend.position = c(0.25, 0.20)) +
  theme(text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12))
theme(legend.background = element_rect(fill = 'white'))



AF_DRC_plot


#Load file

DN2<- read_excel("AC_Donor-Nursery_R.xlsx")

DN2$Geno<- as.factor(DN2$Geno)
str(DN2)

#Subset data
Control_check1<-subset(DN2,  Temp=="31")
Control_check<-subset(Control_check1,Timepoint=="Hold")

DN2_no_field<-subset(DN2, Temp=="31" | Temp=="35" | Temp=="37" | Temp=="40")
DN2_no_field$Site_Source<- as.factor(DN2_no_field$Site_Source)
DN2_no_field$Source<- as.factor(DN2_no_field$Source)

Hold_all<-subset(DN2_no_field, Timepoint=="Hold")


###### Compute and plot T1 Fv/Fm ED50 for acute (CBASS) 

Hold_all$Temp<-as.character(Hold_all$Temp)
Hold_all$Temp<-as.numeric(Hold_all$Temp)
str(Hold_all)

#### Compare curves and ED50s #### Acropora florida first ####
#compromise model

AC_DRC <- drm(PAM ~ Temp, data = Hold_all, curveid = Site_Source, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl= c(100, 0.8, NA),lowerl = c(10, 0.3, 30))
mselect(AC_DRC, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(AC_DRC)
summary(AC_DRC)

WWNAC1_coeff<-AC_DRC$coefficients[9]
WWDAC1_coeff<-AC_DRC$coefficients[10]
LWNAC1_coeff<-AC_DRC$coefficients[11]
LWDAC1_coeff<-AC_DRC$coefficients[12]

### Acropora cytherea

#Run model
WWNAC1<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Nursery",],fct = LL.3())
summary(WWNAC1)
plot(WWNAC1)

### WW donor AF

WWDAC1<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Donor",],fct = LL.3())
summary(WWDAC1)
plot(WWDAC1)

###leeward nursery AF

LWNAC1<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Nursery",],fct = LL.3())
summary(LWNAC1)
plot(LWNAC1)


### leeward donor

#Run model
LWDAC1<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Donor",],fct = LL.3())
summary(LWDAC1)
plot(LWDAC1)

####Acropora cytherea


AC_coeffs<-data.frame(WWNAC1_coeff, WWDAC1_coeff, LWNAC1_coeff, LWDAC1_coeff)

WWNAC1_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWNAC1_pred = as.data.frame(predict(WWNAC1, newdata = WWNAC1_preddata, interval = 'confidence'))
WWNAC1_preddata = data.frame(WWNAC1_preddata, fvfm = WWNAC1_pred$Prediction, Lower = WWNAC1_pred$Lower, Upper = WWNAC1_pred$Upper)


WWDAC1_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWDAC1_pred = as.data.frame(predict(WWDAC1, newdata = WWDAC1_preddata, interval = 'confidence'))
WWDAC1_preddata = data.frame(WWDAC1_preddata, fvfm = WWDAC1_pred$Prediction, Lower = WWDAC1_pred$Lower, Upper = WWDAC1_pred$Upper)


LWNAC1_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWNAC1_pred = as.data.frame(predict(LWNAC1, newdata = LWNAC1_preddata, interval = 'confidence'))
LWNAC1_preddata = data.frame(LWNAC1_preddata, fvfm = LWNAC1_pred$Prediction, Lower = LWNAC1_pred$Lower, Upper = LWNAC1_pred$Upper)

LWDAC1_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWDAC1_pred = as.data.frame(predict(LWDAC1, newdata = LWDAC1_preddata, interval = 'confidence'))
LWDAC1_preddata = data.frame(LWDAC1_preddata, fvfm = LWDAC1_pred$Prediction, Lower = LWDAC1_pred$Lower, Upper = LWDAC1_pred$Upper)



#### PLOT  ####
levels(Hold_all$Site)

AC_DRC_plot<- ggplot() +
  geom_jitter(data = Hold_all, aes(x = Temp, y = PAM, color = Site_Source), size = 1.5, width = 0.4) +
  scale_x_continuous(limits=c(30,42), breaks=c(31,34,37,40)) +
  scale_y_continuous(limits=c(0.0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = LWDAC1_preddata, aes(x = temp, y = fvfm), color = 'red', show.legend = FALSE) +
  geom_ribbon(data = LWDAC1_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red', linetype=2, alpha = 0.2) +
  geom_vline(data = AC_coeffs, aes(xintercept = LWDAC1_coeff), color = 'red', show.legend = FALSE) +
  geom_text(data = AC_coeffs, aes(label = round(LWDAC1_coeff, digits=2), size = 30.0), x = 41.0, y = 0.75, show.legend = FALSE, color = 'red') +
  
  geom_line(data = LWNAC1_preddata, aes(x = temp, y = fvfm), color = 'gold', show.legend = FALSE) +
  geom_ribbon(data = LWNAC1_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'gold', linetype=2, alpha = 0.2) +
  geom_vline(data = AC_coeffs, aes(xintercept = LWNAC1_coeff), color = 'gold' , show.legend = FALSE) +
  geom_text(data = AC_coeffs, aes(label = round(LWNAC1_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.72, show.legend = FALSE, color = 'gold') +  
  
  geom_line(data = WWDAC1_preddata, aes(x = temp, y = fvfm), color = '#66ccfe', show.legend = FALSE) +
  geom_ribbon(data = WWDAC1_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#66ccfe', linetype=2, alpha = 0.2) +
  geom_vline(data = AC_coeffs, aes(xintercept = WWDAC1_coeff), color = '#66ccfe', show.legend = FALSE) +
  geom_text(data = AC_coeffs, aes(label = round(WWDAC1_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.66, show.legend = FALSE, color = '#66ccfe') +
  
  geom_line(data = WWNAC1_preddata, aes(x = temp, y = fvfm), color = 'blue', show.legend = FALSE) +
  geom_ribbon(data = WWNAC1_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'blue', linetype=2, alpha = 0.2) +
  geom_vline(data = AC_coeffs, aes(xintercept = WWNAC1_coeff), color = 'blue', show.legend = FALSE) +
  geom_text(data = AC_coeffs, aes(label = round(WWNAC1_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.69, show.legend = FALSE, color = 'blue') +
  
  scale_color_manual(values=c('#e92000', 'gold', 'darkblue', '#66ccfe')) +
  ylab("Fv/Fm") +
  xlab("Temperature°C") +
  theme_light()+
  theme(legend.position = 'none') +
  theme(text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))


AC_DRC_plot

#Load file number 3

DN3<- read_excel("AC_Donor-Nursery_R2.xlsx")

DN3$Geno<- as.factor(DN3$Geno)
str(DN3)

#Subset data
Control_check1<-subset(DN3,  Temp=="31")
Control_check<-subset(Control_check1,Timepoint=="Hold")

DN3_no_field<-subset(DN3, Temp=="31" | Temp=="35" | Temp=="37" | Temp=="40")
DN3_no_field$Site_Source<- as.factor(DN3_no_field$Site_Source)
DN3_no_field$Source<- as.factor(DN3_no_field$Source)

Hold_all<-subset(DN3_no_field, Timepoint=="Hold")


###### Compute and plot T1 Fv/Fm ED50 for acute (CBASS) 

Hold_all$Temp<-as.character(Hold_all$Temp)
Hold_all$Temp<-as.numeric(Hold_all$Temp)
str(Hold_all)

#### Compare curves and ED50s #### Acropora florida first ####
#compromise model

AC_DRC2 <- drm(PAM ~ Temp, data = Hold_all, curveid = Site_Source, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl= c(100, 0.8, NA),lowerl = c(10, 0.3, 30))
mselect(AC_DRC2, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(AC_DRC2)
summary(AC_DRC2)

WWNAC2_coeff<-AC_DRC2$coefficients[9]
WWDAC2_coeff<-AC_DRC2$coefficients[10]
LWNAC2_coeff<-AC_DRC2$coefficients[11]
LWDAC2_coeff<-AC_DRC2$coefficients[12]

### Acropora cytherea

#Run model
WWNAC2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Nursery",],fct = LL.3())


### WW donor AF

WWDAC2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Donor",],fct = LL.3())


###leeward nursery AF

LWNAC2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Nursery",],fct = LL.3())



### leeward donor

#Run model
LWDAC2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Donor",],fct = LL.3())
summary(LWDAC2)
plot(LWDAC2)

####Acropora cytherea


AC2_coeffs<-data.frame(WWNAC2_coeff, WWDAC2_coeff, LWNAC2_coeff, LWDAC2_coeff)

WWNAC2_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWNAC2_pred = as.data.frame(predict(WWNAC2, newdata = WWNAC2_preddata, interval = 'confidence'))
WWNAC2_preddata = data.frame(WWNAC2_preddata, fvfm = WWNAC2_pred$Prediction, Lower = WWNAC2_pred$Lower, Upper = WWNAC2_pred$Upper)


WWDAC2_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWDAC2_pred = as.data.frame(predict(WWDAC2, newdata = WWDAC2_preddata, interval = 'confidence'))
WWDAC2_preddata = data.frame(WWDAC2_preddata, fvfm = WWDAC2_pred$Prediction, Lower = WWDAC2_pred$Lower, Upper = WWDAC2_pred$Upper)


LWNAC2_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWNAC2_pred = as.data.frame(predict(LWNAC2, newdata = LWNAC2_preddata, interval = 'confidence'))
LWNAC2_preddata = data.frame(LWNAC2_preddata, fvfm = LWNAC2_pred$Prediction, Lower = LWNAC2_pred$Lower, Upper = LWNAC2_pred$Upper)

LWDAC2_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWDAC2_pred = as.data.frame(predict(LWDAC2, newdata = LWDAC2_preddata, interval = 'confidence'))
LWDAC2_preddata = data.frame(LWDAC2_preddata, fvfm = LWDAC2_pred$Prediction, Lower = LWDAC2_pred$Lower, Upper = LWDAC2_pred$Upper)



#### PLOT  ####
levels(Hold_all$Site)


AC_DRC_plot2<- ggplot() +
  geom_jitter(data = Hold_all, aes(x = Temp, y = PAM, color = Site_Source), size = 1.5, width = 0.4) +
  scale_x_continuous(limits=c(30,42), breaks=c(31,34,37,40)) +
  scale_y_continuous(limits=c(0.0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = LWDAC2_preddata, aes(x = temp, y = fvfm), color = 'red', show.legend = FALSE) +
  geom_ribbon(data = LWDAC2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red', linetype=2, alpha = 0.2) +
  geom_vline(data = AC2_coeffs, aes(xintercept = LWDAC2_coeff), color = 'red', show.legend = FALSE) +
  geom_text(data = AC2_coeffs, aes(label = round(LWDAC2_coeff, digits=2), size = 30.0), x = 41.0, y = 0.75, show.legend = FALSE, color = 'red') +
  
  geom_line(data = LWNAC2_preddata, aes(x = temp, y = fvfm), color = 'gold', show.legend = FALSE) +
  geom_ribbon(data = LWNAC2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'gold', linetype=2, alpha = 0.2) +
  geom_vline(data = AC2_coeffs, aes(xintercept = LWNAC2_coeff), color = 'gold' , show.legend = FALSE) +
  geom_text(data = AC2_coeffs, aes(label = round(LWNAC2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.72, show.legend = FALSE, color = 'gold') +  
  
  geom_line(data = WWDAC2_preddata, aes(x = temp, y = fvfm), color = '#66ccfe', show.legend = FALSE) +
  geom_ribbon(data = WWDAC2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#66ccfe', linetype=2, alpha = 0.2) +
  geom_vline(data = AC2_coeffs, aes(xintercept = WWDAC2_coeff), color = '#66ccfe', show.legend = FALSE) +
  geom_text(data = AC2_coeffs, aes(label = round(WWDAC2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.66, show.legend = FALSE, color = '#66ccfe') +
  
  geom_line(data = WWNAC2_preddata, aes(x = temp, y = fvfm), color = 'blue', show.legend = FALSE) +
  geom_ribbon(data = WWNAC2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'blue', linetype=2, alpha = 0.2) +
  geom_vline(data = AC2_coeffs, aes(xintercept = WWNAC2_coeff), color = 'blue', show.legend = FALSE) +
  geom_text(data = AC2_coeffs, aes(label = round(WWNAC2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.69, show.legend = FALSE, color = 'blue') +
  
  scale_color_manual(values=c('#e92000', 'gold', 'darkblue', '#66ccfe')) +
  ylab("") +
  xlab("Temperature°C") +
  theme_light()+
  theme(legend.position = 'none') +
  theme(text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))


AC_DRC_plot2


#Load file number 4

DN4<- read_excel("AF_DRC2.xlsx")

DN4$Geno<- as.factor(DN4$Geno)
str(DN4)

#Subset data
Control_check1<-subset(DN4,  Temp=="31")
Control_check<-subset(Control_check1,Timepoint=="Hold")

DN4_no_field<-subset(DN4, Temp=="31" | Temp=="35" | Temp=="37" | Temp=="40")
DN4_no_field$Site_Source<- as.factor(DN4_no_field$Site_Source)
DN4_no_field$Source<- as.factor(DN4_no_field$Source)

Hold_all<-subset(DN4_no_field, Timepoint=="Hold")


###### Compute and plot T1 Fv/Fm ED50 for acute (CBASS) 

Hold_all$Temp<-as.character(Hold_all$Temp)
Hold_all$Temp<-as.numeric(Hold_all$Temp)
str(Hold_all)

#### Compare curves and ED50s #### Acropora florida first ####
#compromise model

AF_DRC2 <- drm(PAM ~ Temp, data = Hold_all, curveid = Site_Source, fct = LL.3(names = c('hill', 'max', 'ed50')), upperl= c(100, 0.8, NA),lowerl = c(10, 0.3, 30))
mselect(AF_DRC2, list(LL.2(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.4(), LL2.5(), AR.2(), AR.3(), EXD.2(), EXD.3()), icfct = AIC)
modelFit(AF_DRC2)
summary(AF_DRC2)

WWNAF2_coeff<-AF_DRC2$coefficients[9]
WWDAF2_coeff<-AF_DRC2$coefficients[10]
LWNAF2_coeff<-AF_DRC2$coefficients[11]
LWDAF2_coeff<-AF_DRC2$coefficients[12]

### Acropora cytherea

#Run model
WWNAF2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Nursery",],fct = LL.3())


### WW donor AF

WWDAF2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="WW_Donor",],fct = LL.3())


###leeward nursery AF

LWNAF2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Nursery",],fct = LL.3())



### leeward donor

#Run model
LWDAF2<- drm(PAM ~ Temp, data = Hold_all[Hold_all$Site_Source=="LW_Donor",],fct = LL.3())



AF2_coeffs<-data.frame(WWNAF2_coeff, WWDAF2_coeff, LWNAF2_coeff, LWDAF2_coeff)

WWNAF2_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWNAF2_pred = as.data.frame(predict(WWNAF2, newdata = WWNAF2_preddata, interval = 'confidence'))
WWNAF2_preddata = data.frame(WWNAF2_preddata, fvfm = WWNAF2_pred$Prediction, Lower = WWNAF2_pred$Lower, Upper = WWNAF2_pred$Upper)


WWDAF2_preddata = data.frame(temp = seq(30,42, length.out = 100))
WWDAF2_pred = as.data.frame(predict(WWDAF2, newdata = WWDAF2_preddata, interval = 'confidence'))
WWDAF2_preddata = data.frame(WWDAF2_preddata, fvfm = WWDAF2_pred$Prediction, Lower = WWDAF2_pred$Lower, Upper = WWDAF2_pred$Upper)


LWNAF2_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWNAF2_pred = as.data.frame(predict(LWNAF2, newdata = LWNAF2_preddata, interval = 'confidence'))
LWNAF2_preddata = data.frame(LWNAF2_preddata, fvfm = LWNAF2_pred$Prediction, Lower = LWNAF2_pred$Lower, Upper = LWNAF2_pred$Upper)

LWDAF2_preddata = data.frame(temp = seq(30,42, length.out = 100))
LWDAF2_pred = as.data.frame(predict(LWDAF2, newdata = LWDAF2_preddata, interval = 'confidence'))
LWDAF2_preddata = data.frame(LWDAF2_preddata, fvfm = LWDAF2_pred$Prediction, Lower = LWDAF2_pred$Lower, Upper = LWDAF2_pred$Upper)



#### PLOT  ####
levels(Hold_all$Site)


AF_DRC_plot2<- ggplot() +
  geom_jitter(data = Hold_all, aes(x = Temp, y = PAM, color = Site_Source), size = 1.5, width = 0.4) +
  scale_x_continuous(limits=c(30,42), breaks=c(31,34,37,40)) +
  scale_y_continuous(limits=c(0.0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  
  geom_line(data = LWDAF2_preddata, aes(x = temp, y = fvfm), color = 'red', show.legend = FALSE) +
  geom_ribbon(data = LWDAF2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'red', linetype=2, alpha = 0.2) +
  geom_vline(data = AF2_coeffs, aes(xintercept = LWDAF2_coeff), color = 'red', show.legend = FALSE) +
  geom_text(data = AF2_coeffs, aes(label = round(LWDAF2_coeff, digits=2), size = 30.0), x = 41.0, y = 0.75, show.legend = FALSE, color = 'red') +
  
  geom_line(data = LWNAF2_preddata, aes(x = temp, y = fvfm), color = 'gold', show.legend = FALSE) +
  geom_ribbon(data = LWNAF2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'gold', linetype=2, alpha = 0.2) +
  geom_vline(data = AF2_coeffs, aes(xintercept = LWNAF2_coeff), color = 'gold' , show.legend = FALSE) +
  geom_text(data = AF2_coeffs, aes(label = round(LWNAF2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.72, show.legend = FALSE, color = 'gold') +  
  
  geom_line(data = WWDAF2_preddata, aes(x = temp, y = fvfm), color = '#66ccfe', show.legend = FALSE) +
  geom_ribbon(data = WWDAF2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = '#66ccfe', linetype=2, alpha = 0.2) +
  geom_vline(data = AF2_coeffs, aes(xintercept = WWDAF2_coeff), color = '#66ccfe', show.legend = FALSE) +
  geom_text(data = AF2_coeffs, aes(label = round(WWDAF2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.66, show.legend = FALSE, color = '#66ccfe') +
  
  geom_line(data = WWNAF2_preddata, aes(x = temp, y = fvfm), color = 'blue', show.legend = FALSE) +
  geom_ribbon(data = WWNAF2_preddata, aes(x = temp, ymin=Lower, ymax=Upper), color = 'blue', linetype=2, alpha = 0.2) +
  geom_vline(data = AF2_coeffs, aes(xintercept = WWNAF2_coeff), color = 'blue', show.legend = FALSE) +
  geom_text(data = AF2_coeffs, aes(label = round(WWNAF2_coeff, digits = 2), size = 30.0), x = 41.0, y = 0.69, show.legend = FALSE, color = 'blue') +
  

  scale_color_manual(values=c('#e92000', 'gold', 'darkblue', '#66ccfe')) +
  ylab("") +
  xlab("") +  
  theme_light()+
  theme(legend.position = 'none') +
  theme(text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16))


AF_DRC_plot2


#combine indiviudla DRC curves

figure <- ggarrange(AF_DRC_plot, AF_DRC_plot2, AC_DRC_plot, AC_DRC_plot2,
                    labels = c("A", "B", "C", "D", size = 24.0),
                    ncol = 2, nrow = 2)

figure

#wrap_elements(textGrob('LW Nursery', color='red'))
