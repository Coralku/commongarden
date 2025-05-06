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
library(car)

#Load file A. cytherea

DN <- read_excel("Common Garden Stats_Ac.xlsx")

DN$Geno<- as.factor(DN$Geno)
DN$Conditions<- as.factor(DN$Condition)
DN$Site<- as.factor(DN$Site)


str(DN)

#Subset data
Control_check1<-subset(DN,  Temp=="31")
Control_check<-subset(Control_check1,Timepoint=="Hold")

DN_no_field<-subset(DN, Temp=="31" | Temp=="35" | Temp=="37" | Temp=="40")
DN_no_field$Condition<- as.factor(DN_no_field$Condition)
DN_no_field$Time<- as.factor(DN_no_field$Time)

Hold_all<-subset(DN_no_field, Timepoint=="Hold")

#Stats model
PAM_hold<-lmer(PAM ~ Temp + Condition + (1|Time) + (1|Site), data=Hold_all)
summary(PAM_hold)
anova(PAM_hold)


drop1(PAM_hold)
plot(fitted(PAM_hold), residuals(PAM_hold))


residuals <-resid(PAM_hold)
qqnorm(residuals)
qqline(residuals)


dev.off()

sjPlot::plot_model(PAM_hold, type="diag")
step(PAM_hold, reduce.random=FALSE)
PAM_hold_final<-lmer(PAM ~ Temp + Condition + (1|Time) + (1|Site),data=Hold_all)
anova(PAM_hold_final)
