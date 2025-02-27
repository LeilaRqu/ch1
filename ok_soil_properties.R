####Old Kernza####
####Changes in soil properties over time/across treatments####

#load packages
library(ggplot2)
library(tidyverse)
library(nlme)
library(lme4)
library(stringr)

#set wd
setwd('/Users/leilarquibi/Desktop/oldkernza')

####Aggregates####
#read in data
aggdata <- read.csv('okaggs.csv')
aggdata[c('block', 'plot')] <-str_split_fixed(aggdata$plot, '', 2)

aggdata$year <- as.numeric(aggdata$year)
aggdata$mwd <- as.numeric(aggdata$mwd)
aggdata$Lmacro <- as.numeric(aggdata$Lmacro)
aggdata$Smacro <- as.numeric(aggdata$Smacro)
aggdata$micro <- as.numeric(aggdata$micro)
aggdata$block <- as.factor(aggdata$block)
aggdata$plot <- as.factor(aggdata$plot)

#plot MWD for all trts over time
ggplot(aggdata, aes(x=year, y=mwd, color=trt))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE)+
  theme_bw()

#plot MWD for crop trts over time
agg_crop <- aggdata %>% filter(trt != 'prairie' &
                                 trt != 'restored' &
                                 trt != 'NA') %>% 
  filter(mwd!='NA')

ggplot(agg_crop, aes(x=year, y=mwd, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

####model change in MWD over time-designing LMM####
#start w full model
LMagg <- lm(data=agg_crop, mwd ~ year + block + month + trt)
summary(LMagg)
par(mfrow = c(2,2))
plot(LMagg)

LMagg2 <- lm(data=agg_crop, mwd ~ trt + year/month)
summary(LMagg2)
par(mfrow = c(2,2))
plot(LMagg2)

LMagg3 <- lm(data=agg_crop, mwd ~ year/month)
summary(LMagg3)
par(mfrow = c(2,2))
plot(LMagg3)

LMagg4 <- lme(data=agg_crop, mwd~trt+year/month,
              random=~1|block, method = 'ML')
summary(LMagg4)
plot(LMagg4)

LMagg5 <- lme(data=agg_crop, mwd~year/month,
              random=~1|block, method = 'ML')
summary(LMagg5)
plot(LMagg5) #this model works best!! lowest AIC/BIC, all variables sig

####Using best model to also look at changes in agg fracs####
#Lmacro
lmaggL1 <- lme(data=agg_crop, Lmacro~trt*year/month,
              random=~1|block, method = 'ML')
summary(lmaggL1) #trt not significant

lmaggL2 <- lme(data=agg_crop, Lmacro~year/month,
               random=~1|block, method = 'ML')
summary(lmaggL2) #month not significant

lmaggL3 <- lme(data=agg_crop, Lmacro~year,
               random=~1|block, method = 'ML')
summary(lmaggL3) #best model

ggplot(data=agg_crop, aes(x=year, y=Lmacro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()
#Smacro
lmaggS1 <- lme(data=agg_crop, Smacro~trt*year/month,
               random=~1|block, method = 'ML')
summary(lmaggS1) #trt alone not sig, but sig intxns

lmaggS2 <- lme(data=agg_crop, Smacro~trt:year/month,
               random=~1|block, method = 'ML')
summary(lmaggS2) #significant interactive effects for all, except wheat:yr:mo

lmaggS3 <- lme(data=agg_crop, Smacro~year/month,
               random=~1|block, method = 'ML')
summary(lmaggS3) #lower BIC than lmaggs2, technically best fit

ggplot(data=agg_crop, aes(x=year, y=Smacro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

#micro
lmaggmi1 <- lme(data=agg_crop, micro ~ trt*year/month,
                random=~1|block, method = 'ML')
summary(lmaggmi1)

lmaggmi2 <- lme(data=agg_crop, micro ~ trt:year/month,
                random=~1|block, method = 'ML')
summary(lmaggmi2) #everything significant here, better AIC/BIC

ggplot(data=agg_crop, aes(x=year, y=micro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

####Final year comparison####
#mwd between treatments
agg23 <- aggdata %>% filter(year=='2023') %>% na.omit()
aggaov <- aov(data = agg23, mwd ~ trt)
summary(aggaov)

ggplot(agg23, aes(x=trt, y=mwd, fill=trt))+
  geom_boxplot()
#%L between treatments
Laov <- aov(data = agg23, Lmacro ~ trt)
summary(Laov)

ggplot(agg23, aes(x=trt, y=Lmacro, fill=trt))+
  geom_boxplot()
#%M between treatments
Saov <- aov(data = agg23, Smacro ~ trt)
summary(Saov)

ggplot(agg23, aes(x=trt, y=Smacro, fill=trt))+
  geom_boxplot()
#%S between treatments
miaov <- aov(data=agg23, micro ~ trt)
summary(miaov)

ggplot(agg23, aes(x=trt, y=micro, fill=trt))+
  geom_boxplot()

####Nutrients over time####
data <- read.csv('SampleData.csv')
data[c('block', 'plot')] <-str_split_fixed(data$plot, '', 2)
data$block <- as.factor(data$block)
crop <- data %>% filter(treatment=='wheat'|treatment=='KZ'|treatment=='KA')

#####poxc####
cropc <- crop %>% filter(poxc != 'NA')
ggplot(data=cropc, aes(x=year, y=poxc, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

lmc1 <- lme(data=cropc, poxc ~ treatment*year/month,
            random =~1|block)
summary(lmc1) #no crop effect

lmc2 <- lme(data=cropc, poxc~treatment:year/month,
            random =~1|block)
summary(lmc2)

lmc3 <- lme(data=cropc, poxc~treatment:year,
            random =~1|block)
summary(lmc3) #best fit

####no3####
cropni <- crop %>% filter(no3 != 'NA')
ggplot(data=cropni, aes(x=year, y=no3, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

lmni1 <- lme(data=cropni, no3~treatment*year/month,
             random=~1|block)
summary(lmni1) #nothng sig

lmni2 <- lme(data=cropni, no3~treatment:year/month,
             random=~1|block)
summary(lmni2)

lmni3 <- lme(data=cropni, no3~treatment:year,
             random=~1|block)
summary(lmni3) #all sig

####phosphate####
cropp <- crop %>% filter(po4 != 'NA')
ggplot(data=cropp, aes(x=year, y=po4, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

lmp1 <- lme(data=cropp, po4~treatment*year/month,
            random=~1|block)
summary(lmp1) #nothing sig

lmp2 <- lme(data=cropp, po4~treatment:year/month,
            random=~1|block)
summary(lmp2)#month not sig

lmp3 <- lme(data=cropp, po4~treatment:year,
            random=~1|block)
summary(lmp3) #all sig, but same value for every trt

lmp4 <- lme(data=cropp, po4~year,
            random=~1|block)
summary(lmp4) #lower bic

####ammonium####
cropa <- crop %>% filter(nh4!='NA')
ggplot(data=cropa, aes(x=year, y=nh4, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

lma1 <- lme(data=cropa, nh4~year,
            random=~1|block)
summary(lma1)

lma2 <- lme(data=cropa, nh4~year/month,
            random=~1|block)
summary(lma2) #higher aic/bic

lma3 <- lme(data=cropa, nh4~treatment*year,
            random=~1|block)
summary(lma3) #no sig

lma4 <- lme(data=cropa, nh4~treatment:year,
            random=~1|block)
summary(lma4) #significant interactions but higher aic, and same effect size for every trt??

####5 yr comparison####
cropc23 <- cropc %>% filter(year=='2023')
aovc <- aov(data=cropc23, poxc~treatment)
summary(aovc) #no sig diff

cropni23 <- cropni %>% filter(year=='2023')
aovni <- aov(data=cropni, no3~treatment)
summary(aovni) #no sig diff

cropp23 <- cropp %>% filter(year=='2023')
aovp <- aov(data=cropp23, po4~treatment)
summary(aovp) #no sig diff

cropa23 <- cropa %>% filter(year=='2023')
aova <- aov(data=cropa23, nh4~treatment)
summary(aova) #no sig diff
