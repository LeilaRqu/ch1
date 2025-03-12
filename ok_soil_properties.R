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

####model change in MWD, aggfracs over time with LMM####
agg_crop <- aggdata %>% filter(trt != 'prairie' &
                                 trt != 'restored' &
                                 trt != 'NA') %>% 
  filter(mwd!='NA')

mwdmodel <- lme(data=agg_crop, mwd~trt*year,
                random=~1|plot/month, 
                method='ML')
summary(mwdmodel)

lmacmodel <- lme(data=agg_crop, Lmacro~trt*year,
                 random=~1|plot/month, 
                 method='ML')
summary(lmacmodel)

smacmodel <- lme(data=agg_crop, Smacro~trt*year,
                 random=~1|plot/month, 
                 method='ML')
summary(smacmodel)

micmodel <- lme(data=agg_crop, micro~trt*year,
                random=~1|plot/month, 
                method='ML')
summary(micmodel)

####Plot changes in MWD, agg fracs
pmwd <- ggplot(agg_crop, aes(x=year, y=mwd, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()
plmac <- ggplot(data=agg_crop, aes(x=year, y=Lmacro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()
psmac <- ggplot(data=agg_crop, aes(x=year, y=Smacro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()
pmic <- ggplot(data=agg_crop, aes(x=year, y=micro, color=trt))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

ggarrange(nrow=2, ncol=2, common.legend = TRUE,
          lmwd, plmac, psmac, pmic)

####Final year comparison####
#mwd between treatments
agg23 <- aggdata %>% filter(year=='2023') %>% na.omit()
aggaov <- aov(data = agg23, mwd ~ trt)
summary(aggaov)

pmwd <- ggplot(agg23, aes(x=trt, y=mwd, fill=trt))+
  geom_boxplot()
#%L between treatments
Laov <- aov(data = agg23, Lmacro ~ trt)
summary(Laov)

pL <- ggplot(agg23, aes(x=trt, y=Lmacro, fill=trt))+
  geom_boxplot()
#%M between treatments
Saov <- aov(data = agg23, Smacro ~ trt)
summary(Saov)

pM <- ggplot(agg23, aes(x=trt, y=Smacro, fill=trt))+
  geom_boxplot()
#%S between treatments
miaov <- aov(data=agg23, micro ~ trt)
summary(miaov)

pS <-ggplot(agg23, aes(x=trt, y=micro, fill=trt))+
  geom_boxplot()

ggarrange(ncol=2, nrow=2, pmwd, pL, pM, pS,
          common.legend = TRUE)

####Nutrients over time####
data <- read.csv('SampleData.csv')
crop <- data %>% filter(treatment=='wheat'|treatment=='KZ'|treatment=='KA')

#poxc
cropc <- crop %>% filter(poxc != 'NA')
cmodel <- lme(data = cropc, poxc ~ treatment*year,
              random=~1|plot/month)
summary(cmodel)

pc <- ggplot(data=cropc, aes(x=year, y=poxc, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

#nh4
cropa <- crop %>% filter(nh4!='NA')
amodel <- lme(data = cropa, nh4 ~ treatment*year,
              random=~1|plot/month)
summary(amodel)

pa <- ggplot(data=cropa, aes(x=year, y=nh4, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

#no3
cropni <- crop %>% filter(no3 != 'NA')
nimodel <- lme(data = cropni, no3 ~ treatment*year,
              random=~1|plot/month)
summary(nimodel)

pni <- ggplot(data=cropni, aes(x=year, y=no3, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

#po4
cropp <- crop %>% filter(po4 != 'NA')
pmodel <- lme(data = cropp, po4 ~ treatment*year,
              random=~1|plot/month)
summary(pmodel)

pp <- ggplot(data=cropp, aes(x=year, y=po4, color=treatment))+
  geom_point(aes(shape=as.factor(month)))+
  geom_smooth(method='lm')+
  theme_bw()

#plot all together
ggarrange(nrow=2, ncol=2,
          common.legend = TRUE,
          pc, pa, pni, pp)

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

#plot y5 data
data23 <- data %>% filter(year=='2023')
bc <- ggplot(data23, aes(x=treatment, y=poxc, fill=treatment))+
  geom_boxplot()
ba <- ggplot(data23, aes(x=treatment, y=nh4, fill=treatment))+
  geom_boxplot()
bni <- ggplot(data23, aes(x=treatment, y=no3, fill=treatment))+
  geom_boxplot()
bp <- ggplot(data23, aes(x=treatment, y=po4, fill=treatment))+
  geom_boxplot()
ggarrange(nrow=2, ncol=2,
          common.legend=TRUE,
          bc, ba, bni, bp)

####Ignore below this#####
#start w full model-----disregard this. practice! no model selection, formulate
#according to your hypothesis, then see what is/isn't sig. That's the result!
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



#micro
lmaggmi1 <- lme(data=agg_crop, micro ~ trt*year/month,
                random=~1|block, method = 'ML')
summary(lmaggmi1)

lmaggmi2 <- lme(data=agg_crop, micro ~ trt:year/month,
                random=~1|block, method = 'ML')
summary(lmaggmi2) #everything significant here, better AIC/BIC

####Nutrients over time####
data <- read.csv('SampleData.csv')
data[c('block', 'plot')] <-str_split_fixed(data$plot, '', 2)
data$block <- as.factor(data$block)
crop <- data %>% filter(treatment=='wheat'|treatment=='KZ'|treatment=='KA')

#####poxc####
cropc <- crop %>% filter(poxc != 'NA')
pc <- ggplot(data=cropc, aes(x=year, y=poxc, color=treatment))+
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
pni <- ggplot(data=cropni, aes(x=year, y=no3, color=treatment))+
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
pp <- ggplot(data=cropp, aes(x=year, y=po4, color=treatment))+
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
pa <- ggplot(data=cropa, aes(x=year, y=nh4, color=treatment))+
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

