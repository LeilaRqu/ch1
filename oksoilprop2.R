####New draft of analyzing OK soil properties####
library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd('/Users/leilarquibi/Desktop/oldkernza')
data <- read.csv('SampleData.csv')

view(data)
data <- data[,1:20]
data$sample <- as.factor(data$sample)
data$type <- as.factor(data$type)
data$plot <- as.factor(data$plot)
data$treatment <- as.factor(data$treatment)
data$lmacro <- as.numeric(data$lmacro)
data$smacro <- as.numeric(data$smacro)
data$micro <- as.numeric(data$micro)
data$mwd <- as.numeric(data$mwd)
data$percentOM = as.numeric(gsub('%','', data$percentOM)) 
data$percentOC <- as.numeric(gsub('%','', data$percentOC))
data$DOC <- as.numeric(data$DOC)
data$MBC <- as.numeric(data$MBC)

####Plotting change in soil properties over time between trts####
#plot NO3 over time
no3_plot <- ggplot(data=data, aes(x=year, y=no3, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot nh4 over time
nh4_plot <- ggplot(data=data, aes(x=year, y=nh4, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot po4 over time
po4_plot <- ggplot(data=data, aes(x=year, y=po4, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot mwd over time
mwd_plot <- ggplot(data=data, aes(x=year, y=mwd, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot lmacro over time
lmacro_plot <- ggplot(data=data, aes(x=year, y=lmacro, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot smacro over time
smacro_plot <- ggplot(data=data, aes(x=year, y=smacro, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot micro over time
micro_plot <- ggplot(data=data, aes(x=year, y=micro, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot percent OM over time
percentOM_plot <- ggplot(data=data, aes(x=year, y=percentOM, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot percent OC over time
percentOC_plot <- ggplot(data=data, aes(x=year, y=percentOC, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot DOC over time
DOC_plot <- ggplot(data=data, aes(x=year, y=DOC, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot MBC over time
MBC_plot <- ggplot(data=data, aes(x=year, y=MBC, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot POXC over time
poxc_plot <- ggplot(data=data, aes(x=year, y=poxc, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot AMF hyphal length over time
amflength_plot <- ggplot(data=data, aes(x=year, y=AMFlength.prelim, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot sapro hyphal length over time
saprolength_plot <- ggplot(data=data, aes(x=year, y=sapro.hyphae, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

#plot total hyphal length over time
totalhyphae_plot <- ggplot(data=data, aes(x=year, y=total.hyphae, group=treatment, color=treatment))+
  geom_point()+
  stat_smooth(method=lm)+
  theme_classic()

##Combined Plots##
#plot nutrients together
ggarrange(nrow=1, ncol=3, 
          common.legend = TRUE,
          no3_plot, nh4_plot, po4_plot)

#plot aggs together
ggarrange(nrow=2, ncol = 2,
          common.legend = TRUE,
          mwd_plot, lmacro_plot,
          smacro_plot, micro_plot)

#plot carbon metrics together
ggarrange(nrow=2, ncol=2,
          common.legend = TRUE,
          percentOC_plot, DOC_plot,
          poxc_plot, MBC_plot)

#plot hyphal length data together
ggarrange(nrow=1, ncol=3,
          common.legend=TRUE,
          amflength_plot, saprolength_plot, totalhyphae_plot)
