####F26 Committee meeting Figures####
library(tidyverse)
library(ggpubr)
library(phyloseq)
library(adespatial)
library(corrplot)
library(vegan)

setwd('/Users/leilarquibi/Desktop/oldkernza')

data <- read.csv('SampleData.csv')
data$mwd <- as.numeric(data$mwd)
data$treatment <- as.factor(data$treatment)
data$lmacro <- as.numeric(data$lmacro)
data$smacro <- as.numeric(data$smacro)
data$micro <- as.numeric(data$micro)
data$DOC <- as.numeric(data$DOC)
data$MBC <- as.numeric(data$MBC)
data$poxc <- as.numeric(data$poxc)

data$percentOM <- gsub('%', '', data$percentOM)
data$percentOM <- as.numeric(data$percentOM)

data.lms <- data %>% filter(treatment!='restored' & treatment!='prairie')
data.lms <- droplevels(data.lms)

data <- data %>% filter(treatment!='prairie')
data.23 <- data %>% filter(year=='2023')
data.23 <- droplevels(data.23)

lmcols <- c('olivedrab3', 'gold', 'red3')
col23 <- c('olivedrab3', 'gold', 'mediumorchid','red3')
bar_order <- c('wheat', 'KZ', 'KA', 'restored')

####Soil Health Ish!####
####Nutrients Stats/Figures####
no3lm <- lm(no3~year+treatment, data=data.lms)
summary(no3lm)
pno3lm <- ggplot(data.lms, aes(x=year, y=no3, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic() +
  labs(y='NO3 (ug/g soil)')+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

nh4lm <- lm(nh4~year+treatment, data=data.lms)
summary(nh4lm)
pnh4lm <- ggplot(data.lms, aes(x=year, y=nh4, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic() +
  labs(y='NH4 (ug/g soil)')+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

po4lm <- lm(po4~year+treatment, data=data.lms)
summary(po4lm)
ppo4lm <- ggplot(data.lms, aes(x=year, y=po4, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic() +
  labs(y='PO4 (ug/g soil)')+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

aovno3 <- aov(no3~treatment, data=data.23)
summary(aovno3)
TukeyHSD(aovno3)
sumno3 <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(no3),
            stdev = sd(no3),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pno3bar <- ggplot(data=sumno3, aes(x=treatment, y=mean, color=treatment))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('color'))

aovnh4 <- aov(nh4~treatment, data=data.23)
summary(aovnh4)
sumnh4 <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(nh4),
            stdev = sd(nh4),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pnh4bar <- ggplot(data=sumnh4, aes(x=treatment, y=mean, color=treatment))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('color'))

aovpo4 <- aov(po4~treatment, data=data.23)
summary(aovpo4)
TukeyHSD(aovpo4)
sumpo4 <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(po4),
            stdev = sd(po4),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
ppo4bar <- ggplot(data=sumpo4, aes(x=treatment, y=mean, color=treatment))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('color'))

ggarrange(nrow=3, ncol=2, common.legend = TRUE,
          pno3lm, pno3bar,
          pnh4lm, pnh4bar,
          ppo4lm, ppo4bar)


####Aggregates Figures####
mwdlm <- lm(mwd~treatment*year, data = data.lms)
summary(mwdlm)
pmwdlm <- ggplot(data.lms, aes(x=year, y=mwd, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic()+
  labs(y='Mean\nWeight\nDiameter (mm)')+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

lmaclm <- lm(lmacro~treatment*year, data = data.lms)
summary(lmaclm)
plmaclm <- ggplot(data.lms, aes(x=year, y=lmacro, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic()+
  labs(y="Percent Large\nMacroaggregates\n(>2mm)")+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

smaclm <- lm(smacro~treatment*year, data = data.lms)
summary(smaclm)
psmaclm <- ggplot(data.lms, aes(x=year, y=smacro, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic() +
  labs(y="Percent Small\nMacroaggregates\n(2mm>x>250um)")+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

miclm <- lm(micro~treatment*year, data = data.lms)
summary(miclm)
pmiclm <- ggplot(data.lms, aes(x=year, y=micro, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method = 'lm', se=FALSE)+
  theme_classic() +
  labs(y='Percent\nMicroaggregates\n(250um>x>53um)')+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

aovmwd <- aov(mwd~treatment, data=data.23)
summary(aovmwd)
summwd <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(mwd)),
            stdev = sd(na.omit(mwd)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pmwdbar <- ggplot(data=summwd, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       x='treatment',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

aovlmac <- aov(lmacro~treatment, data=data.23)
summary(aovlmac)
sumlmac <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(lmacro)),
            stdev = sd(na.omit(lmacro)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
plmacbar <- ggplot(data=sumlmac, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       x='treatment',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

aovsmac <- aov(smacro~treatment, data=data.23)
summary(aovsmac)
sumsmac <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(smacro)),
            stdev = sd(na.omit(smacro)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
psmacbar <- ggplot(data=sumsmac, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       x='treatment',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

aovmicro <- aov(micro~treatment, data=data.23)
summary(aovmicro)
summic <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(micro)),
            stdev = sd(na.omit(micro)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pmicbar <- ggplot(data=summic, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(y='',
       x='treatment',
       title = '2023')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

ggarrange(nrow=4, ncol=2, common.legend = TRUE,
          pmwdlm, pmwdbar,
          plmaclm, plmacbar,
          psmaclm, psmacbar,
          pmiclm, pmicbar)

####C Figures####
omlm <- lm(percentOM~treatment*year, data = data.lms)
summary(omlm)
pomlm <- ggplot(data=data.lms, aes(x=year, y=percentOM, color=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm', se=FALSE)+
  labs(y='Percent Organic Matter')+
  theme_classic()+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

poxlm <- lm(poxc~treatment*year, data = data.lms)
summary(poxlm)
ppoxlm <- ggplot(data=data.lms, aes(x=year, y=poxc, color = treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm', se=FALSE)+
  labs(y='POXC (ug/g soil)')+
  theme_classic()+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

doclm <- lm(DOC~treatment*year, data = data.lms)
summary(doclm)
pdoclm <- ggplot(data=data.lms, aes(x=year, y=DOC, color = treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm', se=FALSE)+
  labs(y='DOC (ug/g soil)')+
  theme_classic()+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

mbclm <- lm(MBC~treatment*year, data = data.lms)
summary(mbclm)
pmbclm <- ggplot(data=data.lms, aes(x=year, y=MBC, color = treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm', se=FALSE)+
  labs(y='MBC (ug/g soil)')+
  theme_classic()+
  xlim(2019,2023)+
  scale_color_manual(values=lmcols,
                     aesthetics=c('color'))

omaov <- aov(percentOM~treatment, data=data.23)
summary(omaov)
TukeyHSD(omaov)
omsum <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(percentOM)),
            stdev = sd(na.omit(percentOM)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pombar <- ggplot(data=omsum, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(title='2023',
       x='treatment',
       y='')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

poxaov <- aov(poxc~treatment, data=data.23)
summary(poxaov)
TukeyHSD(poxaov)
poxsum <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(poxc)),
            stdev = sd(na.omit(poxc)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
ppoxbar <- ggplot(data=poxsum, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(title='2023',
       x='treatment',
       y='')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

docaov <- aov(DOC~treatment, data=data.23)
summary(docaov)
TukeyHSD(docaov)
docsum <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(DOC)),
            stdev = sd(na.omit(DOC)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pdocbar <- ggplot(data=docsum, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(title='2023',
       x='treatment',
       y='')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

mbcaov <- aov(MBC~treatment, data=data.23)
summary(mbcaov)
TukeyHSD(mbcaov)
mbcsum <- data %>% filter(year=='2023') %>% 
  group_by(treatment) %>% #to get a mean value for each species
  summarize(mean = mean(na.omit(MBC)),
            stdev = sd(na.omit(MBC)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pmbcbar <- ggplot(data=mbcsum, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity',
           aes(x=factor(treatment, level=bar_order)))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(title='2023',
       x='treatment',
       y='')+
  scale_color_manual(values=col23,
                     aesthetics=c('fill'))

ggarrange(nrow=4, ncol=2, common.legend = TRUE,
          pomlm, pombar,
          ppoxlm, ppoxbar,
          pdoclm, pdocbar,
          pmbclm, pmbcbar)

####Microbial Community Get into it YUH####
ps <- readRDS('phyloseq.RDS')
ps.clr <- readRDS('phyloseq.clr.RDS')
view(sample_data(ps.clr))

sample_data(ps)$shannon <- estimate_richness(ps)$Shannon
pscrop <- subset_samples(ps, treatment!='prairie' & treatment!='restored')
pshanlm <- ggplot(sample_data(pscrop), aes(x=year, y=shannon, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(aes(fill=treatment), method='lm')+
  ylab('Soil Fungal Diversity (Shannon)')+
  theme_classic()

shannon <- as.matrix(sample_data(ps))
shannon <- as.data.frame(shannon)
shannon$shannon <- as.numeric(shannon$shannon)
shannon$year <- as.numeric(shannon$year)
shancrop <- shannon %>% filter(treatment!='prairie' & treatment!='restored')
shanlm <- lm(shannon~year*treatment, data=shancrop)
summary(shanlm)

shan.23 <- shannon %>% filter(treatment != 'prairie' & year=='2023')
aovshan <- aov(shannon~treatment, data=shan.23)
summary(aovshan)
TukeyHSD(aovshan)
shansum <- shannon %>% filter(treatment != 'prairie' & year=='2023') %>% 
  group_by(treatment) %>%
  summarize(mean = mean(na.omit(shannon)),
            stdev = sd(na.omit(shannon)),
            sample_size = n(),
            se = stdev/sqrt(sample_size))
pshanbar <- ggplot(data=shansum, aes(x=treatment, y=mean, fill=treatment))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position='identity',
                width=0.2)+
  theme_classic()+
  labs(title='2023',
       y='')

ggarrange(nrow=1, ncol=2, common.legend = TRUE,
          pshanlm, pshanbar)
####PCoA####
#read in for next time
phyloseq.clr <- readRDS('phyloseq.clr.RDS')
#remove prairie plots bc they're nonsense
phyloseq.clr <- subset_samples(phyloseq.clr, treatment!='prairie')

#create ordination object
ordinate(phyloseq.clr, distance='euclidean',method='PCoA') -> phyloseq.clr.euclid.pcoa

#PCoA plot of whole fungal community dissimilarity between all samples using Aitchison distance, colored by treatment
sample_data(phyloseq.clr)$year <- as.factor(sample_data(phyloseq.clr)$year)
levels(sample_data(phyloseq.clr)$year)

ptrt <- plot_ordination(phyloseq.clr, phyloseq.clr.euclid.pcoa, 
                type='samples', 
                color='treatment', 
                shape = 'year')+
  geom_point(size=2.5)+
  stat_ellipse(aes(group=treatment))+
  theme_classic()+
  labs(title='Treatment')+
  scale_color_manual(values=col23,
                     aesthetics='color')
pyr <- plot_ordination(phyloseq.clr, phyloseq.clr.euclid.pcoa, 
                       type='samples', 
                       color='treatment', 
                       shape = 'year')+
  geom_point(size=2.5)+
  stat_ellipse(aes(group=year))+
  theme_classic()+
  labs(title='Year')+
  scale_color_manual(values=col23,
                     aesthetics='color')

ggarrange(ncol=2, nrow=1, common.legend=TRUE, legend='right',
          ptrt, pyr)

#PERMANOVA
asvs.clr <- as(otu_table(phyloseq.clr),'matrix')
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')
smd.clr$treatment <- as.factor(smd.clr$treatment)

set.seed(888)
asvs.clr <- t(asvs.clr)
perm <- how(within=Within(type='free'), 
            plots = Plots(strata=smd.clr$plot), 
            nperm=888)
physeq.clr.euc.adonis1 <- adonis2(asvs.clr~year*treatment,
                                  strata=smd.clr$plot, #look at perm instead potentially (permute package, set permutations)
                                  data=smd.clr,
                                  method='euclidean', 
                                  by='margin') #interaction is not significant
another.one <- adonis2(asvs.clr~year+treatment,
                       strata=smd.clr$plot, #look at perm instead potentially (permute package, set permutations)
                       data=smd.clr,
                       method='euclidean', 
                       by='margin') #year p=0.001, trt p=0.002
#Pairwise
pwa <- pairwise.adonis2(asvs.clr~treatment+year, nperm = perm, data=smd.clr,
                        method='euclidean', by='margin')
pwa
#adjust pvals for multiple comparisons
pvals <- data.frame(comp = c('KAvKZ', 'KAvWheat', 'KAvRestored',
                             'KZvWheat', 'KZvRestored', 'WheatvRestored'),
                    p=c(pwa$KA_vs_KZ$`Pr(>F)`[1],
                        pwa$KA_vs_wheat$`Pr(>F)`[1],
                        pwa$KA_vs_restored$`Pr(>F)`[1],
                        pwa$KZ_vs_wheat$`Pr(>F)`[1],
                        pwa$KZ_vs_restored$`Pr(>F)`[1],
                        pwa$wheat_vs_restored$`Pr(>F)`[1]))
pvals

####RDA####
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')
smd.lim <- smd.clr %>% dplyr::select('treatment','poxc', 'no3', 'po4', 'nh4', 'mwd', 'percentOM', 'DOC', 'MBC') 
#make sure variables are coded correctly
smd.lim$poxc <- as.numeric(smd.lim$poxc)
smd.lim$mwd <- as.numeric(smd.lim$mwd)
smd.lim$percentOM <- gsub('%', '', smd.lim$percentOM)
smd.lim$percentOM <- as.numeric(smd.lim$percentOM)
smd.lim$DOC <- as.numeric(smd.lim$DOC)
#complete cases only
smd.lim <- na.omit(smd.lim)
keep <- row.names(smd.lim)
#make asv samples match complete cases
ps.clr.complete <- subset_samples(phyloseq.clr, sample_data(phyloseq.clr)$sample %in% keep)
asvs.clr <- as(otu_table(ps.clr.complete), 'matrix')
asvs.clr <- t(asvs.clr)
#check for multicollinearity between explanatory variables
m = cor(smd.lim)
corrplot(m)

#only numeric/soil variables in the test df
smd.lim.num <- smd.lim %>% dplyr::select(!'treatment')
mod0 <- rda(asvs.clr ~ 1, smd.lim.num)
mod1 <- rda(asvs.clr ~ ., smd.lim.num)
adjR2.fullmod <- RsquareAdj(mod1)$adj.r.squared

mod <- ordiR2step(mod0, scope=formula(mod1), 
                  R2scope=TRUE)
#fwd selected model
mod$call 
mod.sel <- rda(asvs.clr ~ mwd + DOC + percentOM + no3, data = smd.lim)
RsquareAdj(mod.sel) #Adj R2 = 0.07395224
anova.cca(mod.sel, step = 1000, by = 'term') #p=0.001 whole model, all variables also significant
anova.cca(mod.sel, step = 1000, by = 'axis') #RDA 1, 2, and 3 are significant

levels(as.factor(smd.lim$treatment))
colvec <- c('olivedrab3', 'gold', 'mediumpurple', 'red3')
plot(mod.sel, type='none', scaling = 2)
points(mod.sel, display='sites',
       col=colvec[as.factor(smd.lim$treatment)],
       pch = 19)
text(mod.sel, display='bp')
legend('topright', legend = levels(as.factor(smd.lim$treatment)),
       pch =19,
       bty='n',
       col=colvec)

#check my work
ordiplot(mod.sel, scaling = 2, type = 'text')

####Carbon RDA####
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')
smd.c <- smd.clr %>% dplyr::select('treatment','poxc', 'mwd', 
                                   'percentOM', 'DOC', 'MBC') 
#make sure variables are coded correctly
smd.c$mwd <- as.numeric(smd.c$mwd)
smd.c$percentOM <- gsub('%', '', smd.c$percentOM)
smd.c$percentOM <- as.numeric(smd.c$percentOM)
smd.c$DOC <- as.numeric(smd.c$DOC)

#complete cases only
smd.c <- na.omit(smd.c)
keepc <- row.names(smd.c)
#make asv samples match complete cases
ps.c.complete <- subset_samples(phyloseq.clr, sample_data(phyloseq.clr)$sample %in% keepc)
asvs.c <- as(otu_table(ps.c.complete), 'matrix')
asvs.c <- t(asvs.c)

#rda model
cmod <- rda(asvs.c ~ mwd + DOC + percentOM + MBC + poxc, data = smd.c)

RsquareAdj(cmod) #Adj R2 = 0.07346645
anova.cca(cmod, step = 1000)
anova.cca(cmod, step = 1000, by = 'term') #p=0.001 mwd, DOC, percentOM, p=0.11 MBC, p=0.2 POXC
anova.cca(cmod, step = 1000, by = 'axis') #RDA 1, 2, and 3 are significant

#plot that ish
levels(as.factor(smd.c$treatment))
colvec <- c('olivedrab3', 'gold', 'mediumpurple', 'red3')
plot(cmod, type='none', scaling = 2)
points(cmod, display='sites',
       col=colvec[as.factor(smd.c$treatment)],
       pch = 19)
text(cmod, display='bp')
legend('topright', legend = levels(as.factor(smd.c$treatment)),
       pch =19,
       bty='n',
       col=colvec)
