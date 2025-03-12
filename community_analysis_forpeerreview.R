#######################################################################
####Analysis of alpha and beta diversity of soil fungal communities####
#######################################################################

####Set up####
#load packages
library(phyloseq)
library(stringr)
library(tidyverse)
library(genefilter)
library(vegan)
library(lattice)
library(ape)
library(ALDEx2)
library(pairwiseAdonis)

#setwd
setwd("/Users/leilarquibi/Desktop/oldkernza")

####Phyloseq basics####
#load in ASV table
asvs <- readRDS("seqtab.nochim.RDS")

#load in taxonomy table
tax <- readRDS("taxa_fungi.RDS")

#load in sample metadata
smd <- read.csv("SampleData.csv")

#format data for phyloseq
asvs <- otu_table(asvs, taxa_are_rows = FALSE) 
tax <- tax_table(tax)
smd <- sample_data(smd)

#reformat sample names to match between components
sample_names(smd)
sample_names(asvs)
sample_names(smd) <- smd$sample
sample_names(asvs) <- str_remove(sample_names(asvs), 'filtered_') %>% 
  str_remove('_S1_L001_R1_001.fastq')

#combine data into single phyloseq object
phyloseq <- phyloseq(asvs, tax, smd)

#give ASVs short nicknames
for (i in 1:ntaxa(phyloseq)) { #iteration for each asv
  #create blank dataframe during 1st iteration to store ASV nicknames
  if (i==1) {ASVkey=data.frame('ASV'=character(), 'Sequence'=character())}
  #each following iteration changes the integer value of i
  nickname <- paste0('asv', as.character(i)) #store nickname as asv + integer
  seq <- taxa_names(phyloseq)[i] #pull sequence for the asv (should this be i?)
  thisASV <- data.frame('ASV'=nickname, 'Sequence'=seq) #connect nickname & sequence for this ASV
  ASVkey <- rbind(ASVkey,thisASV) #add this row to ASVkey dataframe, grows w each iteration
  taxa_names(phyloseq)[i] <- nickname #rename ASV in phyloseq object
  if (i==ntaxa(phyloseq)) {saveRDS(ASVkey, 'ASVkey.RDS')} #once every ASV is nicknamed, save key to RDS file
}

#check that names have been changed correctly
taxa_names(phyloseq)

#save phyloseq object to RDS
saveRDS(phyloseq, file = 'phyloseq.RDS')

####Quality control####
#rarefaction curve to confirm adequate sampling effort
rarecurve(as(otu_table(phyloseq), 'matrix'), step=100, label=FALSE)

#threshold dataset, ie remove non-reproducible ASVs
threshold <- kOverA(k=10, A=1) #remove ASVs that appear less than 10 times in at least 1 sample
phyloseq.thresh <- filter_taxa(phyloseq, threshold, prune=TRUE)

#check that thresholding reduced number of taxa observed
ntaxa(phyloseq)
ntaxa(phyloseq.thresh)


####Data exploration
#taxonomic makeup of samples
plot_bar(phyloseq,fill='Phylum')
plot_bar(phyloseq_thresholded,fill='Phylum')

#taxonomic makeup of samples with phyla as proportions of community
phyloseq.prop <- transform_sample_counts(phyloseq, function(asv) asv/sum(asv))
phyloseq.thresh.prop <- transform_sample_counts(phyloseq.thresh, function(asv) asv/sum(asv))
p <- plot_bar(phyloseq.thresh.prop, fill='Phylum', facet_grid=~treatment)
desired_order = c(1:47, 49:83)
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)

pspka <- subset_samples(phyloseq.thresh.prop, treatment=='KA')
pspkz <- subset_samples(phyloseq.thresh.prop, treatment=='KZ')
pspp <- subset_samples(phyloseq.thresh.prop, treatment=='prairie')
pspr <- subset_samples(phyloseq.thresh.prop, treatment=='restored')
pspw <- subset_samples(phyloseq.thresh.prop, treatment=='wheat')

pka <- plot_bar(pspka, fill='Phylum', facet_grid =  ~year)
pkz <- plot_bar(pspkz, fill='Phylum', facet_grid =  ~year)
pp <- plot_bar(pspp, fill='Phylum', facet_grid =  ~year)
pr <- plot_bar(pspr, fill='Phylum', facet_grid =  ~year)
pw <- plot_bar(pspw, fill='Phylum', facet_grid =  ~year)

####Alpha diversity of whole fungal community#####
#calculate alpha diversity of each sample
estimate_richness(phyloseq)

#plot alpha diversity across samples (non-thresholded count data) by year
phyloseq.crop <- phyloseq %>% 
  subset_samples(treatment %in% c('KA', 'KZ', 'wheat'))
richness <- estimate_richness(phyloseq.crop)
richness$sample <- row.names(richness)
richness$sample <- str_replace(richness$sample, "X", "")
year.shannon <- as.factor(sample_data(phyloseq.crop)$year)
pchao <- plot_richness(phyloseq.crop, x='year', measures=c('Chao1'))+
  geom_boxplot(fill=NA, aes(group=year, color=as.factor(year)))+
  facet_wrap(~treatment)+
  ylab('Species richness of soil fungal community (Chao)')
pshannon<- plot_richness(phyloseq.crop, x='year', measures=c('Shannon'))+
  geom_boxplot(fill=NA, aes(group=year, color=as.factor(year)))+
  facet_wrap(~treatment)+
  ylab('Diversity of soil fungal community (Shannon)')
ggarrange(nrow=1, ncol=2, common.legend=TRUE, 
          pchao, pshannon)

#mixed linear model of richness~year
cropsd <- as.matrix(sample_data(phyloseq.crop))
cropsd <- as.data.frame(cropsd)
cropsd$year <- as.numeric(cropsd$year)
cropsd$chao1 <- as.numeric(cropsd$chao1)
richmodel <- lme(data=cropsd, chao1 ~ year*treatment,
                 random = ~1|plot/month)
summary(richmodel)

#mixed linear model of shannon~year
cropsd$shannon <- as.numeric(cropsd$shannon)
shanmodel <- lme(data=cropsd,shannon ~ year*treatment,
                 random = ~1|plot/month)
summary(shanmodel)

#plot alpha diversity across treatments in 2023
physeq.2023 <- subset_samples(phyloseq, year=="2023")
chao23 <- plot_richness(physeq.2023, x='treatment', measures=c('Chao1'))+
  geom_boxplot(fill=NA, aes(color=treatment))+
  ylab('Chao richness of soil fungal community')
shan23 <- plot_richness(physeq.2023, x='treatment', measures=c('Shannon'))+
  geom_boxplot(fill=NA, aes(color=treatment))+
  ylab('Shannon diversity of soil fungal community')
ggarrange(nrow=1, ncol=2, common.legend=TRUE,
          chao23, shan23)

#aov of rich/div across trt in 2023
ps23 <- as.matrix(sample_data(physeq.2023))
ps23 <- as.data.frame(ps23)

aovc23 <- aov(data=ps23, chao1~treatment)
summary(aovc23)
TukeyHSD(aovc23)

aovs23 <- aov(data=ps23, shannon~treatment)
summary(aovs23)
TukeyHSD(aovs23)

#####Beta diversity of whole fungal community####
#input: QC'd ASV table (thresholded, but before zero removal) as matrix (rows=ASVs)
#pull out thresholded asv table
asv_thresholded <- as(otu_table(phyloseq.thresh), 'matrix')
treatment <-  sample_data(phyloseq.thresh)$treatment %>%  as.character
asv_thresholded_clr_aldex <- aldex.clr(t(asv_thresholded), conds=treatment, mc.samples=128,
                                       denom='all',useMC=TRUE)
#function to extract median CLR-transformed values across MC instances
extractMedianCLRs <- function(clr.obj){
  for (sample in getSampleIDs(clr.obj)) {
    MCinstances <- getMonteCarloInstances(clr.obj)[[sample]]
    medianCLRs <- apply(MCinstances,MARGIN = 1, function(x) median(x))
    if (sample == getSampleIDs(clr.obj)[1]) {
      medianCLRs.all <- data.frame(medianCLRs)
      colnames(medianCLRs.all) <- sample}
    else {medianCLRs.all[[sample]] <- medianCLRs}}
  return(medianCLRs.all)
}

asv_thresholded_clr_aldex_median <- extractMedianCLRs(asv_thresholded_clr_aldex)
#^output: CLR transformed ASV table ready to be put back in phyloseq obj

#create new phyloseq object with CLR transformed ASV table
asvs.clr <- otu_table(asv_thresholded_clr_aldex_median, taxa_are_rows = TRUE)
#make asv nicknames match between asvs.clr and tax
for (i in 1:ntaxa(tax)) { #iteration for each asv
  #create blank dataframe during 1st iteration to store ASV nicknames
  if (i==1) {ASVkey.tax=data.frame('ASV'=character(), 'Sequence'=character())}
  #each following iteration changes the integer value of i
  nickname <- paste0('asv', as.character(i)) #store nickname as asv + integer
  seq <- taxa_names(tax)[i] 
  thisASV <- data.frame('ASV'=nickname, 'Sequence'=seq) #connect nickname & sequence for this ASV
  ASVkey.tax <- rbind(ASVkey.tax,thisASV) #add this row to ASVkey dataframe, grows w each iteration
  taxa_names(tax)[i] <- nickname #rename ASV in phyloseq object
  if (i==ntaxa(tax)) {saveRDS(ASVkey.tax, 'ASVkey.RDS')} #once every ASV is nicknamed, save key to RDS file
}
phyloseq.clr <- phyloseq(asvs.clr, tax, smd)

#create ordination object
ordinate(phyloseq.clr, distance='euclidean',method='PCoA') -> phyloseq.clr.euclid.pcoa

#PCoA plot of whole funal community dissimilarity between all samples using Aitchison distance, colored by treatment
sample_data(phyloseq.clr)$year <- as.factor(sample_data(phyloseq.clr)$year)
levels(sample_data(phyloseq.clr)$year)

plot_ordination(phyloseq.clr, phyloseq.clr.euclid.pcoa, type='samples', 
                color='treatment', shape='year')

plot_ordination(phyloseq.clr, phyloseq.clr.euclid.pcoa, type='samples', color='treatment')+
  facet_wrap(~year)+
  stat_ellipse(level=0.95)

  
#Permanova to test for effect of crop trt, year, month, on whole fungal community#
asvs.clr <- as(otu_table(phyloseq.clr),'matrix')
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')

set.seed(8)
asvs.clr <- t(asvs.clr)
physeq.clr.euc.adonis1 <- adonis2(asvs.clr~year*treatment,
                                 strata=smd.clr$plot:smd.clr$month,
                                 data=smd.clr,
                                 method='euclidean', 
                                 by='margin') #all sig

physeq.clr.euc.adonis2 <- adonis2(asvs.clr~year+treatment+month,
                                  strata=smd.clr$plot,
                                  data=smd.clr,
                                  method='euclidean', 
                                  by='margin') #all sig

#PCoA and permanova of 2023 samples
#subset 2023 samples from clr phyloseq object, pull out ASV table and sample data
phyloseq.clr.23 <- subset_samples(phyloseq.clr, year=='2023') #this is wrong, need, to recalc distance matrix for 23 only
asv.clr.23 <- as(otu_table(phyloseq.clr.23), 'matrix')
smd.clr.23 <- as(sample_data(phyloseq.2023), 'data.frame')

#perform permanova
set.seed(8)
asv.clr.23 <- t(asv.clr.23)
physeq23.clr.euc.adonis <- adonis2(asv.clr.23~treatment, data=smd.clr.23, method='euclidean', by='margin')

#plot pcoa
ordinate(phyloseq.clr.23, distance='euclidean', method='PCoA')-> physeq23.clr.euc.pcoa
plot_ordination(phyloseq.clr.23, physeq23.clr.euc.pcoa, type='samples', color='treatment')+
  stat_ellipse(level=0.95)

#conduct pairwise permanovas to ID which treatments differ from one another
#calculate distance matrix
aitchison.phyloseq.23 <- vegdist(asv.clr.23, method='euclidean')
pairwise.perm.23 <- pairwise.adonis2(aitchison.phyloseq.23~treatment, data=smd.clr.23)
#function returns adjusted p-values

####Response of AMF over time to cropping treatment####
#Compare richness of AMF species from year to year within crop treatments
#start with non-thresholded count data, subset by taxa
amf <- subset_taxa(phyloseq, Phylum=="p__Glomeromycota")
amf.crop<- amf %>% subset_samples(treatment %in% c('KA', 'KZ', 'wheat'))
year.amf.crop <- as.factor(sample_data(amf.crop)$year)

amfrich <- plot_richness(amf.crop, x='year', measures='Observed')+
  geom_boxplot(fill=NA, aes(group=year.amf.crop, 
                            color=year.amf.crop))+
  facet_wrap(~treatment)+
    ylab('Number of ASVs: AMF')

#linear mixed model for AMF richness over time
sample_data(amf.crop)$amfrich <- estimate_richness(amf.crop, measures='Observed')
amfrichsd <- as.matrix(sample_data(amf.crop))
amfrichsd <- as.data.frame(amfrichsd)
amfrichsd$amfrich <- as.numeric(amfrichsd$amfrich)
amfrichsd$year <- as.numeric(amfrichsd$year)

amfrichmodel <- lme(data=amfrichsd, amfrich ~ treatment*year,
                    random = ~1|plot/month)
summary(amfrichmodel)

#aov of richness ~ treatment in 2023 samples
amf.23 <- subset_samples(amf, year=='2023')
richness.amf.23 <- estimate_richness(amf.23, measures='Observed')
aov.amf.23 <- aov(richness.amf.23$Observed ~ sample_data(amf.23)$treatment)
summary(aov.amf.23)

amfrich23 <- plot_richness(amf.23, x='treatment', measures='Observed')+
  geom_boxplot(fill=NA, aes(group=treatment, color=treatment))+
  ylab('Number of ASVs: AMF')

#Compare relative abundance of AMF between years within cropping treatments
#starting with clr/thresholded phyloseq object
#combine taxa at phylum level
phyloseq.clr.p <- tax_glom(phyloseq.clr, taxrank=rank_names(phyloseq.clr)[2])
amf <- subset_taxa(phyloseq.clr.p, Phylum=='p__Glomeromycota')

amf.df <- psmelt(amf)
head(amf.df)
amf.df.crop <- amf.df %>% 
  filter(treatment %in% c('KA', 'KZ', 'wheat'))
amf.df.crop$year <- as.factor(amf.df.crop$year)


ggplot(data=amf.df.crop, aes(x=year, y=Abundance, group=year, col=year))+
  geom_boxplot(fill=NA)+
  facet_wrap(~treatment)+
  ylab('Relative abundance of AMF')

#anovas to check for sig differences in amf abundance from year to year within each crop trt
amf.df.KA <- amf.df %>% filter(treatment=="KA")
amf.KA.aov <- aov(amf.df.KA$Abundance ~ amf.df.KA$year)


amf.df.KZ <- amf.df %>% filter(treatment=="KZ") 
amf.KZ.aov <- aov(amf.df.KZ$Abundance ~ amf.df.KZ$year)

amf.df.wheat <- amf.df %>% filter(treatment=="wheat")
amf.wheat.aov <- aov(amf.df.wheat$Abundance ~ amf.df.wheat$year)

#anova of AMF abundance between treatments in 2023
amf.df.23 <- amf.df %>% 
  filter(year=='2023')

summary(aov(amf.df.23$Abundance ~ amf.df.23$treatment))

#plot amf abundance between trts in 2023
ggplot(amf.df.23, aes(x=treatment, y=Abundance, col=treatment))+
  geom_boxplot(fill=NA)+
  ylab('Relative abundance of AMF')


####trying with just proportions, not clr
phyloseq.prop <- tax_glom(phyloseq.prop, taxrank=rank_names(phyloseq.prop)[2])
amfp <- subset_taxa(phyloseq.prop, Phylum=='p__Glomeromycota')
amfpcrop <- subset_samples(amfp, treatment != 'prairie' 
                           & treatment != 'restored')

amfp.df <- psmelt(amfp)
amfp.crop <- psmelt(amfpcrop)
head(amfp.df)

#lmm prop AMF over time
amfpropmodel <- lme(data=amfp.crop, Abundance ~ treatment*year,
                    random= ~1|plot/month)
summary(amfpropmodel)

#aov differences in prop AMF yr 5
amfp.df.23 <- amfp.df %>% 
  filter(year=='2023')

amfp23model <- aov(data=amfp.df.23, Abundance ~ treatment)
summary(amfp23model)

#plot amf prop over time and year 5
amfprop <- ggplot(data=amfp.crop, aes(x=year, y=Abundance, group=year))+
  geom_point()+
  geom_boxplot(fill=NA, aes(color=as.factor(year)))+
  facet_wrap(~treatment)+
  ylab('Relative abundance of AMF')

amfprop23 <- ggplot(data=amfp.df.23, aes(x=treatment, y=Abundance))+
  geom_boxplot(fill=NA, aes(color=treatment))+
  ylab('Relative abundance of AMF')

#plot amfrichnes and prop over time together
ggarrange(nrow=1, ncol=2,
          common.legend = TRUE,
          amfrich, amfprop)

#plot amf rich and prop y5 together
ggarrange(nrow=1, ncol=2, common.legend = TRUE,
          amfrich23, amfprop23)



