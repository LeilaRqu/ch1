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
setwd("/Users/leilarquibi/Desktop/microbiomes")

####Phyloseq basics####
#load in ASV table
asvs <- readRDS("dada2_output/seqtab.nochim.RDS")

#load in taxonomy table
tax <- readRDS("dada2_output/taxa_fungi.RDS")

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
phyloseq.thresh.prop <- transform_sample_counts(phyloseq.thresh, function(asv) asv/sum(asv))
p <- plot_bar(phyloseq.thresh.prop, fill='Phylum')
desired_order = c(1:47, 49:83)
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)

####Alpha diversity of whole fungal community#####
#calculate alpha diversity of each sample
richness <- estimate_richness(phyloseq)

#plot alpha diversity across samples (non-thresholded count data) by year
phyloseq.crop <- phyloseq %>% 
  subset_samples(treatment %in% c('KA', 'KZ', 'wheat'))
year.shannon <- as.factor(sample_data(phyloseq.crop)$year)
plot_richness(phyloseq.crop, x='year', measures=c('Shannon'), color=year.shannon)+
  geom_boxplot(fill=NA, aes(group=year))+
  facet_wrap(~treatment)+
  ylab('Shannon diversity of soil fungal community')

#lm of richness ~ year in KA
physeq.KA <- subset_samples(phyloseq, treatment=="KA")
sample_data(physeq.KA)$year <- as.numeric(sample_data(physeq.KA)$year)
richness.KA <- estimate_richness(physeq.KA)
hist(richness.KA$Shannon) #make sure data are reasonably normally distributed
lm.KA <- lm(richness.KA$Shannon ~ sample_data(physeq.KA)$year)
summary(lm.KA)

#lm of richness ~ year in KZ
physeq.KZ <- subset_samples(phyloseq, treatment=="KZ")
sample_data(physeq.KZ)$year <- as.numeric(sample_data(physeq.KZ)$year)
richness.KZ <- estimate_richness(physeq.KZ)
hist(richness.KZ$Shannon) #make sure data are reasonably normally distributed
lm.KZ <- lm(richness.KZ$Shannon ~ sample_data(physeq.KZ)$year)
summary(lm.KZ)


#lm of richness ~ year in wheat
physeq.wheat <- subset_samples(phyloseq, treatment=="wheat")
sample_data(physeq.wheat)$year <- as.numeric(sample_data(physeq.wheat)$year)
richness.wheat <- estimate_richness(physeq.wheat)
hist(richness.wheat$Shannon) #make sure data are reasonably normally distributed
glm.wheat <- lm(richness.wheat$Shannon ~ sample_data(physeq.wheat)$year)
summary(glm.wheat)


#plot alpha diversity across treatments in 2018
physeq.2018 <- subset_samples(phyloseq, year == '2018')
trt.2018 <- as.factor(sample_data(physeq.2018)$treatment)
plot_richness(physeq.2018, x= 'treatment', measures=c('Shannon'), color=trt.2018)+
  geom_boxplot(fill=NA)+
  ylab('Shannon diversity of soil fungal community')

#plot alpha diversity across treatments in 2023
physeq.2023 <- subset_samples(phyloseq, year=="2023")
trt.23 <- as.factor(sample_data(physeq.2023)$treatment)
plot_richness(physeq.2023, x='treatment', measures=c('Shannon'), color=trt.23)+
  geom_boxplot(fill=NA)+
  ylab('Shannon diversity of soil fungal community')

#aov to compare alpha diversity across treatments in 2018
richness.2018 <- estimate_richness(physeq.2018)
aov.2018 <- aov(richness.2018$Shannon ~ sample_data(physeq.2018)$treatment)
summary(aov.2018)

#aov to compare alpha diversity across treatments in 2023
richness.2023 <- estimate_richness(physeq.2023)
aov.2023 <- aov(richness.2023$Shannon ~ sample_data(physeq.2023)$treatment)
summary(aov.2023)
TukeyHSD(aov.2023)

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
plot_ordination(phyloseq.clr, phyloseq.clr.euclid.pcoa, type='samples', color='treatment')+
  facet_wrap(~year)+
  stat_ellipse(level=0.95)

#Permanova to test for effect of crop trt, year, month, on whole fungal community#
asvs.clr <- as(otu_table(phyloseq.clr),'matrix')
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')

set.seed(8)
asvs.clr <- t(asvs.clr)
physeq.clr.euc.adonis <- adonis2(asvs.clr~treatment+year+month, data=smd.clr, method='euclidean', by='margin')
physeq.clr.euc.adonis.int <- adonis2(asvs.clr~treatment*year*month, data=smd.clr, method='euclidean', by='margin')

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

plot_richness(amf.crop, x='year', measures='Observed', color=year.amf.crop)+
  geom_boxplot(fill=NA, aes(group=year.amf.crop))+
  facet_wrap(~treatment)+
    ylab('AMF species richness')

#lm of richness ~ year in KA
amf.KA <- subset_samples(amf, treatment=='KA')
sample_data(amf.KA)$year <- as.numeric(sample_data(amf.KA)$year)
richness.amf.KA <- estimate_richness(amf.KA, measures = 'Observed')
hist(richness.amf.KA$Observed) #check that distribution is reasonably normal-- very zero inflated
lm.amf.ka <- lm(richness.amf.KA$Observed ~ sample_data(amf.KA)$year)
summary(lm.amf.ka)

#lm of richness ~ year in KZ
amf.KZ <- subset_samples(amf, treatment=='KZ')
sample_data(amf.KZ)$year <- as.numeric(sample_data(amf.KZ)$year)
richness.amf.KZ <- estimate_richness(amf.KZ, measures = 'Observed')
hist(richness.amf.KZ$Observed) #check that distribution is reasonably normal-- very zero inflated
lm.amf.kz <- lm(richness.amf.KZ$Observed ~ sample_data(amf.KZ)$year)
summary(lm.amf.kz)

#lm of richness ~ year in wheat
amf.wheat <- subset_samples(amf, treatment=='wheat')
sample_data(amf.wheat)$year <- as.numeric(sample_data(amf.wheat)$year)
richness.amf.wheat <- estimate_richness(amf.wheat, measures = 'Observed')
hist(richness.amf.wheat$Observed) #check that distribution is reasonably normal-- very zero inflated
lm.amf.wheat <- lm(richness.amf.wheat$Observed ~ sample_data(amf.wheat)$year)
summary(lm.amf.wheat)

#aov of richness ~ treatment in 2023 samples
amf.23 <- subset_samples(amf, year=='2023')
richness.amf.23 <- estimate_richness(amf.23, measures='Observed')
aov.amf.23 <- aov(richness.amf.23$Observed ~ sample_data(amf.23)$treatment)
summary(aov.amf.23)

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










