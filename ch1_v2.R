########################################################
####Analysis of soil properties and fungal community####
########################################################

#install and load packages
if(!require('BiocManger', quietyl = TRUE))
  install.packages('BiocManager')

packages <- c('phyloseq', 'stringr', 'tidyverse', 'genefilter',
              'vegan', 'lattice', 'ape', 'ALDEx2', 'devtools', 'ggpubr')
BiocManager::install(packages)
lapply(packages, require, character.only = TRUE)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#setwd
setwd("/Users/leilarquibi/Desktop/microbiomes")

####Create Phyloseq Object####
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

#add alpha diversity metrics to phyloseq sample data
richness <- estimate_richness(phyloseq)
sample_data(phyloseq)$shannon <- richness$Shannon
sample_data(phyloseq)$chao1 <- richness$Chao1
sample_data(phyloseq)$ace <- richness$ACE
sample_data(phyloseq)$simpson <- richness$Simpson
sample_data(phyloseq)$fisher <- richness$Fisher

#add sequencing depth to sample data
sample_data(phyloseq)$depth <- sample_sums(phyloseq)

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

####Data exploration####
#taxonomic makeup of samples
plot_bar(phyloseq,fill='Phylum')
plot_bar(phyloseq_thresholded,fill='Phylum')

#taxonomic makeup of samples with phyla as proportions of community
phyloseq.thresh.prop <- transform_sample_counts(phyloseq.thresh, function(asv) asv/sum(asv))
p <- plot_bar(phyloseq.thresh.prop, fill='Phylum')
desired_order = c(1:47, 49:83)
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)

####Alpha Diversity####
#boxplots of various alpha diversity metrics
pshan <- ggplot(data=sample_data(phyloseq), 
                aes(x=treatment, y=shannon, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year')

pchao <- ggplot(data=sample_data(phyloseq), 
                aes(x=treatment, y=chao1, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year')

pace <- ggplot(data=sample_data(phyloseq), 
                aes(x=treatment, y=ace, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year')

psimp <- ggplot(data=sample_data(phyloseq), 
                aes(x=treatment, y=simpson, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year')

pfish <- ggplot(data=sample_data(phyloseq), 
                aes(x=treatment, y=fisher, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year')

ggarrange(pshan, pchao, pace, psimp, pfish, common.legend = TRUE)

#linear models of alpha diversity metrics for crop treatments
physeq.crop <- subset_samples(phyloseq, 
                              treatment != 'restored' & treatment != 'prairie')
sdcrop <- as(sample_data(physeq.crop), 'data.frame')

lmshan <- lm(shannon~year+treatment+depth, data=sdcrop)
summary(lmshan)

lmchao <- lm(chao1~year+treatment+depth, data=sdcrop)
summary(lmchao)

lmace <- lm(ace~year+treatment+depth, data=sdcrop)
summary(lmace)

lmsimp <- lm(ace~year+treatment+depth, data=sdcrop)
summary(lmace)

lmfish <- lm(fisher~year+treatment+depth, data=sdcrop)
summary(lmfish)

#plot
cshan <- ggplot(data=sample_data(physeq.crop), 
                aes(x=year, y=shannon, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

cchao <- ggplot(data=sample_data(physeq.crop), 
                aes(x=year, y=chao1, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

cace <- ggplot(data=sample_data(physeq.crop), 
               aes(x=year, y=ace, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

csimp <- ggplot(data=sample_data(physeq.crop), 
                aes(x=year, y=simpson, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

cfish <- ggplot(data=sample_data(physeq.crop), 
                aes(x=year, y=fisher, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

ggarrange(cshan, cchao, cace, csimp, cfish, common.legend = TRUE)

####Whole Fungal Community Composition####
##Whole data set
#clr whole data set
#input: QC'd ASV table (thresholded, but before zero removal) as matrix (rows=ASVs)
#pull out thresholded asv table
asv_thresholded <- as(otu_table(phyloseq.thresh), 'matrix')
treatment <-  sample_data(phyloseq.thresh)$treatment %>%  as.character
asv_thresholded_clr_aldex <- aldex.clr(t(asv_thresholded), 
                                       conds=treatment, 
                                       mc.samples=128,
                                       denom='all', useMC=TRUE)

#define function to extract median CLR-transformed values across MC instances
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

asv_thresholded_clr_aldex_median <- extractMedianCLRs(asv_thresholded_clr_aldex) #output: CLR transformed ASV table ready to be put back in phyloseq obj

#generate CLR transformed ASV table
asvs.clr <- otu_table(asv_thresholded_clr_aldex_median, taxa_are_rows = TRUE)

#create new phyloseq with CLR values
phyloseq.clr <- phyloseq(asvs.clr, tax_table(phyloseq.thresh), sample_data(phyloseq.thresh))

#create ordination object
ordinate(phyloseq.clr, distance='euclidean',method='PCoA') -> phyloseq.pca

#plot pca of whole data set
plot_ordination(phyloseq.clr, phyloseq.pca, type='samples', color='treatment')+
  facet_wrap(~year)+
  stat_ellipse(level=0.95)

#Permanova to test for effect of crop trt, year, month, on whole fungal community#
asvs.clr <- as(otu_table(phyloseq.clr),'matrix')
smd.clr <- as(sample_data(phyloseq.clr), 'data.frame')

set.seed(8)
asvs.clr <- t(asvs.clr)
physeq.clr.euc.adonis <- adonis2(asvs.clr~treatment+year+month, data=smd.clr, method='euclidean', by='margin')
physeq.clr.euc.adonis
physeq.clr.euc.adonis.int <- adonis2(asvs.clr~treatment*year*month, data=smd.clr, method='euclidean', by='margin')
physeq.clr.euc.adonis.int

##2023 samples only
physeq23 <- subset_samples(phyloseq.thresh, year=='2023')

asv23 <- as(otu_table(physeq23), 'matrix') #pull out otu table
treatment <- sample_data(physeq23)$treatment %>%  as.character()
asv23clraldex <- aldex.clr(t(asv23), conds=treatment, 
                      mc.samples=128, denom='all',useMC=TRUE)

asv23clr.median <- extractMedianCLRs(asv23clraldex)

asv23clrtable <- otu_table(asv23clr.median, taxa_are_rows = TRUE)

physeq23.clr <- phyloseq(asv23clrtable, tax_table(physeq23), sample_data(physeq23))

ordinate(physeq23.clr, distance='euclidean', method = 'PCoA') -> physeq23.pca

plot_ordination(physeq23.clr, physeq23.pca, type='samples', color='treatment')+
  stat_ellipse(level=0.95)

#permanova across treatments
asvclr23 <- as(otu_table(physeq23.clr), 'matrix')
smdclr23 <- as(sample_data(physeq23.clr), 'data.frame')

set.seed(8)
asvclr23 <- t(asvclr23)

physeq23.clr.euc.adonis <- adonis2(asvclr23~treatment, data=smdclr23, method='euclidean', by='margin')
physeq23.clr.euc.adonis

#pairwise permanovas to ID which treatments differ from one another
physeq23.dist <- vegdist(asvclr23, method='euclidean')
pairwise.perm23 <- pairwise.adonis2(physeq23.dist~treatment, data=smdclr23)
pairwise.perm23

####AMF####
##Richness
#subset original phyloseq object to only include AMF
physeq.amf <- subset_taxa(phyloseq, Phylum=="p__Glomeromycota")

#add AMF richness data
amf.rich <- estimate_richness(physeq.amf)
sample_data(physeq.amf)$richness <- amf.rich$Observed

#plot richness
ggplot(data=sample_data(physeq.amf), aes(x=treatment, y=richness, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(fill='Year', y='amf richness')

#lm of richness over time
physeq.amf.crop <- physeq.amf %>% 
  subset_samples(treatment != 'prairie' & treatment != 'restored')

sd.amf.crop <- as(sample_data(physeq.amf.crop), 'data.frame')

lm.amf.rich <- lm(richness ~ year+treatment+depth, data=sd.amf.crop)
summary(lm.amf.rich) #not sure how to interpret significant p-val for wheat- lower than other trts?

#plot
ggplot(data=sample_data(physeq.amf.crop), aes(x=year, y=richness, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

##Relative abundance
#agglomerate at phylum level
phyloseq.p <- tax_glom(phyloseq.thresh, taxrank=rank_names(phyloseq.clr)[2])

#clr
asvp <- as(otu_table(phyloseq.p), 'matrix')
treatment <- sample_data(phyloseq.p)$treatment %>% as.character()
asvpclraldex <- aldex.clr(t(asvp), conds=treatment,
                          mc.samples=128, denom='all', useMC=TRUE)
asvpclr.median <- extractMedianCLRs(asvpclraldex)
asvpclrtable <- otu_table(asvpclr.median, taxa_are_rows=TRUE)
physeqp.clr <- phyloseq(asvpclrtable, tax_table(phyloseq.p), sample_data(phyloseq.p))

#subset for AMF
physeqp.amf <- physeqp.clr %>% subset_taxa(Phylum=='p__Glomeromycota')

relamf <- psmelt(physeqp.amf)

#plot rel abundance across treatments
ggplot(data=relamf, aes(x=treatment, y=Abundance, fill=as.character(year)))+
  geom_boxplot(position=position_dodge(1))+
  labs(y="Relative Abundance of AMF", legend='Year')
  
#lm rel abundance over time
relamf.crop <- relamf %>% filter(treatment != 'prairie' & treatment != 'restored')
lm.relamf <- lm(Abundance ~ year + treatment + depth, data=relamf.crop)
summary(lm.relamf)

#plot
ggplot(relamf.crop, aes(x=year, y=Abundance, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

####Reactive C####
#poxc over time
lm.poxc <- lm(poxc ~ year + treatment, data=sdcrop)
summary(lm.poxc)

ggplot(data=sdcrop, aes(x=year, y=poxc, fill=treatment))+
  geom_point(aes(color=treatment))+
  geom_smooth(method='lm')

#poxc v diversity
lmshp <- lm(poxc ~ shannon, data=sdcrop)
summary(lmshp)
lmcp <- lm(poxc ~ chao1, data=sdcrop)
summary(lmcp)
lmap <- lm(poxc ~ ace, data=sdcrop)
summary(lmap)
lmsip <- lm(poxc ~ simpson, data=sdcrop)
summary(lmsip)
lmfp <- lm(poxc ~ fisher, data=sdcrop)
summary(lmfp)

#plot
shanpox <- ggplot(data=sdcrop, aes(x=shannon, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')

chaopox <- ggplot(data=sdcrop, aes(x=chao1, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')

acepox <- ggplot(data=sdcrop, aes(x=ace, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')

simppox <- ggplot(data=sdcrop, aes(x=simpson, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')

fishpox <- ggplot(data=sdcrop, aes(x=fisher, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')

ggarrange(shanpox, chaopox, acepox, simppox, fishpox)

#poxc v amf richness
lmarp <- lm(poxc ~ richness, data=sd.amf.crop)
summary(lmarp)

ggplot(data=sd.amf.crop, aes(x=richness, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x='AMF richness')

#poxc v amf rel abundance
lmarap <- lm(poxc ~ Abundance, data=relamf.crop)
summary(lmarap)

ggplot(data=relamf.crop, aes(x=Abundance, y=poxc))+
  geom_point()+
  geom_smooth(method='lm')+
  labs(x='Relative Abundance of AMF')
