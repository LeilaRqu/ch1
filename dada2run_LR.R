#!/usr/bin/Rscript

####Workflow for fastq/fasta sequnces with adapters trimmed using R on cluster
##use to quality filter, trim, denoise, merge paired reads, infer ASVs, make ASV table, 
##remove chimeras, assign taxonomy, save ASV table and taxonomy table
##Adapted from code provided by Dr. Maggie Wagner

#usage: Rscript dada2run_LR.R


library(dada2)
library(stringr)

#set and display wd
setwd("/home/ler4794/ch1/trimmed")
wdPath <- getwd()
cat("working directory is:",wdPath)


#Load trimmed reads
print("loading trimmed reads...")
trimpath <-paste0(wdPath) #create directory w fwd and rev trimmed reads
ftrimmed <- list.files(trimpath,pattern="_R1_",full.names=TRUE) #fwd reads
rtrimmed <- list.files(trimpath,pattern="_R2_",full.names=TRUE) #rev reads

#check that file names look correct
print("filenames for trimmed reads:")
print(ftrimmed)
print(rtrimmed)

#make a place for R to put the quality filtered files by rewriting the file names
library(stringr)
ffiltered <- str_replace_all(ftrimmed,'trimmed','filtered')
rfiltered <- str_replace_all(rtrimmed,'trimmed','filtered')


#execute command for quality-filtering and trimming
library(dada2)
filterOutput <- filterAndTrim(fwd=ftrimmed, filt=ffiltered,
	rev=rtrimmed, filt.rev=rfiltered,
	maxEE=2, truncQ=2, maxN=0, rm.phix=FALSE, 
	compress=TRUE, verbose=TRUE, multithread=TRUE)

#what do I need to change ab maxEE, truncQ?
#omit trunclen bc using ITS! reqd. overlap of. F/R reads excludes taxa w longer ITS, can’t use expected sequence length for QC
#output of above command: new quality-filtered and trimmed FastQ files located on cluster at the paths stored in ffiltered and rfiltered

#denoising step 1: learn error rates for fwd/rev reads
print('Begin learning error rates for forward reads...')
errf <- learnErrors(ffiltered, nbases=1e8, multithread=TRUE)
print('Begin learning error rates for reverse reads...')
errr <- learnErrors(rfiltered, nbases=1e8, multithread=TRUE)
#output: R objects errf and errr which contain estimated error rates

#denoising step 2: de-replicate filtered input fastq files
print('Begin dereplicating input FASTQ files...')
derepf <- derepFastq(ffiltered)
derepr <- derepFastq(rfiltered)
#output: R object derepF and derepR containing dereplicated sequences

#denoising step 3: use error rates to infer asvs
print('Begin using error rates to infer ASVs...')
ddf <- dada(derepf, err=errf, multithread=TRUE)
ddr <- dada(derepr, err=errr, multithread=TRUE)
#output: R objects ddf and ddr containing de-noised sequences

#merge paired reads
print('Merging paired reads...')
merger <- mergePairs(ddf, derepf, ddr, derepr)
#output: R object called merger containing merged sequences

#count sequences, make ASV table
print('Making sequence table...')
seqtab <- makeSequenceTable(merger)
#output: R object called seqtab (matrix of how many times each ASV was found in each sample)

#remove chimeras from ASV table
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#output: new ASV table w chimeric sequences removed

#taxonomic assignment
taxa_fungi <- assignTaxonomy(seqtab.nochim,"/home/ler4794/unite/unite_9.fasta", multithread=TRUE)
#output: R object called taxa_fungi, df containing taxonomic assignment of all sequences

#save necessary R objects to file
saveRDS(seqtab.nochim,'seqtab.nochim.RDS') #final ASV table
saveRDS(taxa_fungi,'taxa_fungi.RDS') #taxonomy table
