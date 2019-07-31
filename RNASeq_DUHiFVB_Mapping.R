
# INSTALL R PACKAGES

#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library("rhdf5")
#biocLite("tximport",lib="../Library/") 
library("tximport")
## DESeq
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")


###### Analysis on mapping made on classical mouse genome (C57) #############################################################################################

samples <- read.csv("samples.txt", header = TRUE, sep = "\t", row.names=1) 
# Annotation of the samples
     #Strain = FVB or DUHi
     #Jaw = upper jaw : mX or lower jaw : mD
     #Weight = 3 classes based on the weight in mg: Low, Medium, High

samples$Strain=relevel(samples$Strain,'FVB')
samples$Jaw=relevel(samples$Jaw, 'mD')
samples$Weight=relevel(samples$Weight, 'L') 

# Read gene annotations and mapping files
MusId <- read.csv("Biomart_GRCm38.p5.txt",header=TRUE,sep="\t") #EnsembleID et MGI
filesMus <- Sys.glob("Mapping/quant_Mus/SHPC*.fastq.qt/abundance.tsv") 

namesSamples=filesMus #noms des 12 echantillons
for (i in c(1:12)){
     namesSamples[i] <- unlist(strsplit(filesMus[i], "Mapping/quant_Mus/"))[2]
     namesSamples[i] <- unlist(strsplit(namesSamples[i], ".fastq.qt/abundance.tsv", fixed = TRUE))
     f=read.csv(filesMus[i], header = TRUE, sep = "\t")
     f$target_id <- sapply( f$target_id, function(x){strsplit(as.character(x), "\\.")[[1]][1]} )
     write.table(f, file = paste(filesMus[i], "1", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
}

filesMus1 = paste(filesMus, "1", sep = "") 
names(filesMus1)=namesSamples

tx2gene <- MusId[,c("Transcript.stable.ID", "MGI.symbol")]
tx2Spike <- read.csv("CaracteristicsERCC.txt", header=FALSE, sep="\t")
tx2Spike <- tx2Spike[-1,-3]
tx2Spike$V2 <- tx2Spike$V1
colnames(tx2Spike) <- colnames(tx2gene) 
tx2gene <- rbind(tx2gene,tx2Spike) #add spikes


txiMus <- tximport(filesMus1, type = "kallisto", tx2gene = tx2gene) #builds count matrix

## Tximport
txiMusCounts <- data.frame(txiMus[2]) 
colnames(txiMusCounts) <- namesSamples
txiMusCounts <- txiMusCounts[-c(1),] 
txiMusCounts <- as.matrix(round(txiMusCounts[complete.cases(txiMusCounts),]))
# 32646 MGI on classical mouse genome

ddsMusStrain <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain)
DESeqMUSStrain <- DESeq(ddsMusStrain)
DESeqMUSStrainResults <- results(DESeqMUSStrain)
table(DESeqMUSStrainResults$padj < 0.05)
# 3959 genes are DE between strains when mapping is performed on classical mouse genome

ddsMusJaw <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain + Jaw)
DESeqMUSJaw <- DESeq(ddsMusJaw)
DESeqMUSJawResults <- results(DESeqMUSJaw)
table(DESeqMUSJawResults$padj < 0.05)
# 1809 genes are DE between jaws when mapping is performed on classical mouse genome, and difference between strains is taken into account

ddsMusWeight <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain + Jaw + Weight)
DESeqMUSWeight <- DESeq(ddsMusWeight)
DESeqMUSWeightResults <- results(DESeqMUSWeight)
table(DESeqMUSWeightResults$padj < 0.05)
# 144 genes are DE between M/L stages when mapping is performed on classical mouse genome, and difference between strains and jaws is taken into account



###### Analysis on mapping made on FVB strain mouse genome #############################################################################################

FVBId <- read.csv("Biomart_FVB_MGI.txt",header=F,sep="|") #EnsembleID et MGI
names(FVBId)=c("GeneID","Transcript.stable.ID","MGI.symbol")

filesFVB <- Sys.glob("Mapping/quant_FVB/SHPC*.fastq.qt/abundance.tsv") #mappings on FVB strain

tx2geneFVB <- FVBId[,c("Transcript.stable.ID", "MGI.symbol")]
tx2geneFVB <- rbind(tx2geneFVB,tx2Spike)
tx2geneFVB$MGI.symbol=as.character(tx2geneFVB$MGI.symbol)
tx2geneFVB$MGI.symbol[tx2geneFVB$MGI.symbol==""]="NA"

for (i in c(1:12)){
     f=read.csv(filesFVB[i], header = TRUE, sep = "\t")
     ff=f
     ff$target_id <- sapply(f$target_id, function(x){strsplit(as.character(x), "|", fixed = TRUE)[[1]][2]} )
     ff$target_id[grep("ERCC",f$target_id)]=as.character(f$target_id[grep("ERCC",f$target_id)])
     #f <- f[order(f$target_id), ]
     #f <- f[tx2geneFVB$MGI.symbol!="",]
     write.table(ff, file = paste(filesFVB[i], "1", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
}
filesFVB1 = paste(filesFVB, "1", sep = "")
names(filesFVB1)=namesSamples

txiFVB <- tximport(filesFVB1, type = "kallisto", tx2gene = tx2geneFVB) 

## Tximport
txiFVBCounts <- data.frame(txiFVB[2]) # build count matrix
colnames(txiFVBCounts) <- namesSamples
txiFVBCounts <- txiFVBCounts[-c(1),] 
txiFVBCounts <- as.matrix(round(txiFVBCounts[complete.cases(txiFVBCounts),])) #round numbers
# remove counts of transcripts with no corresponding MGIs
txiFVBCounts=txiFVBCounts[row.names(txiFVBCounts)!=c("NA"),]
# 29917 genes with MGI on FVB reference

ddsFVBStrain <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain)
DESeqFVBStrain <- DESeq(ddsFVBStrain)
DESeqFVBStrainResults <- results(DESeqFVBStrain)
table(DESeqFVBStrainResults$padj < 0.05)
#  3623 genes are DE between strains when mapping is performed on FVB mouse genome


ddsFVBJaw <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain + Jaw)
DESeqFVBJaw <- DESeq(ddsFVBJaw)
DESeqFVBJawResults <- results(DESeqFVBJaw)
table(DESeqFVBJawResults$padj < 0.05)
# 1719 genes are DE between jaws when mapping is performed on classical mouse genome, and difference between strains is taken into account


ddsFVBWeight <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain + Jaw + Weight)
DESeqFVBWeight <- DESeq(ddsFVBWeight)
DESeqFVBWeightResults <- results(DESeqFVBWeight)
table(DESeqFVBWeightResults$padj < 0.05)
# 106 genes are DE between M/L stages when mapping is performed on classical mouse genome, and difference between strains and jaws is taken into account



###### Comparaison of the mappings  ###############################################################################################################
## Names of genes expressed in common

namesFVBgenes <- rownames(txiFVBCounts)
namesMusgenes <- rownames(txiMusCounts)

commongenes <- namesFVBgenes[namesFVBgenes%in%namesMusgenes]
# 19202 genes with MGI in common between the 2 genomes.

## Base mean of the genes with a mapping on FVB and classical Mus
Merge <- merge(data.frame(DESeqMUSStrainResults),data.frame(DESeqFVBStrainResults), by="row.names") 


## linear regression
plot(x=Merge$baseMean.x, y=Merge$baseMean.y, log="xy", xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="All genes expression")
abline(a=0, b=1, col="red")





###### Spike-ins counts are rather well correlated with expected concentrations  #############################################################################################################################
##  baseMean
Spikein <- read.csv("CaracteristicsERCC.txt", header=TRUE, sep="\t")
SpikegenesAll <- rownames(Spikein)
#length(SpikegenesAll)
row.names(Merge)=Merge$Row.names
Merge=Merge[,c("baseMean.x","baseMean.y")]
names(Merge)=c("MapMus","MapFVB")
Spikein <- merge(Spikein,Merge[SpikegenesAll,],by=0)
rownames(Spikein) <- Spikein$Row.names

plot(log(Spikein$concentration), log(Spikein$MapMus), ylab="log(baseMean map C57)", xlab="log(concentration)", main="Mapping on C57")
abline(a=0, b=1, col="red")

plot(log(Spikein$concentration), log(Spikein$MapFVB),  ylab="log(baseMean map FVB)", xlab="log(concentration)", main="Mapping on FVB")
abline(a=0, b=1, col="red")

cor.test(Spikein$concentration, Spikein$MapFVB)$estimate^2
cor.test(Spikein$concentration, Spikein$MapMus)$estimate^2
# the R squared are above 0.97 between expected and observed levels of expression in ERCC.



###### DESeq analyses to find genes affected by the mapping on FVB or classical MUS genome ##########################################################################################


MergeForDESeq <- merge(txiFVBCounts,txiMusCounts, by="row.names")
rownames(MergeForDESeq) <- MergeForDESeq$Row.names
MergeForDESeq <- MergeForDESeq[-1]
samplesMerge <- data.frame(colnames(MergeForDESeq))
colnames(samplesMerge) <- "mapping_strain"
rownames(samplesMerge) <- samplesMerge$mapping_strain
samplesMerge$mapping_strain <- c("map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB")

ddsMerge <- DESeqDataSetFromMatrix(MergeForDESeq, colData = samplesMerge, design = ~ mapping_strain)
DESeqMerge <- DESeq(ddsMerge)
DESeqMerge <- DESeqMerge[!is.na(DESeqMerge)]
DESeqMergeResults <- results(DESeqMerge)
table(DESeqMergeResults$padj < 0.05)
# 2234 genes are DE between strains 

SameGenes <- rownames(DESeqMergeResults)[!DESeqMergeResults$padj < 0.05]
# so we work on 16968 genes with no mapping artifact
txiSameCounts <- txiMusCounts[rownames(txiMusCounts)%in%SameGenes,]

ddsSame <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ 1)
DESeqSame <- DESeq(ddsSame)
DESeqSameResults <- results(DESeqSame)

ddsSameStrain <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Strain)
DESeqSameStrain <- DESeq(ddsSameStrain)
DESeqSameStrainResults <- results(DESeqSameStrain)
table(DESeqSameStrainResults$padj < 0.05)
# 3050 genes differ between strains

ddsSameJaw <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Strain+Jaw)
DESeqSameJaw <- DESeq(ddsSameJaw)
DESeqSameJawResults <- results(DESeqSameJaw)
table(DESeqSameJawResults$padj < 0.05)
# 1612 genes differ between jaws once the strain effect is taken into account

ddsSameJS <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Jaw + Strain)
DESeqSameJS <- DESeq(ddsSameJS)
DESeqSameJSResults <- results(DESeqSameJS)
table(DESeqSameJSResults$padj < 0.05)
# 3619 genes differ between strains once the jaw effect is taken into account

write.table(DESeqSameJawResults,file="DESeqSameJawResults.txt",sep="\t",col.names=T,row.names=T,quote=F)
write.table(DESeqSameJSResults,file="DESeqSameStrainResults.txt",sep="\t",col.names=T,row.names=T,quote=F)

bm=counts(DESeqSameJS,norm=T)
bm2=merge(DESeqSameJSResults,bm,by.x=0,by.y=0)
write.table(bm2,file="DESeqSameJSResults.txt",sep="\t",col.names=T,quote=F,row.names=F)


list=c("Axin2", "Kremen1", "Osr2", "Sfrp2","Spry2","Sostdc1")
DESeqSameJSResults[list,]
