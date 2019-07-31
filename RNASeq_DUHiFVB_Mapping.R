
# INSTALL PACKAGES

#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library("rhdf5")
#biocLite("tximport",lib="../Library/") 
library("tximport")
## DESeq
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")


###### DESeq des donnees de mapping sur Mus (C57) #############################################################################################
samples <- read.csv("samples.txt", header = TRUE, sep = "\t", row.names=1) #conditions de chaque echantillon
     #Strain = souche
     #Row = mx ou md
     #Weight = stade developpemental

samples$Strain=relevel(samples$Strain,'FVB')
samples$Weight=relevel(samples$Weight, 'Medium')
samples$Weight=relevel(samples$Weight, 'Low') #On ordonne par ordre d'importance

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

filesMus1 = paste(filesMus, "1", sep = "") #necessaire car IDs (file): [ENSMUST00000177564.1 .... sont en .1
names(filesMus1)=namesSamples #on attribue ces noms aux differentes chaines de caractere

## Creation tx2gene
tx2gene <- MusId[,c("Transcript.stable.ID", "MGI.symbol")]
tx2Spike <- read.csv("Mapping/CaracteristicsERCC.txt", header=FALSE, sep="\t")
tx2Spike <- tx2Spike[-1,-3]
tx2Spike$V2 <- tx2Spike$V1
colnames(tx2Spike) <- colnames(tx2gene) 
tx2gene <- rbind(tx2gene,tx2Spike) #Ajout des noms des spike au tx2gene


txiMus <- tximport(filesMus1, type = "kallisto", tx2gene = tx2gene) #Association des matrices d'abondance aux MGI

## Tximport
txiMusCounts <- data.frame(txiMus[2]) #On prend que les valeurs de counts
colnames(txiMusCounts) <- namesSamples
txiMusCounts <- txiMusCounts[-c(1),] #On enleve la ligne -1 qui est la somme des comptes pour chaque echantillon
txiMusCounts <- as.matrix(round(txiMusCounts[complete.cases(txiMusCounts),]))



ddsMus <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ 1)
DESeqMUS <- DESeq(ddsMus)
DESeqMUSResults <- results(DESeqMUS)
table(DESeqMUSResults$padj < 0.05)

ddsMusStrain <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain)
DESeqMUSStrain <- DESeq(ddsMusStrain)
DESeqMUSStrainResults <- results(DESeqMUSStrain)
table(DESeqMUSStrainResults$padj < 0.05)

ddsMusJaw <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain + Jaw)
DESeqMUSJaw <- DESeq(ddsMusJaw)
DESeqMUSJawResults <- results(DESeqMUSJaw)
table(DESeqMUSJawResults$padj < 0.05)

ddsMusWeight <- DESeqDataSetFromMatrix(txiMusCounts, colData = samples, design = ~ Strain + Jaw + Weight)
DESeqMUSWeight <- DESeq(ddsMusWeight)
DESeqMUSWeightResults <- results(DESeqMUSWeight)
table(DESeqMUSWeightResults$padj < 0.05)





###### DESeq des donnees de mapping sur FVB #############################################################################################

FVBId <- read.csv("Biomart_FVB_MGI.txt",header=TRUE,sep="\t") #EnsembleID et MGI
filesFVB <- Sys.glob("Mapping/quant_FVB2/SHPC*.fastq.qt/abundance.tsv") #chemins des fichiers d'abondances

tx2geneFVB <- FVBId[,c("Transcript.stable.ID", "MGI.symbol")]
tx2geneFVB <- rbind(tx2geneFVB,tx2Spike)

#for (i in c(1:12)){
 #    f=read.csv(filesFVB[i], header = TRUE, sep = "\t")
  #   f$target_id <- sapply( f$target_id, function(x){strsplit(as.character(x), "|", fixed = TRUE)[[1]][2]} )
     #f <- f[order(f$target_id), ]
     #f <- f[tx2geneFVB$MGI.symbol!="",]
#     write.table(f, file = paste(filesFVB[i], "1", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
#}
#filesFVB1 = paste(filesFVB, "1", sep = "") #necessaire car IDs (file): [ENSMUST00000177564.1 .... sont en .1
names(filesFVB)=namesSamples

txiFVB <- tximport(filesFVB, type = "kallisto", tx2gene = tx2geneFVB) #associse les matrices d'abondance aux MGI

## Tximport
txiFVBCounts <- data.frame(txiFVB[2]) #On prend que les valeurs de counts
colnames(txiFVBCounts) <- namesSamples
txiFVBCounts <- txiFVBCounts[-c(1),] #On enleve la ligne -1 qui est la somme des comptes pour chaque echantillon
txiFVBCounts <- as.matrix(round(txiFVBCounts[complete.cases(txiFVBCounts),])) #Nombres entiers

## DESeq
ddsFVB <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ 1)
DESeqFVB <- DESeq(ddsFVB)
DESeqFVBResults <- results(DESeqFVB)
table(DESeqFVBResults$padj < 0.05)

#padj <- DESeqFVBResults$padj
#log2FC <- DESeqFVBResults$log2FoldChange
#baseMean <- DESeqFVBResults$baseMean
#genesDEMapFVB <- data.frame(baseMean, log2FC, padj)
#rownames(genesDEMapFVB) <- rownames(DESeqFVBResults)
#write.table(genesDEMapFVB, file = "Mapping/genes DE Map FVB.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

ddsFVBStrain <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain)
DESeqFVBStrain <- DESeq(ddsFVBStrain)
DESeqFVBStrainResults <- results(DESeqFVBStrain)
table(DESeqFVBStrainResults$padj < 0.05)

ddsFVBJaw <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain + Jaw)
DESeqFVBJaw <- DESeq(ddsFVBJaw)
DESeqFVBJawResults <- results(DESeqFVBJaw)
table(DESeqFVBJawResults$padj < 0.05)

ddsFVBWeight <- DESeqDataSetFromMatrix(txiFVBCounts, colData = samples, design = ~ Strain + Jaw + Weight)
DESeqFVBWeight <- DESeq(ddsFVBWeight)
DESeqFVBWeightResults <- results(DESeqFVBWeight)
table(DESeqFVBWeightResults$padj < 0.05)





###### Comparaison des mapping ###############################################################################################################
## Noms des genes exprimes en commun
namesFVBgenes <- rownames(txiFVBCounts)
namesMusgenes <- rownames(txiMusCounts)

commongenes <- namesFVBgenes[namesFVBgenes%in%namesMusgenes]


## Base mean des genes avec mapping chez FVB et mapping chez Mus
BaseMeanMus <- as.matrix(DESeqMUSResults$baseMean)
row.names(BaseMeanMus) <- rownames(DESeqMUSResults)

BaseMeanFVB <- as.matrix(DESeqFVBResults$baseMean)
row.names(BaseMeanFVB) <- rownames(DESeqFVBResults)

Merge <- merge(BaseMeanMus, BaseMeanFVB, by="row.names", all=TRUE) ## Merge = base mean pour les deux mapping
colnames(Merge) <- c("gene", "BaseMeanMus","BaseMeanFVB")
rownames(Merge) <- Merge[, 1] ## set rownames
Merge <- Merge[, -1]

Merge <- replace(Merge, Merge==0, NA) # Supp des genes au basemean a 0
Merge <- na.omit(Merge) # Supp des genes pas en commun (NA)


## Regression lineaire entre deux mapping
plot(x=Merge$BaseMeanMus, y=Merge$BaseMeanFVB, log="xy", xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="All genes expression")
abline(a=0, b=1, col="red")
lm(log(Merge$BaseMeanFVB) ~ log(Merge$BaseMeanMus))


## Ratios d'expression entre les deux mapping
Merge$MusFVB <- Merge$BaseMeanMus/Merge$BaseMeanFVB 
Merge$FVBMus <- Merge$BaseMeanFVB/Merge$BaseMeanMus
#write.table(Merge, file = "Mapping/basemean FVB Mus.txt", sep="\t", row.names=TRUE, quote=FALSE)

#hist(Merge$MusFVB[Merge$MusFVB>4], breaks=200)   #... bof, log ?
#plot(density(Merge$MusFVB[Merge$MusFVB>4])) 

#plot(Merge$MusFVB[Merge$MusFVB>5])
#idx <- identify(Merge$BaseMeanMus, Merge$BaseMeanFVB)
#rownames(ddsTxiStrainResults)[idx]

## Avoir les basemean pour chaque Strain
baseMeanStrainMus <- sapply( levels(DESeqMUS$Strain), function(lvl) rowMeans( counts(DESeqMUS, normalized=TRUE)[,DESeqMUS$Strain == lvl] ) )
baseMeanStrainFVB <- sapply( levels(DESeqFVB$Strain), function(lvl) rowMeans( counts(DESeqFVB, normalized=TRUE)[,DESeqFVB$Strain == lvl] ) )
MergeStrain <- merge(baseMeanStrainMus, baseMeanStrainFVB, by="row.names", all=TRUE)
#write.table(MergeStrain, file = "Mapping/basemean FVB Mus - Strain.txt", sep="\t", row.names=TRUE, quote=FALSE)




###### Analyse des genes DE #######################################################################################################################
## Recuperation des genes DE
DEmapMus <- rownames(txiMusCounts)[DESeqMUSStrainResults$padj < 0.05]
DEmapMus <- DEmapMus[!is.na(DEmapMus)]

DEmapFVB <- rownames(txiFVBCounts)[DESeqFVBStrainResults$padj < 0.05]
DEmapFVB <- DEmapFVB[!is.na(DEmapFVB)]



## Genes DE en commun pour les deux mapping
outersect <- function(x, y) {
     sort(c(setdiff(x, y),
            setdiff(y, x)))  }

DEAll <- unique(rbind(matrix(DEmapMus), matrix(DEmapFVB)))
MergeDEAll <- Merge[rownames(Merge)%in%DEAll,]
plot(x=MergeDEAll$BaseMeanMus, y=MergeDEAll$BaseMeanFVB, log="xy", xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="All DE genes expression")
abline(a=0, b=1, col="red")
lm(log(MergeDEAll$BaseMeanFVB) ~ log(MergeDEAll$BaseMeanMus))

DEcommon <- matrix(intersect(DEmapMus,DEmapFVB))
MergeDEcommon <- Merge[rownames(Merge)%in%DEcommon,]
plot(x=MergeDEcommon$BaseMeanMus, y=MergeDEcommon$BaseMeanFVB, log="xy", xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="Common DE genes expression")
#idx <- identify(MergeDE$BaseMeanMus, MergeDE$BaseMeanFVB)
abline(a=0, b=1, col="red")
summary(lm(log(MergeDEcommon$BaseMeanFVB) ~ log(MergeDEcommon$BaseMeanMus)))

DENoncommon <- matrix(outersect(DEmapMus,DEmapFVB))
MergeDENoncommon <- Merge[rownames(Merge)%in%DENoncommon,]
plot(x=MergeDENoncommon$BaseMeanMus, y=MergeDENoncommon$BaseMeanFVB, log="xy", xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="Non Common DE genes expression")
abline(a=0, b=1, col="red")
lm(log(MergeDENoncommon$BaseMeanFVB) ~ log(MergeDENoncommon$BaseMeanMus))




###### Analyse Spike-in #############################################################################################################################
## En baseMean
Spikein <- read.csv("Mapping/CaracteristicsERCC.txt", header=TRUE, sep="\t")
SpikegenesAll <- rownames(Spikein)
#length(SpikegenesAll)
Spikein <- merge(Spikein, Merge, by="row.names")
SpikegenesExpr <- Spikein$Row.names
#length(SpikegenesExpr)
rownames(Spikein) <- SpikegenesExpr

plot(log(Spikein$concentration), log(Spikein$BaseMeanMus), col = ifelse(Spikein$lg < 800,'blue','red'), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="Mapping on C57")
regSpikein <- lm(log(Spikein$BaseMeanMus) ~ log(Spikein$concentration))
summary(regSpikein)
abline(regSpikein)

plot(log(Spikein$concentration), log(Spikein$BaseMeanFVB), col = ifelse(Spikein$lg < 800,'blue','red'), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="Mapping on FVB")
regSpikein <- lm(log(Spikein$BaseMeanFVB) ~ log(Spikein$concentration))
summary(regSpikein)
abline(regSpikein)


## baseMean / length
plot(log(Spikein$concentration), log(Spikein$BaseMeanMus/Spikein$lg), col = ifelse(Spikein$lg < 800,'blue','red'), pch = 19, ylab="log(baseMean/lenght)", xlab="log(concentration)", main="Mapping on C57")
regSpikein <- lm(log(Spikein$BaseMeanMus/Spikein$lg) ~ log(Spikein$concentration))
summary(regSpikein)
abline(regSpikein)


## Par echantillon normalises
CountsSampleMUSNorm <- round(counts(DESeqMUS, normalized = TRUE))
SpikeinMUSNorm <- merge(CountsSampleMUSNorm, Spikein, by="row.names")
SpikeinMUSNorm <- SpikeinMUSNorm[,-c(14,18,19,20)]

#CountsSampleFVB <- round(counts(DESeqFVB, normalized = TRUE))
#SpikeinFVB <- merge(CountsSampleFVB, Spikein, by="row.names")
#SpikeinFVB <- SpikeinFVB[,-c(14,17,19,20)]

par(mfrow=c(3,4))
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC25), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC25")
lm(log(SpikeinMUSNorm$SHPC25) ~ log(SpikeinMUSNorm$concentration))
abline(a=0, b=1, col="red")
#pb de log(0)
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC26), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC26")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC27), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC27")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC28), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC28")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC29), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC29")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC30), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC30")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC31), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC31")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC32), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC32")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC33), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC33")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC34), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC34")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC35), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC35")
abline(a=0, b=1, col="red")
plot(log(SpikeinMUSNorm$concentration), log(SpikeinMUSNorm$SHPC36), pch = 19, ylab="log(baseMean)", xlab="log(concentration)", main="SHPC36")
abline(a=0, b=1, col="red")




## Pourcentage de spike-in
#sumSpike <- sum(Spikein$BaseMeanMus)
#sumAll <- sum(Merge$BaseMeanMus)
#sumSpike/sumAll

ratioSpikeNorm <- rbind(sum(SpikeinMUSNorm$SHPC25)/sum(CountsSampleMUSNorm[,"SHPC25"]),
                             sum(SpikeinMUSNorm$SHPC26)/sum(CountsSampleMUSNorm[,"SHPC26"]),
                             sum(SpikeinMUSNorm$SHPC27)/sum(CountsSampleMUSNorm[,"SHPC27"]),
                             sum(SpikeinMUSNorm$SHPC28)/sum(CountsSampleMUSNorm[,"SHPC28"]),
                             sum(SpikeinMUSNorm$SHPC29)/sum(CountsSampleMUSNorm[,"SHPC29"]),
                             sum(SpikeinMUSNorm$SHPC30)/sum(CountsSampleMUSNorm[,"SHPC30"]),
                             sum(SpikeinMUSNorm$SHPC31)/sum(CountsSampleMUSNorm[,"SHPC31"]),
                             sum(SpikeinMUSNorm$SHPC32)/sum(CountsSampleMUSNorm[,"SHPC32"]),
                             sum(SpikeinMUSNorm$SHPC33)/sum(CountsSampleMUSNorm[,"SHPC33"]),
                             sum(SpikeinMUSNorm$SHPC34)/sum(CountsSampleMUSNorm[,"SHPC34"]),
                             sum(SpikeinMUSNorm$SHPC35)/sum(CountsSampleMUSNorm[,"SHPC35"]),
                             sum(SpikeinMUSNorm$SHPC36)/sum(CountsSampleMUSNorm[,"SHPC36"]))
plot(c(25:36), ratioSpikeNorm, xlab="SHPC", ylab="Spike-in counts / Total counts", main="Normalized")
#text(c(25:36), ratioSpikeNorm, namesSamples)
summary(lm(ratioSpikeNorm ~ c(25:36)))

SpikeinMUS <- merge(txiMusCounts, Spikein, by="row.names")
rownames(SpikeinMUS) <- SpikeinMUS$Row.names

SpikeinMUS <- SpikeinMUS[,-c(1,14,18,19,20)]


ratioSpike <- rbind(sum(SpikeinMUS[,"SHPC25"])/sum(txiMusCounts[,"SHPC25"]),
                    sum(SpikeinMUS[,"SHPC26"])/sum(txiMusCounts[,"SHPC26"]),
                    sum(SpikeinMUS[,"SHPC27"])/sum(txiMusCounts[,"SHPC27"]),
                    sum(SpikeinMUS[,"SHPC28"])/sum(txiMusCounts[,"SHPC28"]),
                    sum(SpikeinMUS[,"SHPC29"])/sum(txiMusCounts[,"SHPC29"]),
                    sum(SpikeinMUS[,"SHPC30"])/sum(txiMusCounts[,"SHPC30"]),
                    sum(SpikeinMUS[,"SHPC31"])/sum(txiMusCounts[,"SHPC31"]),
                    sum(SpikeinMUS[,"SHPC32"])/sum(txiMusCounts[,"SHPC32"]),
                    sum(SpikeinMUS[,"SHPC33"])/sum(txiMusCounts[,"SHPC33"]),
                    sum(SpikeinMUS[,"SHPC34"])/sum(txiMusCounts[,"SHPC34"]),
                    sum(SpikeinMUS[,"SHPC35"])/sum(txiMusCounts[,"SHPC35"]),
                    sum(SpikeinMUS[,"SHPC36"])/sum(txiMusCounts[,"SHPC36"]))
plot(c(25:36), ratioSpike, xlab="SHPC", ylab="Spike-in counts / Total counts", main="Non normalized")
summary(lm(ratioSpike ~ c(25:36)))


## Analyse differences entre echantillons 

SpikeinMUSsample25 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC25"], "SHPC25", "mD")
colnames(SpikeinMUSsample25) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample26 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC26"], "SHPC26", "mD")
colnames(SpikeinMUSsample26) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample27 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC27"], "SHPC27", "mD")
colnames(SpikeinMUSsample27) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample28 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC28"], "SHPC28", "mD")
colnames(SpikeinMUSsample28) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample29 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC29"], "SHPC29", "mD")
colnames(SpikeinMUSsample29) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample30 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC30"], "SHPC30", "mD")
colnames(SpikeinMUSsample30) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample31 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC31"], "SHPC31", "mX")
colnames(SpikeinMUSsample31) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample32 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC32"], "SHPC32", "mX")
colnames(SpikeinMUSsample32) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample33 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC33"], "SHPC33", "mX")
colnames(SpikeinMUSsample33) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample34 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC34"], "SHPC34", "mX")
colnames(SpikeinMUSsample34) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample35 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC35"], "SHPC35", "mX")
colnames(SpikeinMUSsample35) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")
SpikeinMUSsample36 <- data.frame(rownames(SpikeinMUSNorm), SpikeinMUSNorm$concentration, SpikeinMUS[,"SHPC36"], "SHPC36", "mX")
colnames(SpikeinMUSsample36) <- c("Spike_in", "Concentration", "BaseMean", "Echantillon", "Jaw")

SpikeinMUSsample <- rbind(SpikeinMUSsample25, SpikeinMUSsample26, SpikeinMUSsample27, SpikeinMUSsample28, SpikeinMUSsample29, SpikeinMUSsample30, 
                          SpikeinMUSsample31, SpikeinMUSsample32, SpikeinMUSsample33, SpikeinMUSsample34, SpikeinMUSsample35, SpikeinMUSsample36)
#pour mettre a la suite d'un même tableau tous les counts trouves pour chaque spike in de chaque echantillon

par(mfrow=c(1,1))
plot(log(SpikeinMUSsample$Concentration), log(SpikeinMUSsample$BaseMean))

reglin <- lm(SpikeinMUSsample$BaseMean ~ SpikeinMUSsample$Jaw + SpikeinMUSsample$Echantillon : SpikeinMUSsample$Concentration)
summary(reglin)








###### Analyse DESeq : genes non biaises par la souche de mapping ##########################################################################################

#pour voir les genes qui sont differements exprimes selon le mapping, on merge les txi des deux mappings,
#et on effectue un DESeq avec la condition sur le mapping effectue

MergeForDESeq <- merge(txiFVBCounts,txiMusCounts, by="row.names")
rownames(MergeForDESeq) <- MergeForDESeq$Row.names
MergeForDESeq <- MergeForDESeq[-1]
samplesMerge <- data.frame(colnames(MergeForDESeq))
colnames(samplesMerge) <- "mapping_souche"
rownames(samplesMerge) <- samplesMerge$mapping_souche
samplesMerge$mapping_souche <- c("map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_Mus", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB", "map_FVB")

ddsMerge <- DESeqDataSetFromMatrix(MergeForDESeq, colData = samplesMerge, design = ~ mapping_souche)
DESeqMerge <- DESeq(ddsMerge)
DESeqMerge <- DESeqMerge[!is.na(DESeqMerge)]
DESeqMergeResults <- results(DESeqMerge)
table(DESeqMergeResults$padj < 0.05)


SameGenes <- rownames(DESeqMergeResults)[DESeqMergeResults$padj > 0.05]
SameGenes <- SameGenes[!is.na(SameGenes)]


##DESeq avec les genes exprimes de la meme maniere
txiSameCounts <- txiMusCounts[rownames(txiMusCounts)%in%SameGenes,]

ddsSame <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ 1)
#ddsSame <- ddsSame[ rowSums(counts(ddsSame)) > 1, ]
DESeqSame <- DESeq(ddsSame)
DESeqSameResults <- results(DESeqSame)
table(DESeqSameResults$padj < 0.05)

ddsSameStrain <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Strain)
#ddsSameStrain <- ddsSameStrain[ rowSums(counts(ddsSameStrain)) > 1, ]
DESeqSameStrain <- DESeq(ddsSameStrain)
DESeqSameStrainResults <- results(DESeqSameStrain)
table(DESeqSameStrainResults$padj < 0.05)

ddsSameJaw <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Strain+Jaw)
#ddsSameJaw <- ddsSameJaw[ rowSums(counts(ddsSameJaw)) > 1, ]
DESeqSameJaw <- DESeq(ddsSameJaw)
DESeqSameJawResults <- results(DESeqSameJaw)
table(DESeqSameJawResults$padj < 0.05)

ddsSameWeight <- DESeqDataSetFromMatrix(txiSameCounts, colData = samples, design = ~ Strain+Jaw+Weight)
#ddsSameWeight <- ddsSameWeight[ rowSums(counts(ddsSameWeight)) > 1, ]
DESeqSameWeight <- DESeq(ddsSameWeight)
DESeqSameWeightResults <- results(DESeqSameWeight)
table(DESeqSameWeightResults$padj < 0.05)


plotMA(DESeqSameStrainResults, ylim=c(-10,10))



###### Plot biais mapping sur genes DE #######################################################################################################################
DEmapMus <- rownames(txiMusCounts)[DESeqMUSResults$padj < 0.05]
DEmapMus <- DEmapMus[!is.na(DEmapMus)]

DEmapFVB <- rownames(txiFVBCounts)[DESeqFVBResults$padj < 0.05]
DEmapFVB <- DEmapFVB[!is.na(DEmapFVB)]

DEmapSame <- rownames(txiSameCounts)[DESeqSameResults$padj < 0.05]
DEmapSame <- DEmapSame[!is.na(DEmapSame)]

par(mfrow=c(1,2))

DEAll <- unique(rbind(matrix(DEmapMus), matrix(DEmapFVB))) #rbind permet de garder l'ensemble des genes DE, mappes sur C57 OU FVB
MergeDEAll <- Merge[rownames(Merge)%in%DEAll,]
plot(x=MergeDEAll$BaseMeanMus, y=MergeDEAll$BaseMeanFVB, log="xy", xlim=c(0.1,100000), ylim=c(0.1,100000), xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="All DE genes expression")
abline(a=0, b=1, col="red")
summary(lm(log(MergeDEAll$BaseMeanFVB) ~ log(MergeDEAll$BaseMeanMus)))

MergeDEsame <- Merge[rownames(Merge)%in%DEmapSame,] #on garde les genes DE mappes sur C57, seuelement si ils ont la meme empression pour le mapping sur FVB
plot(x=MergeDEsame$BaseMeanMus, y=MergeDEsame$BaseMeanFVB, log="xy", xlim=c(0.1,100000), ylim=c(0.1,100000), xlab="log(Base mean) mapping on C57", ylab="log(Base mean) mapping on FVB", main="DE genes same expression")
#idx <- identify(MergeDE$BaseMeanMus, MergeDE$BaseMeanFVB)
abline(a=0, b=1, col="red")
summary(lm(log(MergeDEsame$BaseMeanFVB) ~ log(MergeDEsame$BaseMeanMus)))

par(mfrow=c(1,1))

# 
# 
# ###### Creation des tables avec genes DE, counts, logfoldchange, padj, pathway ##################################################################################################
# ## counts norm par echantillon, en gardant que les genes non biaises par le mapping
# OConnell <- read.csv(file = "OConnell_list_pathway.txt", header=TRUE, sep="\t")
# OConnellR <- read.csv(file = "OConnell_list_regulatory.txt", header=TRUE, sep="\t")
# Jernvall <- read.csv(file = "Jernvall_list.txt", header=TRUE, sep="\t")
# 
# OConnell <- OConnell[c(1:500),]
# OConnell <- unique(OConnell[,c(1:2)]) #permet de retirer les elements en double, prenant en compte les elements des colonnes 1 et 2
# OConnellR <- unique(OConnellR[,c(1:2)])
# Jernvall <- unique(Jernvall[,c(1:2)])
# 
# OConnell <- na.omit(OConnell)
# OConnellR <- na.omit(OConnellR)
# Jernvall <- na.omit(Jernvall)
# 
# OConnell <- OConnell[!duplicated(OConnell$Gene), ]
# rownames(OConnell) <- OConnell$Gene
# OConnell <- OConnell[-1]
# 
# OConnellR <- OConnellR[!duplicated(OConnellR$Gene), ]
# # pb genes a la fois dans regulator, a la fois dans target
# # pas de genes en double dans Jernvall
# 
# CountsSampleSame <- CountsSampleMUSNorm[rownames(CountsSampleMUSNorm)%in%SameGenes,]
# ddsSameStrain <- ddsSameStrain[ rowSums(counts(ddsSameStrain)) > 1, ]
# 
# 
# padj <- DESeqSameStrainResults$padj
# log2FC <- DESeqSameStrainResults$log2FoldChange
# baseMean <- DESeqSameStrainResults$baseMean #base mean ne change pas d'un dds a l'autre
# genesDESameStrain <- data.frame(CountsSampleSame, baseMean, log2FC, padj)
# genesDESameStrain <- na.omit(genesDESameStrain)
# genesDESameStrain <- merge(genesDESameStrain, OConnell, by="row.names", all = TRUE)
# colnames(genesDESameStrain)[colnames(genesDESameStrain)=="Row.names"] <- "Gene"
# genesDESameStrain <- merge(genesDESameStrain, OConnellR, by="Gene", all = TRUE)
# genesDESameStrain <- merge(genesDESameStrain, Jernvall, by="Gene", all = TRUE)
# genesDESameStrain <- genesDESameStrain[!is.na(genesDESameStrain$baseMean),]
# #write.table(genesDESameStrain, file = "Mapping/Same genes DE Strain.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# padj <- DESeqSameJawResults$padj
# log2FC <- DESeqSameJawResults$log2FoldChange
# baseMean <- DESeqSameJawResults$baseMean #base mean ne change pas d'un dds a l'autre
# genesDESameJaw <- data.frame(CountsSampleSame, baseMean, log2FC, padj)
# genesDESameJaw <- na.omit(genesDESameJaw)
# genesDESameJaw <- merge(genesDESameJaw, OConnell, by="row.names", all = TRUE)
# colnames(genesDESameJaw)[colnames(genesDESameJaw)=="Row.names"] <- "Gene"
# genesDESameJaw <- merge(genesDESameJaw, OConnellR, by="Gene", all = TRUE)
# genesDESameJaw <- merge(genesDESameJaw, Jernvall, by="Gene", all = TRUE)
# genesDESameJaw <- genesDESameJaw[!is.na(genesDESameJaw$baseMean),]
# #write.table(genesDESameJaw, file = "Mapping/Same genes DE Jaw.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# padj <- DESeqSameWeightResults$padj
# log2FC <- DESeqSameWeightResults$log2FoldChange
# baseMean <- DESeqSameWeightResults$baseMean #base mean ne change pas d'un dds a l'autre
# genesDESameWeight <- data.frame(CountsSampleSame, baseMean, log2FC, padj)
# genesDESameWeight <- na.omit(genesDESameWeight)
# genesDESameWeight <- merge(genesDESameWeight, OConnell, by="row.names", all = TRUE)
# colnames(genesDESameWeight)[colnames(genesDESameWeight)=="Row.names"] <- "Gene"
# genesDESameWeight <- merge(genesDESameWeight, OConnellR, by="Gene", all = TRUE)
# genesDESameWeight <- merge(genesDESameWeight, Jernvall, by="Gene", all = TRUE)
# genesDESameWeight <- genesDESameWeight[!is.na(genesDESameWeight$baseMean),]
# #write.table(genesDESameWeight, file = "Mapping/Same genes DE Weight.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# 
# 
# ## Test de probabilite des genes appartenant a une voie de la liste OConnell
# pOConnell=222/15263
# 
# 
# chisq.test(x=c(54,3068-54),p=c(pOConnell,1-pOConnell))
# chisq.test(x=c(45,1606-45),p=c(pOConnell,1-pOConnell))
# chisq.test(x=c(16,622-16),p=c(pOConnell,1-pOConnell))
# 
# 
# 
# 
# ###### Overlap des genes DE entre conditions #############################################################################################################################################
# genesDEoverlapStrainJaw <- merge(genesDESameStrain[genesDESameStrain$padj<0.05,], genesDESameJaw[genesDESameJaw$padj < 0.05,], by="row.names")
# dim(genesDEoverlapStrainJaw)
# genesDEoverlapJawWeight <- merge(genesDESameJaw[genesDESameJaw$padj<0.05,], genesDESameWeight[genesDESameWeight$padj < 0.05,], by="row.names")
# dim(genesDEoverlapJawWeight)
# genesDEoverlapWeightStrain <- merge(genesDESameWeight[genesDESameWeight$padj<0.05,], genesDESameStrain[genesDESameStrain$padj < 0.05,], by="row.names")
# dim(genesDEoverlapWeightStrain)
# 
# genesDEoverlapAll <- merge(genesDEoverlapStrainJaw, genesDESameWeight[genesDESameWeight$padj<0.05,], by="row.names")
# dim(genesDEoverlapAll)
# 
# genesDEoverlap <- merge(genesDEoverlapStrainJaw$Row.names, genesDEoverlapJawWeight$Row.names, by="row.names", all=TRUE)
# genesDEoverlap <- merge(genesDEoverlap[,-1], genesDEoverlapWeightStrain$Row.names, by="row.names", all=TRUE)[,-1]
# colnames(genesDEoverlap) <- c("Strain & Jaw", "Jaw & Weight", "Weight & Strain")
# #write.table(genesDEoverlap, file = "Mapping/Same genes DE Overlap.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# table(genesDEoverlapStrainJaw$log2FC.x > 0, genesDEoverlapStrainJaw$log2FC.y < 0)
# 
# 
# ## Test de probabilite de genes communs
# pt=3068/15263
# pj=1606/15263
# pw=622/15263
# 
# chisq.test(x=c(432,15263-432),p=c(pt*pj,1-pt*pj))
# chisq.test(x=c(142,15263-142),p=c(pj*pw,1-pj*pw))
# chisq.test(x=c(189,15263-189),p=c(pw*pt,1-pw*pt))
# 
# 
# 
# 
# 
# ###### Liste de genes pour Gorilla et String ###############################################################################################################
# sumSame <- data.frame(rowSums (txiSameCounts, na.rm = FALSE, dims = 1))
# 
# GOSameAll <- rownames(txiSameCounts)[sumSame > 0]
# #write.table(GoSameAll, file = "Mapping/Gorilla same genes DE All.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 
# 
# GoSameDE <- rownames(txiSameCounts)[DESeqSameResults$padj < 0.05 & (DESeqSameResults$log2FoldChange < -1 | DESeqSameResults$log2FoldChange > 1) ]
# GoSameDE <- na.omit(GoSameDE)
# #write.table(GoSameDE, file = "Mapping/Gorilla same genes DE 1 lfc.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 
# GoSameDEStrain <- rownames(txiSameCounts)[DESeqSameStrainResults$padj < 0.05 & (DESeqSameResults$log2FoldChange < -1 | DESeqSameResults$log2FoldChange > 1)]
# GoSameDEStrain <- na.omit(GoSameDEStrain)
# #write.table(GoSameDEStrain, file = "Mapping/Gorilla same genes DE Strain lfc.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 
# GoSameDEJaw <- rownames(txiSameCounts)[DESeqSameJawResults$padj < 0.05 & (DESeqSameResults$log2FoldChange < -1 | DESeqSameResults$log2FoldChange > 1)]
# GoSameDEJaw <- na.omit(GoSameDEJaw)
# #write.table(GoSameDEJaw, file = "Mapping/Gorilla same genes DE Jaw lfc.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 
# GoSameDEWeight <- rownames(txiSameCounts)[DESeqSameWeightResults$padj < 0.05 & (DESeqSameResults$log2FoldChange < -1 | DESeqSameResults$log2FoldChange > 1)]
# GoSameDEWeight <- na.omit(GoSameDEWeight)
# #write.table(GoSameDEWeight, file = "Mapping/Gorilla same genes DE Weight lfc.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 
# 
# 
# 
# 
# ###### Analyse PCA des genes non biaises par le mapping ####################################################################################################
# 
# plotPCA.mystyle2 <-  function(object, colorgroup="condition", shapegroup="condition2",ntop, returnData=FALSE, pcs = c(1,2),scale)
# {  stopifnot(length(pcs) == 2)    ## added this to check number of PCs 
#   # calculate the variance for each gene
#   rv <- rowVars(assay(object))
#   # select the ntop genes by variance
#   select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
#   # perform a PCA on the data in assay(x) for the selected genes
#   ### ICI difference avec DESeq2 : scale =T
#   pca <- prcomp(t(assay(object)[select,]),scale. = scale)
#   # the contribution to the total variance for each component
#   percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
#   if (!all(colorgroup %in% names(colData(object)))) {
#     stop("the argument 'colorgroup' should specify columns of colData(dds)")
#   }
#   colorgroup.df <- as.data.frame(colData(object)[, colorgroup, drop=FALSE])
#   shapegroup.df <- as.data.frame(colData(object)[, shapegroup, drop=FALSE])
#   
#   # add the colorgroup factors together to create a new grouping factor
#   Weight <- if (length(colorgroup) > 1) {
#     factor(apply( colorgroup.df, 1, paste, collapse=" : "))
#   } else {
#     colData(object)[[colorgroup]]
#   }
#   Jaw <- if (length(shapegroup) > 1) {
#     factor(apply(shapegroup.df, 1, paste, collapse=" : "))
#   } else {
#     colData(object)[[shapegroup]]
#   }
#   # assembly the data for the plot
#   ### Here we just use the pcs object passed by the end user
#   d <- data.frame(PCa=pca$x[,pcs[1]], PCb=pca$x[,pcs[2]], Color=Weight, Shape=Jaw, colorgroup.df, name=colnames(object))
#   plot = ggplot(data = d, aes_string(x = "PCa", y = "PCb", color = "Weight", shape="Jaw")) + 
#     geom_point(size = 5) + xlab(paste0("PC", pcs[1],": ", round(percentVar[pcs[1]] * 100), "% variance")) + 
#     ylab(paste0("PC", pcs[2],": ", round(percentVar[pcs[2]] * 100), "% variance")) + coord_fixed() +
#     theme(axis.text.x= element_text(size=25),
#           axis.text.y= element_text(size=25),
#           axis.title.x= element_text(size=25, margin = margin(t = 20, r = 0, b = 0, l = 0)), 
#           axis.title.y= element_text(size=25),
#           legend.text = element_text(size=20),
#           legend.title = element_text(size=25),
#           legend.key.size = unit(3, 'lines'))#+ geom_text(aes_string(x = "PCa", y = "PCb", label = "name"), color = "black") + theme_bw()
#   print(plot)
#   if (returnData) {
#     attr(d, "percentVar") <- percentVar
#     return(list(d, plot))
#   }
# }
# 
# ## Plot PCA
# library(ggplot2)
# rldSame <- rlogTransformation(DESeqSame[sumSame>0,], blind=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,2), scale=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,3), scale=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,4), scale=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,5), scale=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,6), scale=TRUE)
# plotPCA.mystyle2(rldSame, colorgroup=c("Weight"), shapegroup=c("Jaw"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,7), scale=TRUE)
# 
# 
# 
# ## DUDI PCA
# library(ade4)
# countrldSame=assay(rldSame)
# dudiSamePCA=dudi.pca(t(countrldSame), scale=T, scannf=F, nf=5)
# 
# dudilambda=dudiSamePCA$eig/sum(dudiSamePCA$eig)
# 
# 
# 
# ###### Noms des genes qui jouent le plus sur les axes PCA ####################################################################
# library("factoextra")
# 
# #Coloration en fonction du cos2 (qualite de representation). Les individus similaires sont groupes ensemble.
# fviz_eig(dudiSamePCA) # Pourcentage de variance expliquee par chaque axe de PCA
# fviz_pca_ind(dudiSamePCA,
#              col.ind = "cos2", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     
# ) 
# 
# #Graphique des variables. Coloration en fonction de la contribution des variables. 
# #Les variables correlees positivement sont du meme cote du graphique. Les variables correlees negativement sont sur des cotes opposes du graphique.
# #fviz_pca_var(dudiSamePCA,
# #col.var = "contrib", 
# #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
# #repel = TRUE)    
# 
# ### NE FONCTIONNE PAS ??
# 
# 
# res.var <- get_pca_var(dudiSamePCA)
# #res.var$coord          # Coordonnees
# res.var$contrib        # Contributions aux axes
# #res.var$cos2           # Qualite de representation 
# 
# contribPCA <- res.var$contrib
# contribPCA <- merge(contribPCA, OConnell, by="row.names", all = TRUE)
# colnames(contribPCA)[colnames(contribPCA)=="Row.names"] <- "Gene"
# contribPCA <- merge(contribPCA, OConnellR, by="Gene", all = TRUE)
# contribPCA <- merge(contribPCA, Jernvall, by="Gene", all = TRUE)
# contribPCA <- contribPCA[!is.na(contribPCA$Dim.1),]
# 
# #write.table(contribPCA, file = "Mapping/Participation des genes au axes PCA.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 
# res.ind <- get_pca_ind(dudiSamePCA)
# #res.ind$coord
# res.ind$contrib
# #res.ind$cos2
# 
# 
# 
# ###### Interactions Strain:Jaw #############################################################################################################################################
# 
# ddsSameInter <- DESeqDataSetFromMatrix(txiSameCounts[sumSame>1,], colData = samples, design = ~ Strain + Strain:Jaw)
# DESeqSameInter <- DESeq(ddsSameInter)
# DESeqSameInterResults <- results(DESeqSameInter)
# mcols(DESeqSameInterResults)
# 
# resultsNames(DESeqSameInter)
# 
# InterSameStrain <- results(DESeqSameInter, name = 'Strain_DUHi_vs_FVB')
# InterSameStrainJaw <- results(DESeqSameInter, name = "StrainDUHi.JawmX" )
# 
# InterSamelfc <- InterSameStrain$log2FoldChange[InterSameStrain$padj < 0.05]
# InterSamelfc <- InterSamelfc[!is.na(InterSamelfc)]
# InterSamenames <- rownames(InterSameStrain)[InterSameStrain$padj < 0.05]
# InterSamenames <- InterSamenames[!is.na(InterSamenames)]
# InterSameStrainTable <- data.frame(InterSamenames,InterSamelfc)
# 
# InterSamelfc <- InterSameStrainJaw$log2FoldChange[InterSameStrainJaw$padj < 0.05]
# InterSamelfc <- InterSamelfc[!is.na(InterSamelfc)]
# InterSamenames <- rownames(InterSameStrainJaw)[InterSameStrainJaw$padj < 0.05]
# InterSamenames <- InterSamenames[!is.na(InterSamenames)]
# InterSamepadj <- InterSameStrainJaw$padj[InterSameStrainJaw$padj < 0.05]
# InterSamepadj <- InterSamepadj[!is.na(InterSamepadj)]
# InterSameStrainJawTable <- data.frame(InterSamenames,InterSamelfc,InterSamepadj)
# rownames(InterSameStrainJawTable) <- InterSameStrainJawTable$InterSamenames
# InterSameStrainJawTable <- InterSameStrainJawTable[-1]
# InterSameStrainJawTable <- merge(InterSameStrainJawTable, OConnell, by="row.names", all = TRUE)
# colnames(InterSameStrainJawTable)[colnames(InterSameStrainJawTable)=="Row.names"] <- "Gene"
# InterSameStrainJawTable <- merge(InterSameStrainJawTable, OConnellR, by="Gene", all = TRUE)
# InterSameStrainJawTable <- merge(InterSameStrainJawTable, Jernvall, by="Gene", all = TRUE)
# InterSameStrainJawTable <- InterSameStrainJawTable[!is.na(InterSameStrainJawTable$InterSamelfc),]
# write.table(InterSameStrainJawTable, file = "Mapping/Same genes DE Inter Strain Jaw.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# 
# 
# ## Comprendre les interactions
# cn=counts(DESeqSameInter,norm=T)
# head(cn)
# tapply(cn[1,],paste(samples$Strain,samples$Jaw),mean)
# results(DESeqSameInter)[1,]
# summary(cn[1,])
# resultsNames(DESeqSameInter)
# results(DESeqSameInter,name="Intercept")[1,]
# #DESeqSameInterResults <- results(DESeqSameInter, contrast=c("FVB","mD", "mX"))
# #test=ref*2^(-log2foldchange)
# 
# 
# ## PCA sur genes avec interaction
# 
# 
# rldInter <- rldSame[rownames(rldSame)%in%InterSameStrainJawTable$Gene,]
# plotPCA.mystyle2(rldInter, colorgroup=c("Weight"), shapegroup=c("Strain"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,2), scale=TRUE)
# plotPCA.mystyle2(rldInter, colorgroup=c("Weight"), shapegroup=c("Strain"), ntop=nrow(rldSame), returnData=TRUE, pcs = c(1,3), scale=TRUE)
# 
# 
# countrldInter=assay(rldInter)
# dudiSamePCAInter=dudi.pca(t(countrldInter), scale=T, scannf=F, nf=5)
# fviz_eig(dudiSamePCAInter)
# 
# 
# 
# 
# 
# 
# ###### Histogramme genes avec pathway #################################################################################################################################
# 
# 
# samplesplot <- samples
# samplesplot$StrainJaw <- paste(samplesplot$Jaw,samplesplot$Strain)
# samplesplot$StrainJaw <- factor(samplesplot$StrainJaw, levels = c("mD FVB", "mD DUHi", "mX FVB", "mX DUHi"))
# 
# CountsPathway <- merge(CountsSampleSame, OConnell, by="row.names", all = TRUE)
# colnames(CountsPathway)[colnames(CountsPathway)=="Row.names"] <- "Gene"
# CountsPathway <- merge(CountsPathway, OConnellR, by="Gene", all = TRUE)
# CountsPathway <- merge(CountsPathway, Jernvall, by="Gene", all = TRUE)
# CountsPathway <- CountsPathway[!is.na(CountsPathway$SHPC25),]
# CountsPathway <- CountsPathway[!(is.na(CountsPathway$Pathway) & is.na(CountsPathway$Regulatory) & is.na(CountsPathway$Jernvall )),]
# #On supprime les genes qui ne sont dans aucune des listes
# rownames(CountsPathway) <- CountsPathway$Gene
# CountsPathway <- CountsPathway[-1]
# 
# 
# plotdvt = function(name){samplesplot$counts <- t(CountsPathway[name,c(1:12)])
# ggplot(samplesplot,aes(x = StrainJaw, y = counts, fill = Weight))+
#   geom_bar(stat="identity",position='dodge') + theme_bw() + ggtitle(name) + xlab("")
# 
# }
# 
# ## Pathway
# plotgenes = rownames(CountsPathway)
# 
# #pdf("Pathway_Hist.pdf")
# #lapply(plotgenes, plotdvt)
# #dev.off()
# 
# 
# 
# ###### Plot genes avec interaction ####################################################################################
# 
# 
# df <- samplesplot
# 
# plotInterAll = function(name){df$counts <- as.numeric(counts(DESeqSameInter[name,],norm=T))
# pv=format.pval(InterSameStrainJaw[name,"padj"])
# lfc = InterSameStrainJaw[name,"log2FoldChange"]
# d=ggplot(df, aes(x=Jaw, y=counts, color=Strain, shape=Weight)) + geom_point() + theme_bw() + ggtitle(paste(name,pv,lfc))
# d + aes(x=Jaw, color=Strain, group=Strain) + stat_summary(fun.y = mean, geom="line") +
#   geom_point(size = 3) +
#   theme(axis.text.x= element_text(size=14, colour = "black", margin = margin(t = 7, r = 0, b = 0, l = 0)),
#         axis.text.y= element_text(size=18, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
#         axis.title.x= element_text(size=20), 
#         axis.title.y= element_text(size=20),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=20),
#         plot.title = element_text(size=20))
# }
# 
# #pdf("Inter_Plot.pdf")
# #lapply(InterSameStrainJawTable$Row.names, plotInterAll)
# #dev.off()
# 
# 
# 
# plotdvtAll = function(name){samplesplot$counts <- CountsSampleSame[name,]
# ggplot(samplesplot,aes(x = StrainJaw, y = counts, fill = Weight))+
#   geom_bar(stat="identity", position='dodge') + theme_bw() + ggtitle(name) + xlab("") +
#   theme(axis.text.x= element_text(size=14, colour = "black", margin = margin(t = 7, r = 0, b = 0, l = 0)),
#         axis.text.y= element_text(size=18, colour = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
#         axis.title.x= element_text(size=20), 
#         axis.title.y= element_text(size=20),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=20),
#         plot.title = element_text(size=20)
#         )
# }
# 
# #pdf("Inter_Hist.pdf")
# #lapply(InterSameStrainJawTable$Row.names, plotdvtAll)
# #dev.off()
# 
# 
# 
# 
# 
# ###### Acide rhetinoique ################################################################################################
# 
# AcReth <- c("Aldh1a1", "Aldh1a2", "Aox3", "Cyp26a1", "Cyp26c1", "Cyp27a1", "Crabp1", "Rarb", "Dhrs3", "Nr2f1", "Sox10")
# 
# #pdf("AcReth_Hist.pdf")
# #lapply(AcReth, plotdvtAll)
# #dev.off()
# 
# #pdf("AcReth_Plot.pdf")
# #lapply(AcReth, plotInterAll)
# #dev.off()
# 
# 
# 
# 
# 
# ###### Diagramme de Venn ###########################################################################################################
# 
# install.packages("VennDiagram")
# library("VennDiagram")
# 
# draw.triple.venn(area1 = 3068, area2 = 1606, area3 = 622, n12 = 432, n23 = 142, n13 = 189, n123 = 13, 
#                  category = c("Strain", "Jaw", "Weight"), fill=c("green","red","yellow"),
#                  cex = 2, cat.cex = 2, scaled=TRUE, lty = "blank")
# 
# 
# 
# AcReth <- c("Aldh1a1", "Aldh1a2", "Aox3", "Cyp26a1", "Cyp26c1", "Cyp27a1", "Crabp1", 
#             "Rarb", "Dhrs3", "Nr2f1", "Sox10")
# 
# name = "Cyp26a1"
# png(paste(name,"_gros.png"))
# plotdvtAll(name)
# dev.off()
# 






