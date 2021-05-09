# Required packages
library(S4Vectors)
library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(caret)
library(ROCR)
library(pROC)
library(GenomicRanges)
library(dplyr)
library(MLmetrics)
source(file = '/home/yuxin/home/yuxin/Nanopore/Script/One_Hot.R')

# The Tombo and basecall error features at each position are stored in grange object
# Note that not all features fore site on the genome can be found, some positions may not be covered by nanopore reads.

Coverage_Grange <- readRDS("/home/yuxin/home/yuxin/Nanopore/DATA/Coverage_Grange.rds") # Tombo coverage level
Fraction_Grange <- readRDS("/home/yuxin/home/yuxin/Nanopore/DATA/Fraction_Grange.rds") # Tombo modification fraction
Error_Grange <- readRDS("/home/yuxin/home/yuxin/Nanopore/DATA/Error_Grange.rds") # Base call error

# Extract feature for sites of interest
# sample: a grange object for target site (single base resolution)
# window: an integer, specify the length of flanking windows along the target site
FeatureExtraction <- function(sample,window){
  Fraction <- Fraction_Grange[Fraction_Grange%over%(sample+window)]
  Coverage <- Coverage_Grange[Coverage_Grange%over%sample]
  Error <- Error_Grange[Error_Grange%over%sample]
  Index <- intersect(intersect(which(sample%over%Fraction),
                               which(sample%over%Coverage)),
                     which(sample%over%Error))
  sample <- sample[Index]
  
  Seq <- as.character(DNAStringSet(Views(Hsapiens,sample+20)))
  mcols(sample) <- NULL
  OneHot <- ONE_HOT(Seq)
  
  Features <- matrix(NA,nrow = length(sample),ncol = (2*window+11))
  for (i in 1:nrow(Features)){
    Features[i,1] <- mcols(Fraction[Fraction%over%sample[i]])[,1] 
    if(length(mcols(Fraction[Fraction%over%(sample[i]+window)])[,1])==2*window+1){
      Features[i,2:(2*window+2)] <- mcols(Fraction[Fraction%over%(sample[i]+window)])[,1]
    }
    Features[i,(2*window+3)] <- mean(mcols(Coverage[Coverage%over%sample[i]])[,1])
    Features[i,(2*window+4)] <- mean(width(Coverage[Coverage%over%sample[i]]))
    Features[i,(2*window+5):(2*window+11)] <-colMeans(as.matrix(mcols(Error[Error%over%sample[i]])))
  }
  Features <- cbind(Features,OneHot)
  return(Features)
}

# Modification site path: for negative site, just replace 'Positive' with 'Negative'
Positive <- readRDS('/data/kunqidir/dpWHISTLE/hm1A/FullPostive.rds') # m1A
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/hAm/hAm_more/FullPostive.rds') # Am
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/pseuoduridine/hPsi/FullPostive.rds') # Pseudouridine
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/hCm/hCm_more/FullPostive.rds') # Cm  
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/hGm/hGm_more/FullPostive.rds') # Gm
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/hTm/hTm_more/FullPostive.rds') # Tm
Positive <- readRDS('/data/kunqidir/dpWHISTLE/hm7G/FullPostive.rds') # m7G
Positive <- readRDS('/data/yuxin/Nanopore/mod_site/hm5U/hm5U/FullPostive.rds') # m5U

# m6A modification site
# Positive site
Mix <- '/home/yuxin/home/yuxin/Nanopore/DATA/modSites/m6ASites/Graphic/Positive.rds'
High_Positive <- '/data/kunqidir/Nanopore/high_sites.rds'
Medium_Positive <- '/data/kunqidir/Nanopore/medium_sites.rds'
Low_Positive <- '/data/kunqidir/Nanopore/low_sites.rds'
Verylow_Positive <- '/data/kunqidir/Nanopore/Verylow_sites.rds'

# Extract negative site
Transcripts <- unique((transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)))
Transcripts <- Transcripts[Transcripts%over%readRDS(Mix)]
DRACH_motif <- m6ALogisticModel::sample_sequence("DRACH",Transcripts,Hsapiens)
Negative_total <-DRACH_motif[which(vcountPattern("NNN",DNAStringSet(Views(Hsapiens,DRACH_motif)))==0)]
Negative_total <- Negative_total[-which(Negative_total%over%readRDS(Mix))]
has <- intersect(intersect(which((Negative_total-2)%over%Coverage_Grange),
                           which((Negative_total-2)%over%Fraction_Grange)),
                 which((Negative_total-2)%over%Error_Grange))
Neg_Has <- Negative_total[has]

Positive <- readRDS(High_Positive)
pos_has <- intersect(intersect(which(Positive%over%Coverage_Grange),
                               which(Positive%over%Fraction_Grange)),
                     which(Positive%over%Error_Grange))
Positive <- Positive[pos_has]

Index <- createDataPartition(y = seqnames(Neg_Has), p = (length(Positive)/length(Neg_Has)), list = FALSE)
Negative  <- Neg_Has[Index]-2

# Training
pos <- FeatureExtraction(Positive,window)
neg <- FeatureExtraction(Negative,window)
Whole <- as.data.frame(rbind(pos,neg))
Whole$Class <- factor(c(rep('Positive',nrow(pos)),rep('Negative',nrow(neg))), 
                      levels = c('Positive','Negative'))

Whole <- na.omit(Whole)
trainInd <- createDataPartition(y = Whole$Class, p = 0.8, list = FALSE)
Training <- Whole[trainInd,]
Testing <- Whole[-trainInd,]

fitControl <- trainControl(method = "cv",
                           number = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)

Model <- train(Class~., data = Training,
               method= 'svmRadial',
               preProc = c('center','scale'),
               trControl = fitControl,
               metric = 'ROC')

# Evaluation
Testing$predicted_class <- predict(Model,Testing,type = 'prob')
perf <- prediction(Testing$predicted_class$Positive,Testing$Class)
perf_auc <- performance(perf,'auc')
auc <- attr(perf_auc,'y.values')[[1]][1]
prauc <- PRAUC(Testing$predicted_class$Positive,Testing$Class)
roc <- max(Model$results$ROC)
