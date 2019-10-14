library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggfortify)

setwd("C:/Users/awilczyn/OneDrive - University of Glasgow/Embo workshop/counts_files_edgeR")
#read in targets
targets <- readTargets("C:/Users/awilczyn/OneDrive - University of Glasgow/Embo workshop/counts_files_edgeR/targets_KD.txt", sep="\t")
targets
samples<- gsub("htseq-count_on_", "",targets$files)
samples<- gsub(".tabular", "",samples)
samples

#read in count tables (skip and comment remove comments from analysis)
d <- readDGE(targets, skip=1, comment.char = "!",group=gsub("_.{1}_","",samples), labels=samples)
dim(d$counts)
View(d$counts)

#see new sample table
d$samples
head(d$counts)
summary(d$counts)
dim(d)

#discard all genes where at least 2 samples don't have at least 2 reads per million
keep <- rowSums(cpm(d) > 4) >= 2
d <- d[keep,]
dim(d)
#check how many discarded
table(keep)

#Having filtered, reset the library sizes:
d$samples$lib.size <- colSums(d$counts)

#TMM normalization is applied to this dataset to account for compositional difference between
#the libraries.
d <- calcNormFactors(d)
d$samples

#plot in which distances between samples correspond to leading biological coeffcient of variation (BCV) between those samples
x<-d$samples$batch
x<- replace(x,x=="1","red")
x<- replace(x,x=="2","blue")
x<- replace(x,x=="4","green")
x<- replace(x,x=="5","purple")
##MDS plot without labels###
plotMDS(d,col=x,pch=c(1,16,2, 17,5, 18)[as.numeric(d$samples$group)])
legend("topleft", legend=c("1", "2","3","4"),
       col=c("red", "blue","green","purple"), pch=20, cex=2)

##MDS plot with labels###
plotMDS(d,col=x, labels =rownames(d$samples))


##PCA plot
## Transform count data
rld <- rlog(d$counts)
## Perform PCA analysis and make plot
autoplot(prcomp(t(rld)),data= t(rbind(rld, t(d$samples))),label=T, colour='batch')+theme_bw()+ylim(-0.3,0.3)+
  xlim(-0.3,0.5)

#create design matrix for additive linear model with batch as blocking factor
batch <- factor(d$samples$batch)
condition <- factor(d$samples$group)


#create design matrix for additive linear model with batch as blocking factor
design <- model.matrix(~0+condition+batch)
design

#check correlation of each batch
logFC <- predFC(d,design,prior.count = 1,dispersion = 0.05)
cor(logFC[,1:9])


#estimate the overall dispersion for the dataset, to get an idea of the overall level of
#biological variability:
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)

#estimate gene-wise dispersion estimates, allowing a possible trend with averge count size
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
plotBCV(d)

#determine differentially expressed genes. Fit genewise glms:
fit <- glmFit(d, design)


#Conduct likelihood ratio tests for total RNA KD vs ctrl differences and show the top genes:
colnames(design)
lrt <- glmLRT(fit, contrast=c(0,0,-1,1,0,0,0,0,0))
topTags(lrt)


#counts-per-million in individual samples for the top genes:
o <- order(lrt$table$PValue)
cpm(d)[o[1:10],]

#total number of differentially expressed genes at 5% FDR is given by:
summary(de <- decideTestsDGE(lrt))

## MA plots #####
#plot log-fold change against log-counts per million, with DE genes highlighted:
detags <- rownames(d)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#save file
sub_tab <- topTags(lrt, n = Inf)   
write.table(sub_tab, file="Sub_KD_vs_ctrl.txt",sep="\t")
cpm<-cpm(d)
write.table(cpm,file="cpm_for_workshop.txt",sep="\t")

### same for polysomes####
#Conduct likelihood ratio tests for total RNA KD vs ctrl differences and show the top genes:
colnames(design)
lrt <- glmLRT(fit, contrast=c(-1,1,0,0,0,0,0,0,0))
topTags(lrt)

#total number of differentially expressed genes at 5% FDR is given by:
summary(de <- decideTestsDGE(lrt))

## MA plots #####
#plot log-fold change against log-counts per million, with DE genes highlighted:
detags <- rownames(d)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

#save file
poly_tab <- topTags(lrt, n = Inf)   
write.table(poly_tab, file="Poly_KD_vs_ctrl.txt",sep="\t")

###scatterplots for sub vs poly###

#combine two datasets#
poly_sub<- merge(sub_tab$table, poly_tab$table, by=0, suffixes=c(".sub",".poly"))
head(poly_sub)

dot<-ggplot() + 
  geom_point(data=poly_sub,aes(x=logFC.sub, y=logFC.poly),colour="gray",size=2)+
  geom_point(data=poly_sub[poly_sub$FDR.sub< 0.2 & poly_sub$FDR.poly< 0.2,],aes(x=logFC.sub, y=logFC.poly),
             colour="pink",size=3)

dot+
  theme_bw()+
  theme(legend.text=element_text(size=12),axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20,face="bold"),
        panel.grid.major = element_line(color = "gray",linetype = "dotted"), legend.title=element_blank())+
  xlim(-2,2)+ylim(-2,2)+ylab("POLYSOME")+xlab("SUBPOLYSOME")+labs(title="significant genes")+
  geom_abline(size=1, linetype="dashed")

###calculate difference between logFCs##
poly_sub$delta<- poly_sub$logFC.poly - poly_sub$logFC.sub

### mark genes that show large differences in opposite directions ###

dot<-ggplot() + 
  geom_point(data=poly_sub,aes(x=logFC.sub, y=logFC.poly),colour="gray",size=2)+
  geom_point(data=poly_sub[poly_sub$FDR.sub< 0.2 & poly_sub$FDR.poly< 0.2,],aes(x=logFC.sub, y=logFC.poly),colour="pink",size=3)+
  geom_point(data=poly_sub[poly_sub$FDR.sub< 0.2 & poly_sub$FDR.poly< 0.2 & abs(poly_sub$delta)>0.5,],
             aes(x=logFC.sub, y=logFC.poly),colour="deeppink",size=3)

dot+
  theme_bw()+
  theme(legend.text=element_text(size=12),axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20,face="bold"),
        panel.grid.major = element_line(color = "gray",linetype = "dotted"), legend.title=element_blank())+
  xlim(-2,2)+ylim(-2,2)+ylab("POLYSOME")+xlab("SUBPOLYSOME")+labs(title="differential genes")+
  geom_abline(size=1, linetype="dashed")

###make table with translationally affected genes##
edger_trans<- poly_sub[poly_sub$FDR.sub< 0.2 & poly_sub$FDR.poly< 0.2 & abs(poly_sub$delta)>0.5,]

########### POLY vs Total with KD vs ctrl ######

#Conduct likelihood ratio tests for Poly vs total, RNA KD vs ctrl differences and show the top genes:
colnames(design)
lrt <- glmLRT(fit, contrast=c(-1,1,0,0,1,-1,0,0,0))
topTags(lrt)
