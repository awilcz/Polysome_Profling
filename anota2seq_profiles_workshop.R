library(anota2seq)
library(ggfortify)
library(ggplot2)


setwd("C:/Users/awilczyn/OneDrive - University of Glasgow/Embo workshop/counts_files_edgeR")
filelist = list.files(pattern = ".*TOT|POLY.*.tabular")
filelist


##read in files and make sample names###
datalist = lapply(filelist, function(x) read.table(x, sep="\t",header=F)) 
samples<- gsub("htseq-count_on_","" ,filelist)
samples<- gsub(".tabular","" ,samples)
samples<- c("gene",samples)
samples

datalist<-lapply(datalist, function(x) as.data.frame(x))

#create counts matrix### 
countframe<- Reduce(function(...) merge(..., all=T, by=1), datalist)
colnames(countframe)<- samples
rownames(countframe)<- countframe$gene
countframe$gene<- NULL
head(countframe)

###we need to check for outliers, since the model assumes there are no significant outliers###
###this can be done by PCA and is also tested for within anota itself####

## Transform count data
rld <- rlog(as.matrix(countframe))
## Perform PCA analysis and make plot
autoplot(prcomp(t(rld)),data= t(rld),label=T, colour=as.factor(c("1","1","2","2","3","3","4","4","1","1","2","2","3","3","4","4")))+theme_bw()

##perform analysis###
ads <- anota2seqDataSetFromMatrix(dataT = countframe[,c(9:16)],
                                  data = countframe[,c(1:8)],
                                  phenoVec = c("ctrl","KD","ctrl","KD","ctrl","KD","ctrl","KD"),
                                  batchVec = c("1","1","2","2","3","3","4","4"),
                                  dataType = "RNAseq",
                                  normalize = TRUE)
ads_def <- anota2seqRun(ads)

###make differential expression scatterplot
anota2seqPlotFC(ads_def, selContrast = 1, plotToFile = FALSE)
##and save it###
pdf('KD_anota_defaults.pdf')

###get gene table###
ads_def_genes<- anota2seqGetOutput(ads_def,output="singleDf",selContrast=1)
write.table(ads_def_genes, "anota2seq_output_table.txt", sep="\t")

##get only translational genes###
trans_ads<- ads_def_genes[ads_def_genes$singleRegMode=="translation",]

### plot edgeR results vs anota2seq results##
edger_anota_merge<- merge(edger_trans,ads_def_genes[ads_def_genes$translation.apvRvmPAdj<0.2,], by=1, all=F)

dot<-ggplot() + 
  geom_point(data=edger_anota_merge,aes(x=delta, y=translation.apvEff),colour="gray",size=2)+
  geom_point(data=edger_anota_merge[edger_anota_merge$singleRegMode=="translation",],aes(x=delta, y=translation.apvEff),
             colour="gold2",size=3)
  
dot+
  theme_bw()+
  theme(legend.text=element_text(size=12),axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20,face="bold"),
        panel.grid.major = element_line(color = "gray",linetype = "dotted"), legend.title=element_blank())+
  xlim(-2,2)+ylim(-2,2)+ylab("translation effect anota")+xlab("polysom-subpolysome erdgeR")+labs(title="comparison of the two approaches")+
  geom_abline(size=1, linetype="dashed")




###########################################################
#### let's see how it correlates with mass spec results ###

silac<- read.table("SILAC_data_for_KD.txt",sep="\t",header = T)

silac_anota_merge<- merge(silac,ads_def_genes[ads_def_genes$translation.apvRvmPAdj<0.2,], by=1, all=F)

dot<-ggplot() + 
  geom_point(data=silac_anota_merge,aes(x=average_FC, y=translation.apvEff),colour="gray",size=2)+
  geom_point(data=silac_anota_merge[silac_anota_merge$singleRegMode=="translation",],aes(x=average_FC, y=translation.apvEff),
             colour="gold2",size=3)

dot+
  theme_bw()+
  theme(legend.text=element_text(size=12),axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20,face="bold"),
        panel.grid.major = element_line(color = "gray",linetype = "dotted"), legend.title=element_blank())+
  xlim(-2,2)+ylim(-2,2)+ylab("translation effect anota")+xlab("SILAC protein amounts")+labs(title="comparison of the two approaches")+
  geom_abline(size=1, linetype="dashed")


silac_edger_merge<- merge(silac,edger_trans, by=1, all=F)

dot<-ggplot() + 
  geom_point(data=silac_edger_merge,aes(x=average_FC, y=delta),colour="gray",size=2)

dot+
  theme_bw()+
  theme(legend.text=element_text(size=12),axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 20,face="bold"),
        panel.grid.major = element_line(color = "gray",linetype = "dotted"), legend.title=element_blank())+
  xlim(-2,2)+ylim(-2,2)+ylab("edger polysome-subpolysome")+xlab("SILAC protein amounts")+labs(title="comparison of the two approaches")+
  geom_abline(size=1, linetype="dashed")