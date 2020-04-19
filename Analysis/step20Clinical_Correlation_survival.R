library(readr)
library(plyr)
library('biomaRt')

#w matrix from home_made pipeline 

w_matrix = read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step17data/change_context/N3Signatures/all_w_matrix.csv',sep=',',header=TRUE)

Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
clinical = read.table(Clinical_files,sep='\t',header=TRUE)
clinical_sort<-clinical[match(gsub('\'','',w_matrix[,1]),clinical$Study_ID),]

#set triple negative cancers 

temp<-cbind.data.frame(w_matrix,clinical_sort)

#test canse/control
fit<-glm(as.factor(CaseControl)~V1+V2+V3+Meno+Benign_Age+Class,data=temp,family = 'binomial')
summary(fit)

#test basal
temp$Basal<-ifelse(temp$ER_Sta=='Negative'&temp$PR_Sta=='Negative'&temp$HER_Sta=='Negative',1,0)
tem_cancer<-temp[temp$CaseControl=='Case',]
tem_cancer<-temp[temp$Basal==1 | temp$CaseControl=='Control' ,]
tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step20data/HomeMade_W_V3_log.tif', width = 800, height = 800, units = 'px', res = 180)
hist(log(temp$V3*100+1))
#hist(temp$V3)
dev.off()

fit<-glm(as.factor(Basal)~V3,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V3+Meno+Benign_Age+Class,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~log(V3*100+1),data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~log(V3*100+1)+Meno+Benign_Age+Class,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V3,data=tem_cancer,family = 'binomial')
fit<-glm(as.factor(Basal)~V3+Meno+Benign_Age+Class,data=tem_cancer,family = 'binomial')
summary(fit)
mean(temp[temp$Basal==1,'V3'])
mean(temp[temp$Basal==0,'V3'])


temp$HER2<-ifelse(temp$HER_Sta=='Positive',1,0)
tem_cancer<-temp[temp$CaseControl=='Case',]
fit<-glm(as.factor(HER2)~V2,data=temp,family = 'binomial')
fit<-glm(as.factor(HER2)~V3+Meno+Benign_Age+Class,data=temp,family = 'binomial')
fit<-glm(as.factor(HER2)~V1,data=tem_cancer,family = 'binomial')
fit<-glm(as.factor(HER2)~V2+Meno+Benign_Age+Class,data=tem_cancer,family = 'binomial')
summary(fit)

#w matrix from another pipeline 

w_matrix = read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step14data/wMatrix.csv',sep='\t',header=TRUE)
colnames(w_matrix)<-c('V1','V2','V3')
rownames(w_matrix)<-gsub('all','',rownames(w_matrix))
Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
clinical = read.table(Clinical_files,sep='\t',header=TRUE)
clinical_sort<-clinical[match(rownames(w_matrix),clinical$Study_ID),]
temp<-cbind.data.frame(w_matrix,clinical_sort)
temp$Basal<-ifelse(temp$ER_Sta=='Negative'&temp$PR_Sta=='Negative'&temp$HER_Sta=='Negative',1,0)
tem_cancer<-temp[temp$CaseControl=='Case',]
tem_cancer<-temp[temp$Basal==1 | temp$CaseControl=='Control' ,]
#test if it is normal distribution 
tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step20data/NMF_W_V3_log.tif', width = 800, height = 800, units = 'px', res = 180)
#hist(exp(temp$V2*100+1))
hist(log(temp$V3*1500+1))
dev.off()

fit<-glm(as.factor(CaseControl)~V2,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~log(V3*1000+1)+Meno+Benign_Age+Class,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V1,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V1+Meno+Benign_Age+Class,data=tem_cancer,family = 'binomial')
summary(fit)




#w matrix from ALE W Matrix pipeline 
w_matrix = read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step23data/N5Signatures/all_w_matrix.csv',sep=',',header=TRUE,row.names=1)
#colnames(w_matrix)<-c('V1','V2','V3')
rownames(w_matrix)<-gsub('all','',rownames(w_matrix))
Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
clinical = read.table(Clinical_files,sep='\t',header=TRUE)
clinical_sort<-clinical[match(rownames(w_matrix),clinical$Study_ID),]
temp<-cbind.data.frame(w_matrix,clinical_sort)
temp$Basal<-ifelse(temp$ER_Sta=='Negative'&temp$PR_Sta=='Negative'&temp$HER_Sta=='Negative',1,0)
tem_cancer<-temp[temp$CaseControl=='Case',]
tem_cancer<-temp[temp$Basal==1 | temp$CaseControl=='Control' ,]
#test if it is normal distribution 
tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step20data/NMF_W_V3_log.tif', width = 800, height = 800, units = 'px', res = 180)
#hist(exp(temp$V2*100+1))
hist(log(temp$V3*1500+1))
dev.off()

fit<-glm(as.factor(CaseControl)~V1,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~V2,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~V3,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~V4,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~V5,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(CaseControl)~V6,data=temp,family = 'binomial')
summary(fit)


fit<-glm(as.factor(Basal)~V1,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V2,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V3,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V4,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V5,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V6,data=tem_cancer,family = 'binomial')
summary(fit)


fit<-glm(as.factor(CaseControl)~log(V1*1000+1)+Meno+Benign_Age+Class,data=temp,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V3,data=tem_cancer,family = 'binomial')
summary(fit)
fit<-glm(as.factor(Basal)~V3+Meno+Benign_Age+Class,data=tem_cancer,family = 'binomial')
summary(fit)




#Now study copy number variation
CNA_files<-'/projects/p30007/Zexian/Alignment/Germline_37/CNV/VarScan_Process/All/all_thresholded.by_genes.txt'
counts = read.delim(CNA_files, as.is = T, row.names=1,sep = "\t", skip = 0)
attache_later<-counts[,1:2]
attache_later$geneName<-sapply(strsplit(rownames(attache_later),"\\|"), `[`, 1)  
counts<-counts[,-1]
counts<-counts[,-1]
counts<-t(counts)
rownames(counts)<-gsub('X','',rownames(counts))

Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
clinical = read.table(Clinical_files,sep='\t',header=TRUE)
clinical_sort<-clinical[match(rownames(counts),clinical$Study_ID),]
clinical_sort$HER2<-ifelse(clinical_sort$HER_Sta=='Positive',1,0)
clinical_sort$CaseControl=as.factor(ifelse(clinical_sort$CaseControl=='Case',1,0))

dataset<-c()
for (Ensembl_ID in colnames(counts)){
  gene_expre<-counts[,Ensembl_ID]
  temp_data<-cbind.data.frame(input=gene_expre,output=clinical_sort$CaseControl,HER2=clinical_sort$HER2)
  cor_resu<- try(glm(output~input,data=temp_data,family='binomial'),silent=TRUE)
  if(class(cor_resu) != "try-error"){
        print(Ensembl_ID)
  }else{
    print('wrong')
  }
  if (typeof(cor_resu) =='list' & dim(coef(summary(cor_resu)))[1] >1 ){
      result<-c(mean(temp_data$input,na.rm=TRUE),mean(temp_data$input[clinical_sort$CaseControl==1],na.rm=TRUE),mean(temp_data$input[clinical_sort$CaseControl==0],na.rm=TRUE),coef(summary(cor_resu))[2,1],coef(summary(cor_resu))[2,4])
      dataset<-rbind(dataset,result)
      rownames(dataset)[NROW(dataset)]<-Ensembl_ID
    }else{
      result<-c(mean(temp_data$input,na.rm=TRUE),mean(temp_data$input[clinical_sort$CaseControl==1],na.rm=TRUE),mean(temp_data$input[clinical_sort$CaseControl==0],na.rm=TRUE),0,1)
      dataset<-rbind(dataset,result)
      rownames(dataset)[NROW(dataset)]<-Ensembl_ID
    }
}
colnames(dataset)<-c('Mean_Copy','mean_in_Case','mean_in_Control','Effect','Pvalue')
dataset<-as.data.frame.matrix(dataset)

#start some annotations
gene_anno_filtered = sapply(strsplit(rownames(dataset),"\\|"), `[`, 1)    

mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
dat = getBM(
  values = gene_anno_filtered,
  filters = c("external_gene_name"),
  attributes = c("ensembl_transcript_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position"),
  mart = mart
)
idx_gene = match(gene_anno_filtered, dat$external_gene_name)
gene_anno_filtered<-as.data.frame(gene_anno_filtered)
gene_anno_filtered$external_gene_name=dat$external_gene_name[idx_gene]
gene_anno_filtered$Description=dat$description[idx_gene]
gene_anno_filtered$chromosome_name=dat$chromosome_name[idx_gene]
gene_anno_filtered$start_position=dat$start_position[idx_gene]
gene_anno_filtered$end_position=dat$end_position[idx_gene]

#gene_anno_sorted<-gene_anno_filtered[match(temp$Transcript_ID,gene_anno_filtered$gene_anno_filtered),]
annotations_sorted<-gene_anno_filtered[match(sapply(strsplit(rownames(dataset),"\\|"), `[`, 1),gene_anno_filtered$external_gene_name),]
out_sub<-cbind.data.frame(dataset,annotations_sorted)

out_sub$FDR<-p.adjust(out_sub[,'Pvalue'], method = "fdr", n = NROW(out_sub))
out_sub<-out_sub[ order(out_sub$FDR), ]
attache_later_sort<-attache_later[match(out_sub$gene_anno_filtered,attache_later$geneName),]
#take every three and save to matrix 
out_sub<-cbind.data.frame(out_sub,attache_later_sort)

write.csv(out_sub,paste('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step20data/','CNV_Correlation_result_CASE_Control.csv',sep=''))





#summary the clinical cohort distribution 
library(ggpubr)
library(ggplot2)
library(gmodels)


#end of the story
Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
comparison = read.table(Clinical_files,sep='\t',header=TRUE)
table(comparison$CaseControl)
CASE_group<-comparison[comparison$CaseControl=='Case',]
Control_group<-comparison[comparison$CaseControl=='Control',]
case_control<-comparison$CaseControl

t.test(CASE_group$Benign_Age,Control_group$Benign_Age)    
sd(CASE_group$Benign_Age) 
sd(Control_group$Benign_Age)


CrossTable(comparison$Meno,comparison$CaseContro,prop.t=FALSE, prop.r=FALSE, prop.c=TRUE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(comparison$Class,comparison$CaseContro,prop.t=FALSE, prop.r=FALSE, prop.c=TRUE,chisq=TRUE,prop.chisq=FALSE )
table(CASE_group$ER_Sta)

CrossTable(comparison$has_Germline,comparison$CaseContro,prop.t=FALSE, prop.r=FALSE, prop.c=TRUE,chisq=TRUE,prop.chisq=FALSE )
CrossTable(comparison$has_Cancer,comparison$CaseContro,prop.t=FALSE, prop.r=FALSE, prop.c=TRUE,chisq=TRUE,prop.chisq=FALSE )



t.test(CASE_group$interval_days/365,Control_group$interval_days/365)    
sd(CASE_group$interval_days)/365
sd(Control_group$interval_days)/365














#laod clinical 
clinical_final <- read_delim(paste(Clinical_files,'/data_clinical_patient_atlas2017.txt',sep=''), "\t", escape_double = FALSE, trim_ws = TRUE,skip=4)
rownames(clinical_final)<-unlist(clinical_final[,1])
clinical_sort<-clinical_final[match(gsub('\'','',w_matrix[,1]),rownames(clinical_final)),]

basal<-ifelse(clinical_sort$SUBTYPE=='BRCA_Basal',1,0)
data<-cbind.data.frame(basal,sig=w_matrix[,5])
fit<-glm(as.factor(basal)~sig,data=data,family = 'binomial')
summary(fit)

#fit survival model
library(survminer)
library(survival)
temp1<-cbind.data.frame(status=clinical_sort$`OS event`,time=clinical_sort$`OS Time`,type=signatureB$sig  )
temp1<-temp1[which(!is.na(temp1$time)),]
coxph_obj = coxph(Surv(temp1$time, temp1$status) ~ temp1$type)
summary(coxph_obj)
#try if the signature is is correlate with age
temp1<-cbind.data.frame(age=clinical_sort$age_at_diagnosis,type=signatureB$sig  )
cor.test(temp1$age,temp1$type)
fit<-lm(type~age,data=temp1)
summary(fit)
#############


#make numbers for figure 2a basal
temp<-cbind.data.frame(w_matrix,sub=clinical_sort$SUBTYPE)
tt1<-temp[which(temp$sub=='BRCA_Basal'),'V1']
tt2<-temp[which(temp$sub=='BRCA_LumA'),'V1']
tt3<-temp[which(temp$sub=='BRCA_LumB'),'V1']
tt4<-temp[which(temp$sub=='BRCA_Her2'),'V1']

tt5<-temp[which(temp$sub=='BRCA_Basal'),'V2']
tt6<-temp[which(temp$sub=='BRCA_LumA'),'V2']
tt7<-temp[which(temp$sub=='BRCA_LumB'),'V2']
tt8<-temp[which(temp$sub=='BRCA_Her2'),'V2']

tt9<-temp[which(temp$sub=='BRCA_Basal'),'V4']
tt10<-temp[which(temp$sub=='BRCA_LumA'),'V4']
tt11<-temp[which(temp$sub=='BRCA_LumB'),'V4']
tt12<-temp[which(temp$sub=='BRCA_Her2'),'V4']

tt13<-temp[which(temp$sub=='BRCA_Basal'),'V3']
tt14<-temp[which(temp$sub=='BRCA_LumA'),'V3']
tt15<-temp[which(temp$sub=='BRCA_LumB'),'V3']
tt16<-temp[which(temp$sub=='BRCA_Her2'),'V3']

temp<-cbind.data.frame(tt1*1000,tt5*1000,tt9*1000,tt13*1000)
write.table(temp,'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step34codes_some_logistics/Basal.csv',col.names=FALSE,row.names=FALSE,sep=',')
temp<-cbind.data.frame(tt2*1000,tt6*1000,tt10*1000,tt14*1000)
write.table(temp,'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step34codes_some_logistics/LumA.csv',col.names=FALSE,row.names=FALSE,sep=',')
temp<-cbind.data.frame(tt3*1000,tt7*1000,tt11*1000,tt15*1000)
write.table(temp,'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step34codes_some_logistics/LumB.csv',col.names=FALSE,row.names=FALSE,sep=',')
temp<-cbind.data.frame(tt4*1000,tt8*1000,tt12*1000,tt16*1000)
write.table(temp,'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step34codes_some_logistics/Her2.csv',col.names=FALSE,row.names=FALSE,sep=',')

#make numbers for figure 2b ERPRHER2
#laod clinical 
clinical_final2 <- read_delim(paste(Clinical_files,'/data_clinical_sample.txt',sep=''), "\t", escape_double = FALSE, trim_ws = TRUE,skip=4)
clinical_sort2<-clinical_final2[match(gsub('\'','',w_matrix[,1]),clinical_final2$ID),]
temp<-cbind.data.frame(w_matrix,ER=clinical_sort2$ER_STATUS_BY_IHC,PR=clinical_sort2$PR_STATUS_BY_IHC,HER2=clinical_sort2$HER2_STATUS)

#get V1 statistics
tt1<-temp[which(temp$ER=='Positive'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Positive'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Positive'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))

tt1<-temp[which(temp$ER=='Negative'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Negative'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Negative'),'V1']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))



#get V2 statistics
tt1<-temp[which(temp$ER=='Positive'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Positive'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Positive'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))

tt1<-temp[which(temp$ER=='Negative'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Negative'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Negative'),'V2']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))



#get V4 statistics
tt1<-temp[which(temp$ER=='Positive'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Positive'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Positive'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))

tt1<-temp[which(temp$ER=='Negative'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Negative'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Negative'),'V4']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))



#get V3 statistics
tt1<-temp[which(temp$ER=='Positive'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Positive'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Positive'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))

tt1<-temp[which(temp$ER=='Negative'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$PR=='Negative'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))
tt1<-temp[which(temp$HER2=='Negative'),'V3']
print(mean(tt1)*1000)
print(sd(tt1)*1000)
print(length(tt1))







for (line in tt){
	print(line*1000)
}

#fit linear
for( i in 2:NCOL(clinical_sort)){
temp1<-cbind.data.frame(outcome=clinical_sort[,i],type=w_matrix[,4]  )
colnames(temp1)<-c('outcome','type')
temp1<-temp1[which(!is.na(temp1$outcome)),]
P_value<-1
fit<-lm(c(1,1)~c(1,1))
tryCatch(  fit<-lm(outcome~type,data=temp1),
           warning = function(w) {print('w')},
           error = function(e) {print('e');NaN}
)
P_value<-summary(fit)$coefficients[,4][2] 
if( is.na(P_value) ){
  P_value<-1
}
if (P_value<0.05) {
  print('bingo')
  print(i)
  print(colnames(clinical_sort)[i])
  print(P_value)
}
tryCatch(  fit<-glm(as.factor(temp1$outcome)~temp1$type,data=temp1,family=binomial(link='logit')),
           warning = function(w) {print('w')},
           error = function(e) {print('e');NaN}
)

P_value<-summary(fit)$coefficients[,4][2] 
if( is.na(P_value) ){
  P_value<-1
}
if (P_value<0.05) {
  print('bingo')
  print(i)
  print(colnames(clinical_sort)[i])
  print(P_value)
}
}

#save the two boxplots to txt files 
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,type=signatureB$sig  )
basal<-temp1[which(temp1$outcome=='Basal'),2]*10000
LumA<-temp1[which(temp1$outcome=='LumA'),2]*10000
LumB<-temp1[which(temp1$outcome=='LumB'),2]*10000
Her2<-temp1[which(temp1$outcome=='Her2'),2]*10000
Norm<-temp1[which(temp1$outcome=='Normal'),2]*10000
write.table(basal,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/SigB/Basal.txt',col.names = FALSE,row.names = FALSE,sep = ',')
write.table(LumA,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/SigB/LumA.txt',col.names = FALSE,row.names = FALSE,sep = ',')
write.table(LumB,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/SigB/LumB.txt',col.names = FALSE,row.names = FALSE,sep = ',')
write.table(Her2,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/SigB/Her2.txt',col.names = FALSE,row.names = FALSE,sep = ',')
write.table(Norm,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/SigB/Norm.txt',col.names = FALSE,row.names = FALSE,sep = ',')

#save the clusters
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,SigA=signatureA$sig,SigD=signatureD$sig  )
library('ggplot2')
write.table(temp1,'/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/sig1Sig4Cluster.txt',col.names = FALSE,row.names = FALSE,sep = ',')
ggplot(temp1, aes(x = SigA*1000, y = SigD, col = outcome)) +
  geom_point() 

#all four together to predict pam
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,SigA=signatureA$sig,SigB=signatureB$sig,SigC=signatureC$sig,SigD=signatureD$sig )
fit<-glm(as.factor(temp1$outcome)~temp1$SigA+temp1$SigB+temp1$SigC+temp1$SigD,data=temp1,family=binomial(link='logit'))
summary(fit)
#draw the pCA plot 
require(genefilter)
require(calibrate)
require(RColorBrewer)
library(ggplot2)
library(ggfortify)
library(ggrepel)
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,SigA=signatureA$sig,SigB=signatureB$sig,SigC=signatureC$sig,SigD=signatureD$sig )
pdfname <-paste("/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/PCA_RPKM.pdf",sep="")
print(pdfname)
pca=prcomp(temp1[,2:5])
pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
pdf(pdfname)
autoplot(pca,data=temp1,colour='outcome')+
  geom_text_repel(aes_string(x = pc1lab, y =pc2lab, label = "Comparisons"))
dev.off()
fit<-glm(as.factor(temp1$outcome)~temp1$SigD,data=temp1,family=binomial(link='logit'))
summary(fit)


#draw heatmap 
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,SigA=signatureA$sig,SigB=signatureB$sig,SigC=signatureC$sig,SigD=signatureD$sig )
temp1<-temp1[ order( temp1[,1]), ]
temp1<-temp1[which(!is.na(temp1$outcome)),]
matrix<-as.matrix(temp1[,2:5]*1000)
rownames(matrix)<-as.character(temp1$outcome)
tiff('/Users/zza847/Dropbox/AminoAcidAnalysis/Text2/WorkingFolder/Figure2/CodesData/HeatmapForSig.tif', width = 1500, height = 2800, units = 'px', res = 180)
hh<-heatmap.2(matrix, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x)hclust(x,method="ward.D"), col=hmcols, trace= "none", main = "correaltion distance",scale = c('row'))
dev.off()
sorted <- matrix[match(rev(labels(hh$rowDendrogram)), rownames(matrix)), ]


#fit and test
temp1<-cbind.data.frame(outcome=clinical_sort$`PR Status.y`,type=signatureD$sig  )
colnames(temp1)<-c('outcome','type')
temp1<-temp1[which(!is.na(temp1$outcome)),]
fit<-glm(as.factor(temp1$outcome)~temp1$type,data=temp1,family=binomial(link='logit'))
summary(fit)

#fit and test
temp1<-cbind.data.frame(outcome=clinical_sort$PAM50,type=signatureD$sig  )
temp1 <- within(temp1, outcome <- relevel(outcome, ref = 'LumA'))
colnames(temp1)<-c('outcome','type')
temp1<-temp1[which(!is.na(temp1$outcome)),]
fit<-lm(temp1$type~as.factor(temp1$outcome))
summary(fit)

boxplot(type~outcome,
        data=temp1,
        main="Different boxplots",
        xlab="Month Number",
        ylab="Degree Fahrenheit",
        col="orange",
        border="brown"
)
means <- tapply(temp1$type,temp1$outcome,mean)
points(means,col="red",pch=18)
means*1000


for(input in inputs) {
       tryCatch(  print(paste("log of", input, "=", log(input))),
                     warning = function(w) {print('w')},
                     error = function(e) {print('e');NaN}
              )
}



#correlate the signature with MRNA data. 
signatureA<-cbind.data.frame(rownames(w),sig=w[,1])
signatureB<-cbind.data.frame(rownames(w),sig=w[,2])
signatureC<-cbind.data.frame(rownames(w),sig=w[,3])
signatureD<-cbind.data.frame(rownames(w),sig=w[,4])

TCGA_RNA<- read_csv("~/Dropbox/AminoAcidAnalysis/WorkingFolder/Data/TCGA_RNA/TCGA_RNA_Big_Table_WithGeneName.csv",na = 'NA')
GeneList<-unique(TCGA_RNA$gene_sub)
mean_sigA<-mean(signatureA$sig)
mean_sigB<-mean(signatureB$sig)
mean_sigC<-mean(signatureC$sig)
mean_sigD<-mean(signatureD$sig)
results<-c('gene_name','average RNA expression','Average signature_A','correlatoinA','pvalue_percentageA','Average signature_B','correlatoinB','pvalue_percentageB','Average signature_C','correlatoinC','pvalue_percentageC','Average signature_D','correlatoinD','pvalue_percentageD')
for (gene in GeneList){
  if (!is.na(gene)){
  #print(gene)
  if (length(which(TCGA_RNA$gene_sub==gene))==1){
  TCGA_RNA_sub<-t(TCGA_RNA[which(TCGA_RNA$gene_sub==gene),][,6:NCOL(TCGA_RNA)])
  TCGA_RNA_sub<-data.frame(Samplename=row.names(TCGA_RNA_sub),TCGA_RNA_sub)
  TCGA_RNA_have_Score<-as.data.frame(TCGA_RNA_sub[match(gsub("\'","",rownames(signatureA)),row.names(TCGA_RNA_sub)),])
  if( sum(TCGA_RNA_have_Score$TCGA_RNA_sub,na.rm=TRUE)>0.256*1092  ) {
  test_percentageA<-cor.test(signatureA$sig,TCGA_RNA_have_Score$TCGA_RNA_sub,use = "complete.obs")
  test_percentageB<-cor.test(signatureB$sig,TCGA_RNA_have_Score$TCGA_RNA_sub,use = "complete.obs")
  test_percentageC<-cor.test(signatureC$sig,TCGA_RNA_have_Score$TCGA_RNA_sub,use = "complete.obs")
  test_percentageD<-cor.test(signatureD$sig,TCGA_RNA_have_Score$TCGA_RNA_sub,use = "complete.obs")
  result<-c(gene,mean(TCGA_RNA_have_Score$TCGA_RNA_sub,na.rm=TRUE),mean_sigA,test_percentageA$estimate,test_percentageA$p.value,mean_sigB,test_percentageB$estimate,test_percentageB$p.value,mean_sigC,test_percentageC$estimate,test_percentageC$p.value,mean_sigD,test_percentageD$estimate,test_percentageD$p.value)
  results<-rbind.data.frame(results,result)
  }
  }
  }
}
results<-as.data.frame.matrix(results)
results <- results[-1, ] 
colnames(results) <- results[3,]
dat = getBM(
  values = as.character(results$),
  filters = c("external_gene_name"),
  attributes = c("external_gene_name", "description"),
  mart = mart
)
idx_gene = match(results[,1], dat$external_gene_name)
results$Description=dat$description[idx_gene]

write.table(results,"/Users/zza847/Dropbox/AminoAcidAnalysis/Text/WorkingFolder/GeneList/signature_with_RNA_List/signature_with_RNA_List.csv", col.names = FALSE,row.names = FALSE,sep = ',' )






