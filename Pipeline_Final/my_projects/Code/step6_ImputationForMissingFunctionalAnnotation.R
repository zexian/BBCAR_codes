library(plyr)
library(stringr)
library(mice)
BBCAR<-'/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic/'
GERMline37<-'/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Mutect/Germline_Somatic_VCF/'

OutFigure='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step6data/'

train_files_ls = list.files(path=GERMline37, pattern="*.txt") 
BBCAR_files_ls = list.files(path=BBCAR, pattern="*.txt") 

test_set<-BBCAR_files_ls[!BBCAR_files_ls %in% train_files_ls]
train_set<-train_files_ls

train_samples<-gsub('.hg19_multianno.txt','',train_set)
test_samples<-gsub('.hg19_multianno.txt','',test_set)

txt_files_ls<-c(paste(GERMline37,train_set,sep=''),paste(BBCAR,test_set,sep=''))

#read in train set: 
for (file in txt_files_ls) {
    Sample=gsub('.hg19_multianno.txt','',file)
    Sample=gsub(GERMline37,'',Sample)
    Sample=gsub(BBCAR,'',Sample)
    print(file)
    if (!exists('titles')){
        titles=as.character(unname(unlist(read.csv(file,sep='\t',nrows=1,header=FALSE))))
        titles<-c('Sample',titles,'V1','V2','V3','V4','rs_number','from','to','V5','class','SNPEFF','GT_name','GT_tumor')
      } 
    if (!exists('dfs')){
          print('in') 
      dfs = read.csv(file,sep='\t',skip=1,header=FALSE)
      dfs<-dfs[,1:114]
      print(dim(dfs))
      dfs<-cbind(Sample,dfs)
      } 
    if(exists("dfs")){
      df_temp = read.csv(file,sep='\t',skip=1,header=FALSE)
      df_temp<-df_temp[,1:114]
      print('here')
      print(dim(df_temp))
      df_temp<-cbind(Sample,df_temp)
      dfs = rbind(dfs, df_temp)
    }
}
colnames(dfs)<-titles
dfs$key<-paste(dfs$Chr,dfs$Start,sep='_')
colnames(dfs)[NCOL(dfs)]<-'key'
count_talbe<-as.data.frame(count(dfs, "key"))
count_talbe<-count_talbe[match(dfs$key,count_talbe$key),]
dfs$fre<-count_talbe$freq
dfs$ale<-sapply(strsplit(as.character(dfs$GT_tumor),":"), `[`, 3)
dfs$ref_depth<-sapply(strsplit(as.character(dfs$GT_tumor),":"), `[`, 2)
dfs$ref_depth<-sapply(strsplit(as.character(dfs$ref_depth),","), `[`, 1)
dfs$GERP<-dfs$'GERP++_RS'
print('the number of rows are: ')
print(dim(dfs))

dfs_select<-dfs[,c('Sample','class','ale','ref_depth','fre','ExAC_ALL','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred','MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','MetaSVM_score','MetaSVM_pred','MetaLR_score','MetaLR_pred','VEST3_score','CADD_raw','CADD_phred','GERP','phyloP20way_mammalian','phyloP100way_vertebrate','SiPhy_29way_logOdds','snp138NonFlagged','cosmic80')]
dfs_select$snp138NonFlagged<-ifelse(dfs_select$snp138NonFlagged=='.',0,1)
dfs_select$cosmic80<-ifelse(dfs_select$cosmic80=='.',0,1)
dfs_select[dfs_select=='.']<-NA
dfs_select <- droplevels(dfs_select)
#dfs_select<-dfs_select[1:300,]
dfs_select$ale<-as.numeric(dfs_select$ale)
dfs_select$ref_depth<-as.numeric(dfs_select$ref_depth)

if (!is.null(levels(dfs_select$ExAC_ALL))){
  dfs_select$ExAC_ALL<-as.numeric(levels(dfs_select$ExAC_ALL))[dfs_select$ExAC_ALL] 
}
if (!is.null(levels(dfs_select$SIFT_score))){
  dfs_select$SIFT_score<-as.numeric(levels(dfs_select$SIFT_score))[dfs_select$SIFT_score] 
}
if (!is.null(levels(dfs_select$Polyphen2_HDIV_score))){
  dfs_select$Polyphen2_HDIV_score<-as.numeric(levels(dfs_select$Polyphen2_HDIV_score))[dfs_select$Polyphen2_HDIV_score] 
}
if (!is.null(levels(dfs_select$Polyphen2_HVAR_score))){
  dfs_select$Polyphen2_HVAR_score<-as.numeric(levels(dfs_select$Polyphen2_HVAR_score))[dfs_select$Polyphen2_HVAR_score] 
}
if (!is.null(levels(dfs_select$LRT_score))){
  dfs_select$LRT_score<-as.numeric(levels(dfs_select$LRT_score))[dfs_select$LRT_score] 
}
if (!is.null(levels(dfs_select$MutationTaster_score))){
  dfs_select$MutationTaster_score<-as.numeric(levels(dfs_select$MutationTaster_score))[dfs_select$MutationTaster_score] 
}
if (!is.null(levels(dfs_select$MutationAssessor_score))){
  dfs_select$MutationAssessor_score<-as.numeric(levels(dfs_select$MutationAssessor_score))[dfs_select$MutationAssessor_score] 
}
if (!is.null(levels(dfs_select$FATHMM_score))){
  dfs_select$FATHMM_score<-as.numeric(levels(dfs_select$FATHMM_score))[dfs_select$FATHMM_score] 
}
if (!is.null(levels(dfs_select$MetaSVM_score))){
  dfs_select$MetaSVM_score<-as.numeric(levels(dfs_select$MetaSVM_score))[dfs_select$MetaSVM_score] 
}
if (!is.null(levels(dfs_select$MetaLR_score))){
  dfs_select$MetaLR_score<-as.numeric(levels(dfs_select$MetaLR_score))[dfs_select$MetaLR_score] 
}
if (!is.null(levels(dfs_select$VEST3_score))){
  dfs_select$VEST3_score<-as.numeric(levels(dfs_select$VEST3_score))[dfs_select$VEST3_score] 
}
if (!is.null(levels(dfs_select$CADD_raw))){
  dfs_select$CADD_raw<-as.numeric(levels(dfs_select$CADD_raw))[dfs_select$CADD_raw] 
}
if (!is.null(levels(dfs_select$CADD_phred))){
  dfs_select$CADD_phred<-as.numeric(levels(dfs_select$CADD_phred))[dfs_select$CADD_phred] 
}
if (!is.null(levels(dfs_select$GERP))){
  dfs_select$GERP<-as.numeric(levels(dfs_select$GERP))[dfs_select$GERP] 
}
if (!is.null(levels(dfs_select$phyloP20way_mammalian))){
  dfs_select$phyloP20way_mammalian<-as.numeric(levels(dfs_select$phyloP20way_mammalian))[dfs_select$phyloP20way_mammalian] 
}
if (!is.null(levels(dfs_select$phyloP100way_vertebrate))){
  dfs_select$phyloP100way_vertebrate<-as.numeric(levels(dfs_select$phyloP100way_vertebrate))[dfs_select$phyloP100way_vertebrate] 
}
if (!is.null(levels(dfs_select$SiPhy_29way_logOdds))){
  dfs_select$SiPhy_29way_logOdds<-as.numeric(levels(dfs_select$SiPhy_29way_logOdds))[dfs_select$SiPhy_29way_logOdds] 
}
#imp<-mice(dfs_select,meth=c('','','','','','pmm','logreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','logreg','pmm','logreg','pmm','logreg','pmm','pmm','pmm','pmm','pmm','pmm','pmm'),m=5,maxit=50,seed=500) 
#method<-as.data.frame(imp$meth)
library(parallel)
cores_2_use<-20
cl <- makeCluster(cores_2_use,outfile="")
clusterSetRNGStream(cl, 9956)
clusterExport(cl, "dfs_select")
clusterEvalQ(cl, library(mice))

set.seed(500)
imp_pars <- 
  parLapply(cl = cl, X = 1:cores_2_use, fun = function(no){
    #mice(dfs_select,meth=c('','','','','','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm','pmm'),m=1,maxit=5) 
    mice(dfs_select,meth=c('','','','','','pmm','pmm','logreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','polyreg','pmm','logreg','pmm','logreg','pmm','logreg','pmm','pmm','pmm','pmm','pmm','pmm','pmm','logreg','logreg'),m=1,maxit=5) 
  })

imp <- imp_pars[[1]]
method<-unlist(imp['method'])
for (n in 2:length(imp_pars)){
  imp <- 
    ibind(imp,
          imp_pars[[n]])
}

stopCluster(cl)

names<-gsub('method.','',names(method))
#get continuous number 
IMP_Conti_Pridictors_ave<-(complete(imp, 1)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 2)[,c(names[which(method=='pmm' )])]  
                            +complete(imp, 3)[,c(names[which(method=='pmm' )])]  
                           +complete(imp, 4)[,c(names[which(method=='pmm' )])]  
                           +complete(imp, 5)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 6)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 7)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 8)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 9)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 10)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 11)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 12)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 13)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 14)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 15)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 16)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 17)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 18)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 19)[,c(names[which(method=='pmm' )])]
                           +complete(imp, 20)[,c(names[which(method=='pmm' )])])/20


Categorical_Data <- complete(imp, 'long')[,c('.id',names[which(method=='logreg' | method=='polyreg' )])]
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
IMP_Cate_Pridictors_ave<-aggregate(  .~Categorical_Data$.id, Categorical_Data, Mode)
IMP_Cate_Pridictors_ave <- IMP_Cate_Pridictors_ave[, -c(1,2)] 
complete_data<-cbind.data.frame(dfs_select[,1:5],IMP_Conti_Pridictors_ave,IMP_Cate_Pridictors_ave)

complete_data$class<-ifelse(complete_data$class=='PASS',1,0)

train_set<-complete_data[complete_data$Sample %in% train_samples, ]
write.table(train_set,paste(OutFigure,"train5_update_.txt",sep=''),sep="\t",row.names=FALSE)

test_set<-complete_data[complete_data$Sample %in% test_samples, ]
write.table(test_set,paste(OutFigure,"test5_update_.txt",sep=''),sep="\t",row.names=FALSE)

write.table(complete_data,paste(OutFigure,"complete5_update_.txt",sep=''),sep="\t",row.names=FALSE)


