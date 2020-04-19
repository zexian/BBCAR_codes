library(SKAT)
library(readr)
library('biomaRt')
library(survminer)
library(survival)

GermlinMatrix_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/Somatic_Matrix/'
Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
outputFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step13data/'


mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#load other databases
direct='/projects/p30007/Zexian/Alignment/TCGA_Germline_Breast/'
Gtex<-read_delim(paste(direct,'Gtex_se.txt',sep=''), "\t", escape_double = FALSE, trim_ws = TRUE)

clinical = read.table(Clinical_files,sep='\t',header=TRUE)
readinFile=paste(GermlinMatrix_files,'TTN/number.csv',sep='')
sort_hel<-read_delim(readinFile, ",", escape_double = FALSE, trim_ws = TRUE)
sort_clinical<-clinical[match(unlist(unname(sort_hel[,1])),clinical$Study_ID),]

sample_list = sort_clinical$Study_ID
allGenes=list.files(GermlinMatrix_files)

outcome<-ifelse(sort_clinical$CaseControl=='Case',1,0)
results<-c('gene_name','N_mutation','N_mutation_in_Case','N_mutation_in_Control','P_value_SKAT','Effect(glm)','sig_sum_pvalue_glm','HR(survival)','Survival_P_value_survival')
obj<-SKAT_Null_Model(outcome ~ sort_clinical$Meno+sort_clinical$Benign_Age+sort_clinical$Class, out_type="D")

for (gene in allGenes){
  file<-paste(GermlinMatrix_files,gene,'/number.csv',sep='')
  if (file.exists(file)) {
    print(gene)
    geneMatrix<-read_delim(file, ",", escape_double = FALSE, trim_ws = TRUE)
    geneMatrix_ordered<-geneMatrix[match(sort_clinical$Study_ID,unlist(unname(geneMatrix[,1]))),]
    geneMatrix_ordered<-geneMatrix_ordered[,2:NCOL(geneMatrix_ordered)]
    geneMatrix_ordered<-as.data.frame(geneMatrix_ordered)
    #geneMatrix_ordered<-as.data.frame(geneMatrix_ordered[,colSums(geneMatrix_ordered) < 0.2*NROW(geneMatrix_ordered)])
    #calculate number in case/control
    if(NCOL(geneMatrix_ordered) >=1 ) {
      M_Number=sum(geneMatrix_ordered)
      idx=which(sort_clinical$CaseControl=='Case')
      idx_control=which(sort_clinical$CaseControl=='Control')

      print(NCOL(geneMatrix_ordered))
      
      if (NCOL(geneMatrix_ordered)==1){
        M_Number_case=sum(geneMatrix_ordered[idx,])
        M_Number_control=M_Number-M_Number_case
      }else{
        M_Number_case=length(which(rowSums(geneMatrix_ordered[idx,])>0))
        M_Number_control=length(which(rowSums(geneMatrix_ordered[idx_control,])>0))      
      }

      result<-c(gene,M_Number,M_Number_case,M_Number_control)
      results<-rbind(results,result)

    }
  }
}


  write.table(results,paste(outputFolder,'GeneStudy_number.csv',sep=''), col.names = FALSE,row.names = FALSE,sep = ',' )
  results<-read_delim(paste(outputFolder,'GeneStudy_number.csv',sep=''), ",", escape_double = FALSE, trim_ws = TRUE)
  