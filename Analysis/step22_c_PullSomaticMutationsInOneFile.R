library(readr)


clinical <- read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' ,header = TRUE,sep='\t')
combined_df<-read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step22data/allMutationsInOneFile.txt' ,header = TRUE,sep='\t')



clinical_sort<-clinical[match(combined_df$sample,clinical$Study_ID),]
combined_df$CaseControl<-clinical_sort$CaseControl
#attach Usag
write.table(combined_df,paste('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step22data/','allMutationsInOneFile_withCaseControl.csv',sep=''), col.names = TRUE,row.names = FALSE,sep = ',' )




combined_df<-read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step22data/allMutationsInOneFile_withCaseControl.txt' ,header = TRUE,sep='\t')
