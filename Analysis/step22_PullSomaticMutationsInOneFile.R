library(SKAT)
library(readr)

args<-commandArgs(TRUE)

GermlinMatrix_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step15data/'
signature_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step17data/change_context/N3Signatures/' 
context_type='change_context' 

outputFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step22data/'


#read what mutaion falls into what signature
sigFile<-paste(signature_files,'all_h_matrix.csv',sep='')
sigNumber<-read.table(sigFile,sep=',',header=TRUE)
sigNumber$maxNum <- apply(sigNumber, 1, which.max)
i<-sigNumber$maxNum
sigNumber$sig=ifelse((i-1)==1,'SigA',ifelse((i-1)==2,'SigB',ifelse((i-1)==3,'SigC',ifelse((i-1)==4,'SigD',ifelse((i-1)==5,'SigE',ifelse((i-1)==6,'SigF','wrong'))))))
sigNumber<-sigNumber[,c('X','sig')]

#read key change number
aa_number_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step33data/AA_change_rate_named.csv'
aa_number<-read.table(aa_number_file,sep=',',header=TRUE)
aa_number<-aa_number[,c('key','number')]
aa_number$number=aa_number$number*100/sum(aa_number$number)

#code to aa
code_to_AA_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step33data/AA_Change.txt'
code_to_AA<-read.table(code_to_AA_file,sep='\t',header=FALSE)


#Polor changed
Polor_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step33data/index2_ChargeGroup.txt'
Polor<-read.table(Polor_file,sep='\t',header=FALSE)

#Usage changed
Usage_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step33data/Codan_USage_Human_two_in_one.txt'
Usage<-read.table(Usage_file,sep='\t',header=FALSE)

#Usage changed
PAM50_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/PAM50Gene_sig.txt'
PAM50<-as.character(unname(unlist(read.table(PAM50_file,sep='\t',header=FALSE))))

#Code score changed
Code_Score_file<-'/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/Step33data/CODE_SCORE.txt'
Code_Score<-read.table(Code_Score_file,sep='\t',header=FALSE)

#load case or control
clinical <- read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' ,header = TRUE,sep='\t')


#read mutation files
files_ls = list.files(path=GermlinMatrix_files, pattern=".hg19_multianno.csv")
txt_files_ls <-paste(GermlinMatrix_files,files_ls,sep='')
txt_files_df <- lapply(txt_files_ls, function(x) {read.table(file = x, header = T, sep ="\t")})
for (i in 1:length(txt_files_df)){txt_files_df[[i]]<-cbind(sample=gsub('.hg19_multianno.csv','',files_ls[i]),txt_files_df[[i]])}

# Combine them

combined_df <- do.call("rbind", lapply(txt_files_df, as.data.frame))
combined_df$change<-paste(combined_df$Ref_after_reverse,combined_df$Alt_after_reverse,combined_df$AA_Code,combined_df$Mutated_AA_Position,combined_df$upstream,combined_df$downstream,sep='_')
combined_df$change<-paste(combined_df$Ref_after_reverse,combined_df$Alt_after_reverse,combined_df$upstream,combined_df$downstream,sep='_')

#guss what signature the mutation falls in
sigNumber_sort<-sigNumber[match(combined_df$change,sigNumber$X),]
combined_df$sig<-sigNumber_sort$sig

#get alt amino acid 
combined_df$ALT_Code<-as.character(combined_df$AA_Code)
substr(combined_df$ALT_Code,as.integer(combined_df$Mutated_AA_Position),as.integer(combined_df$Mutated_AA_Position)) <- as.character(combined_df$Alt_after_reverse)
#now get the matched alt aa
code_to_AA_sort<-code_to_AA[match(combined_df$ALT_Code,code_to_AA$V1),]
combined_df$alt_AA<-code_to_AA_sort$V2
combined_df$key_AA<-paste(combined_df$AA_Name,combined_df$alt_AA,sep='_')

#append the change socre sent by dr. clare
aa_number_sort<-aa_number[match(combined_df$key_AA,aa_number$key),]
combined_df$number_ng345<-aa_number_sort$number

#calculate our own score
counts<-as.data.frame(table(combined_df$key_AA))
counts$percent<-counts$Freq*100/(sum(counts$Freq))
counts_sort<-counts[match(combined_df$key_AA,counts$Var1),]

combined_df$InthisCohortAA_Change_percent<-counts_sort$percent


#attach polor
Polor_sort<-Polor[match(combined_df$AA_Name,Polor$V1),]
combined_df$ref_Polor<-Polor_sort$V2
Polor_sort<-Polor[match(combined_df$alt_AA,Polor$V1),]
combined_df$alt_Polor<-Polor_sort$V2


#attach Usage
Usage_sort<-Usage[match(combined_df$AA_Code,Usage$V1),]
combined_df$ref_Usage<-Usage_sort$V2
Usage_sort<-Usage[match(combined_df$ALT_Code,Usage$V1),]
combined_df$alt_Usage<-Usage_sort$V2


clinical_sort<-clinical[match(combined_df$sample,clinical$Study_ID),]
combined_df$CaseControl<-clinical_sort$CaseControl
#attach Usag
write.table(combined_df,paste(outputFolder,'merged_mutation_AA_FRE.csv',sep=''), col.names = TRUE,row.names = FALSE,sep = ',' )




