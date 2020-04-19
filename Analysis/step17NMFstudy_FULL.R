library(doParallel)
cl<-makeCluster(20)
registerDoParallel(cl)
library(readr)
library(NMF)
#nmf.options(shared.memory=FALSE)
args<-commandArgs(TRUE)

MatrixFile='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative//Step16data/change_context.csv'
SignatureSubFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative//Step17data/change_context/'
MatrixFile=args[1]
SignatureSubFolder=args[2]

temp <- read.csv(MatrixFile)
InputMatrix <- data.frame(temp[,-1], row.names=unlist(temp[,1]))

#InputMatrix<-InputMatrix[rowSums(InputMatrix) > quantile(rowSums(InputMatrix), 0.02), ]
#InputMatrix<-InputMatrix[rowSums(InputMatrix) < quantile(rowSums(InputMatrix), 0.98), ]
InputMatrix = InputMatrix[,colSums(InputMatrix) > 0]

comparison <- read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' ,header = TRUE,sep='\t')
print(dim(comparison))
classes<-comparison$CaseControl
u2 <- unique(classes)
u3<-as.factor(c(as.character(u2),'all'))


for (class_type in u3){
	print('start')
	print(class_type)
	if ( class_type %in% c('Case','Control')){
		sub_group<-comparison[comparison$CaseControl==class_type,]$Study_ID
	}	
	if ( class_type == 'all'){
		sub_group<-comparison$Study_ID
	}
	sub_InputMatrix<-InputMatrix[which(gsub('\'','',rownames(InputMatrix))%in%sub_group),]
	estim.r <- nmf(sub_InputMatrix, 2:6, "nsnmf", nrun = 30, seed='nndsvd',.opt='p', .pbackend=NULL)
	tiff(paste(SignatureSubFolder,class_type,'_Signature_Number.tif',sep=''), width = 1900, height = 1200, units = 'px', res = 250)
	plot(estim.r)
	dev.off()

	for (i in 2:6){
		SignatureSubSubFolder=paste(SignatureSubFolder,'/N',i,'Signatures/',sep='')
		 res <- nmf(sub_InputMatrix,i,"nsnmf",seed='nndsvd',nrun = 100,.opt='p', .pbackend=NULL) # sparse method (NMF)
		 saveRDS(res,paste(SignatureSubSubFolder,class_type,'saved_res.rds',sep=''))
		
		#res<-readRDS(paste(SignatureSubSubFolder,'saved_res.rds',sep=''))
		h <- coef(res) # H  subject feature matrix
		w <- basis(res)
		write.csv(w,paste(SignatureSubSubFolder,class_type,'_w_matrix.csv',sep=''))
		write.csv(t(h),paste(SignatureSubSubFolder,class_type,'_h_matrix.csv',sep=''))

	}
}





