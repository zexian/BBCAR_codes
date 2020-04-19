#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#source("https://bioconductor.org/biocLite.R")
#biocLite("SomaticSignatures")
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38")
#install.packages('ggdendro')
#source("https://bioconductor.org/biocLite.R")
#biocLite("signeR")
library(readr)
library(VariantAnnotation)
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggdendro)
library(ggplot2)
library(signeR)
library(rtracklayer)
library(readr)
library(LaplacesDemon)
print('chicken0')
library(readr)
X30_signatures <- read_csv("/projects/p30007/Zexian/tools/DNAtools/30_signatures.csv")[,4:33]

Analysis <- '/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step14data/'
VCF<-'/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
print('chicken1')

file<-paste('/projects/p30007/Zexian/tools/DNAtools/exome_count.txt',sep="")
#opp_file<-read_csv(file,locale=locale(tz="Australia/Sydney"))
opp_file<-read.table(file,header = TRUE,sep='\t')
comparison <- read.table('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' ,header = TRUE,sep='\t')
print(dim(comparison))
classes<-comparison$CaseControl
u2 <- unique(classes)
u3<-as.factor(c(as.character(u2),'all'))
#u2<-1  #need to delete

print('chicken3')
for (class_type in u3){
	print('start')
	print(class_type)
	if ( class_type %in% c('Case','Control')){
		sub_group<-comparison[comparison$CaseControl==class_type,]$Study_ID
	}	
	if ( class_type == 'all'){
		sub_group<-comparison$Study_ID
	}	
	print(length(sub_group))
	all_vcf<-GRanges()
	for (individual in sub_group) {
		if (!is.na(individual) ) {
			print(individual)
		filename <- paste(VCF,'/',individual,'.hg19_multianno.vcf',sep="")
		print(filename)
		in_vcf <- readVcf(filename, "hg19")
		gvcf<-rowRanges(in_vcf)
		gvcf$patient_id<-factor(replicate(length(gvcf), individual))
		all_vcf<-c(all_vcf,gvcf)}
	}
	all_vcf$REF<-as.factor(all_vcf$REF)
	all_vcf$ALT<-as.factor(unstrsplit(CharacterList(all_vcf$ALT), sep = ","))
	all_vcf$QUAL<-as.factor(all_vcf$QUAL)
	all_vcf$FILTER<-as.factor(all_vcf$FILTER)

	print('ok1')
	vvcf = VRanges(
  					seqnames = seqnames(all_vcf),
  					ranges = ranges(all_vcf),
  					ref = all_vcf$REF,
  					alt = all_vcf$ALT,
  					study = paste(class_type,as.character(all_vcf$patient_id),sep = '')
  					)
	idx_snv<- ref(vvcf) %in% DNA_BASES & alt(vvcf) %in% DNA_BASES
	vvcfR<-vvcf[idx_snv]
	chrome<-c('chr1','chr2','chr3', 'chr4','chr5', 'chr6','chr7', 'chr8','chr9', 'chr10','chr11', 'chr12','chr13', 'chr14','chr15', 'chr16','chr17', 'chr18','chr19', 'chr20','chr21', 'chr22','chrX', 'chrY','chrM' )                    
	idx_snv<- as.character(seqnames(vvcfR)) %in% chrome
	
	vvcfR2<-vvcfR[idx_snv]

	sca_motifs = mutationContext(vvcfR2,k=3,BSgenome.Hsapiens.UCSC.hg19,strand=FALSE,unify=TRUE,check=TRUE)
	print('ok8')
	savefile<-paste(Analysis,'/sca_motifs_class_',class_type,'.rda',sep='')
	saveRDS(sca_motifs,savefile)

	count_person<- as.numeric(as.matrix(table(sca_motifs$study)))
	sca_mm = motifMatrix(sca_motifs, group = "study",normalize = TRUE)
	sca_mm_count<-t(t(sca_mm)*count_person)
	
	savefile<-paste(Analysis,'/matrix_class_',class_type,'_.txt',sep='')
	write.table(sca_mm_count,savefile,sep='\t')
	savefile<-paste(Analysis,'/NMF_MutationSpectrum_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 20, height = 120)
	print(dim(sca_motifs))
	tempplot<-plotMutationSpectrum(sca_motifs, "study")
	tempplot<-tempplot+ scale_fill_manual(values = rep("darkred", ncol(sca_mm)))
	print(tempplot)
	dev.off()
	print('ok7')
	opp <- opp_file[rep(seq_len(nrow(opp_file)), ncol(sca_mm_count)),]
	mut<-t(sca_mm_count)

	max_sig<-round(ncol(sca_mm),digits=0)-1
	if (max_sig>50){
		max_sig=20
	}

	print('max_sig')
	print(max_sig)
	signatures <- signeR(M=mut, Opport=opp, nlim=c(2,max_sig))
	BIC_score<-lapply(signatures$Test_BICs,mean)
	
	n_signature<-signatures$tested_n[which.max(BIC_score)]
	
	savefile<-paste(Analysis,'/EMu_BICboxplot_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-BICboxplot(signatures)
	print(tempplot)
	dev.off()
	
	n_sigs = 2:max_sig
	gof_nmf = assessNumberSignatures(sca_mm, n_sigs, nReplicates = 5)
	savefile<-paste(Analysis,'/NMF_NumberSignatures_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-plotNumberSignatures(gof_nmf)
	print(tempplot)
	dev.off()

	sigs_nmf = identifySignatures(sca_mm, n_signature, nmfDecomposition)
	write.table(samples(sigs_nmf),'/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step14data/wMatrix.csv',sep='\t')
	w_sig <- signatures(sigs_nmf)
	w_norm <- t(t(w_sig) / colSums(w_sig))   #check column
	results<-c()
	for (kg in 1:dim(w_norm)[2]){
    result<-c()
    for (line in 1:dim(X30_signatures)[2]){
      number<-KLD(unname(w_norm[,kg]),unname(unlist(X30_signatures[,line])))$mean.sum.KLD
      result<-c(result,number)
    }
    results<-rbind.data.frame(results,result)
	}
	colnames(results)<-colnames(X30_signatures)
	rownames(results)<-colnames(w_norm)
	savefile=paste(Analysis,'/NMF_sig_validate',class_type,'.csv',sep='')
	write.table(results, file = savefile, sep = ",", qmethod = "double")	

	savefile=paste(Analysis,'/NMF_data_sig_',class_type,'.rda',sep='')
	saveRDS(sigs_nmf,savefile)
	savefile=paste(Analysis,'/EMu_data_sig_',class_type,'.rda',sep='')
	saveRDS(signatures,savefile)

	savefile<-paste(Analysis,'/NMF_SignatureMap_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-plotSignatureMap(sigs_nmf) + ggtitle(" Signatures: NMF - Heatmap")
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/EMu_SignatureMap_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-SignHeat(signatures$SignExposures) 
	print(tempplot)
	dev.off()

	savefile<-paste(Analysis,'/NMF_Signatures_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 7)
	tempplot<-plotSignatures(sigs_nmf) + ggtitle(" Signatures: NMF - Barchart")+ylim(0,1)
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/EMu_Signatures_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 7)
	tempplot<-SignPlot(signatures$SignExposures)
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/NMF_ObservedSpectrum_class_',class_type,'.pdf',sep='')
	pdf(savefile,width = 20, height = 120)
	tempplot<-plotObservedSpectrum(sigs_nmf)
	tempplot<-tempplot+ scale_fill_manual(values = rep("darkred", ncol(sca_mm)))
	print(tempplot)
	dev.off()

	savefile<-paste(Analysis,'/NMF_FittedSpectrum_class_',class_type,'.pdf',sep='')
	pdf(savefile,width = 20, height = 120)
	tempplot<-plotFittedSpectrum(sigs_nmf)
	tempplot<-tempplot+ scale_fill_manual(values = rep("darkred", ncol(sca_mm)))
	print(tempplot)
	dev.off()

	savefile<-paste(Analysis,'/NMF_SampleMap_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 20, height = 80)
	tempplot<-plotSampleMap(sigs_nmf)
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/EMu_SampleMap_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 20, height = 80)
	tempplot<-ExposureHeat(signatures$SignExposures)
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/EMu_SampleExposure_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-ExposureBoxplot(signatures$SignExposures)
	print(tempplot)
	dev.off()
	print(class_type)
	savefile<-paste(Analysis,'/NMF_Samples_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 30, height = 10)
	tempplot<-SomaticSignatures::plotSamples(sigs_nmf)
	print(tempplot)
	dev.off()
	
	savefile<-paste(Analysis,'/EMu_Samples_class_',class_type,'_.pdf',sep='')
	pdf(savefile,width = 30, height = 10)
	tempplot<-ExposureBarplot(signatures$SignExposures)
	dev.off()

	clu_motif = clusterSpectrum(sca_mm, "motif")

	savefile<-paste(Analysis,'/NMF_ggdendrogram_class_',class_type,'.pdf',sep='')
	pdf(savefile,width = 10, height = 20)
	tempplot<-ggdendrogram(clu_motif, rotate = TRUE)
	print(tempplot)
	dev.off()
	print(class_type)
	print('end')
}
