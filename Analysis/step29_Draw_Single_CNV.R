library(readr)
library(plyr)
library("GenVisR")
library("maftools")


Microarry_file<-'/projects/p30007/Zexian/Alignment/MicroArray/ByPatient/'
files<-list.files(Microarry_file)


bandName<-'1p36.33 '
positions<-c('chr1',0,10000000)
pos_Inclu<-c(positions[1],as.integer(positions[2])-3*(as.integer(positions[3])-as.integer(positions[2])),as.integer(positions[3])+3*(as.integer(positions[3])-as.integer(positions[2])))


for(file in files){
	read_file<-paste(Microarry_file,file,sep='')
	print(read_file)
	Mic<-read.csv(read_file,sep='\t',header=TRUE)
	Mic_filter<-Mic[which(Mic$GC.Score>0.15),c('Chr','Position','CNV.Value')]
	Mic_fil_probe<-Mic_filter[ which(as.character(Mic_filter$Chr)==gsub('chr','',pos_Inclu[1] ) & Mic_filter$Position>as.integer(pos_Inclu[2]) & Mic_filter$Position<as.integer(pos_Inclu[3]) ),]
	Mic_fil_probe$Chr<-paste('chr',Mic_fil_probe$Chr,sep='')
	colnames(Mic_fil_probe)<- c("chromosome", "coordinate", "cn")
	Mic_fil_probe$chromosome<-as.factor(Mic_fil_probe$chromosome)

	MAF = '/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step29data/'
	outputFileName<-paste(MAF,gsub('.txt','',file),'_',bandName,'_.tif',sep='')

	tiff(outputFileName,width = 1500, height = 1200, units = 'px', res = 180)
	cnView(Mic_fil_probe, chr = positions[1], genome = "hg19", ideogram_txtSize = 4)
  	dev.off()
}



  jpeg(paste(objects[i], ".jpg", sep=""))
    plot(get(objects[i]))
  dev.off()






MAF = '/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step28data/merged.maf_format.txt'

#
data <- brcaMAF_f[brcaMAF_f$Hugo_Symbol == "CTNNA2", c("Hugo_Symbol", "HGVSp_Short","Group")]
data <- as.data.frame(cbind(data, "ENST00000466387"))
colnames(data) <- c("gene", "amino_acid_change", "Group","transcript_name")
data$Group<-ifelse(data$Group==1,'Case','Control')
tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step27data/CTNNA2.tif', width = 1500, height = 900, units = 'px', res = 180)
lolliplot(data,fillCol = "Group")
dev.off()























