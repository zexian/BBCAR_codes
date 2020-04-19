require(stringr)
library(ggplot2)
library(ggrepel)


direc='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step25data/Array_WES_Cross.txt'
OutFigure='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step26Figures/'

data_in = read.table(direc,sep='\t',header=TRUE)
#data_in<-data_in[data_in$sampleID!=1085,]
#data_in<-data_in[data_in$WES_Depth>90,]
#data_in<-data_in[(data_in$WES_Alt_frequency<0.05 & data_in$WES_Depth>100) |data_in$WES_Alt_frequency>0.05 ,]

data_in<-data_in[data_in$GCScore>0.15,]

WES<-data_in$WES_Alt_frequency
Array<-data_in$Array_Alt_frequency
cor(WES,Array)
data_in$diff<-abs(data_in$WES_Alt_frequency-data_in$Array_Alt_raw_frequency)

dim(data_in)

R2 <- 1 - (sum((WES-Array )^2)/sum((WES-mean(Array))^2))


  tiff(paste(OutFigure,'ALF_ScatterPlot.tif',sep=''), width = 1000, height = 1000, units = 'px', res = 180)
  p<- ggplot(data_in, aes(x=Array_Alt_raw_frequency, y=WES_Alt_frequency)) +
  geom_point(size=1, shape=18)+
  theme_classic() +
  xlim(0, 1)+
  ylim(0, 1)+ geom_abline(intercept = 0, slope = 1,linetype="dashed")+ 
  labs(title = "Comparison of Array and WES (MuTtect2)", x = "Array ALF", y = "WES ALF")
  #theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

results<-c()
for (num in c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)){
  data_in$group<-ifelse(data_in$diff>num,'Different Call','Same Call')
  re<-table(data_in$group)[2]/dim(data_in)[1]*100
  result<-c(num,re)
  results<-rbind(results,result)
}
results<-data.frame(results)
print(results)
