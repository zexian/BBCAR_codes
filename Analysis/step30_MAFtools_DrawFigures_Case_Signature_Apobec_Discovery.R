library(readr)
library(plyr)
library("GenVisR")
library("maftools")

#w matrix from home_made pipeline 
MAF = '/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step28data/merged.maf_format.txt'
brcaMA<-read.csv(MAF,sep='\t',header=TRUE)
brcaMA$i_TumorVAF_WU<-brcaMA$t_alt_count*100/(brcaMA$t_ref_count + brcaMA$t_alt_count )

Clinical_files='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt' 
clinical = read.table(Clinical_files,sep='\t',header=TRUE)
clinical$Basal<-ifelse(clinical$ER_Sta=='Negative'&clinical$PR_Sta=='Negative'&clinical$HER_Sta=='Negative',1,0)
clinical$CaseControl_num<-ifelse(clinical$CaseControl=='Case',1,0)
clinical<-clinical[,c('Study_ID','CaseControl_num','CaseControl','ER_Sta','PR_Sta','HER_Sta','Basal','interval_days')]

colnames(clinical)<-c('Tumor_Sample_Barcode','Type_num','Group','ER','PR','HER2','Basal','interval_days')

pos_sample<-clinical$Tumor_Sample_Barcode[clinical$Type_num==1  & clinical$interval_days>365*6]
neg_sample<-clinical$Tumor_Sample_Barcode[clinical$Type_num==0 & clinical$interval_days<365*3]


brcaMAF_Case<-brcaMA[brcaMA$Tumor_Sample_Barcode %in% pos_sample ,]
clinical_Case<-clinical[clinical$Tumor_Sample_Barcode %in% pos_sample,]

vcNames<-c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins")
vcNames<-as.character(unique(brcaMAF_Case$Variant_Classification))
maf=read.maf(brcaMAF_Case, clinicalData = clinical_Case,vc_nonSyn=vcNames)
maf_all=read.maf(brcaMA, clinicalData = clinical,vc_nonSyn=vcNames)



#start the plot 

rmOutlier = TRUE
dashboard = TRUE
showBarcodes = FALSE 
fontSize=fs = 1
titleSize = c(1, 0.8)
color = color=NULL
titv.color = titvColor=NULL
sfs = statFontSize =0.8
n = top=10
donut = pie=FALSE
rawcount = titvRaw= TRUE
stat = addStat='mean'
barcodes = showBarcodes=FALSE
barcodeSize = textSize =0.8



#start figure 1

col = c(RColorBrewer::brewer.pal(12, name = "Paired"), RColorBrewer::brewer.pal(11, 
        name = "Spectral")[1:3], "black", "violet", "royalblue")
col = grDevices::adjustcolor(col = col, alpha.f = alpha)
names(col) = names = c("Nonstop_Mutation", "Frame_Shift_Del", "IGR", "Missense_Mutation", "Silent", "Nonsense_Mutation","RNA", "Splice_Site", "Intron", "Frame_Shift_Ins", "Nonstop_Mutation","In_Frame_Del", "ITD", "In_Frame_Ins", "Translation_Start_Site","Multi_Hit", "Amp", "Del")
vcs = getSampleSummary(maf)
vcs = vcs[, colnames(vcs)[!colnames(x = vcs) %in% c("total", "Amp", "Del", "CNV_total")], with = FALSE]
vcs = vcs[, c(1, order(colSums(x = vcs[, 2:(ncol(vcs)), with = FALSE]),decreasing = TRUE) + 1), with = FALSE]
vcs.m = data.table::melt(data = vcs, id = "Tumor_Sample_Barcode")
colnames(vcs.m) = c("Tumor_Sample_Barcode", "Variant_Classification", "N")
data.table::setDF(vcs)
rownames(x = vcs) = vcs$Tumor_Sample_Barcode
vcs = vcs[, -1]
vcs = t(vcs)
lo = matrix(data = 1:6, nrow = 2, byrow = TRUE)


#extractSignatures
library('NMF')
#extract Case
matts<-trinucleotideMatrix(maf, '/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta', prefix = NULL, add = TRUE, useSyn = TRUE)
laml.sign = extractSignatures(mat = matts, n=3, plotBestFitRes = FALSE)
laml.sign<-readRDS('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/lamlsign_Case.rds')

tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signatures_CASE_Apobect.tif', width = 1600, height = 800, units = 'px', res = 160)
plotSignatures_modi(laml.sign, title_size = 0.8,show_title = TRUE,yaxisLim = 0.35 )
dev.off()

tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signatures_CASE_no_text.tif', width = 1600, height = 800, units = 'px', res = 160)
plotSignatures_modi(laml.sign, title_size = 0.8 ,show_title = FALSE,yaxisLim = 0.35)
dev.off()


#Signature enrichment
fileName<-'/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signature_Cluster_CASE.tif'
fn='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signature_Cluster_CASE_file_'
laml.se = signatureEnrichment(maf = maf, sig_res = laml.sign,fn =fn,SaveFig=fileName)
tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/EnrichmentResults_CASE.tif', width = 800, height = 800, units = 'px', res = 160)
plotEnrichmentResults(enrich_res = laml.se, pVal = 0.05)
dev.off()


#extract the all
laml.sign<-readRDS('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/lamlsign_All.rds')

tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signatures_ALL.tif', width = 1600, height = 800, units = 'px', res = 160)
plotSignatures_modi(laml.sign, title_size = 0.8,show_title = TRUE,yaxisLim = 0.35 )
dev.off()

tiff('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Signatures_ALL_no_text.tif', width = 1600, height = 800, units = 'px', res = 160)
plotSignatures_modi(laml.sign, title_size = 0.8, show_title = FALSE,yaxisLim = 0.35)
dev.off()




plotSignatures_modi<-function(nmfRes = NULL, contributions = FALSE, color = NULL, 
    patient_order = NULL, font_size = 1.2, show_title = NULL, 
    axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = NULL) {
    conv.mat.nmf.signatures = nmfRes$signatures
    contrib = nmfRes$contributions
    coSineMat = nmfRes$coSineSimMat
    if (contributions) {
        contribt = t(contrib)
        if (!is.null(patient_order)) {
            contribt = contribt[patient_order, ]
        }
        else {
            contribt = contribt[order(contribt[, ncol(contribt)]), 
                ]
        }
        contrib = t(contribt[, 1:(ncol(contribt))])
        cols = RColorBrewer::brewer.pal(n = 8, name = "Set2")
        if (show_barcodes) {
            lo = layout(mat = matrix(data = c(1, 2), nrow = 2), 
                heights = c(6, 2))
            par(mar = c(6, 4, 2, 1))
            b = barplot(contrib, axes = FALSE, horiz = FALSE, 
                col = cols, border = NA, names.arg = rep("", 
                  ncol(contrib)))
            axis(side = 1, at = b, labels = colnames(contrib), 
                lwd = 2, cex.axis = font_size, las = 2, line = 0.2, 
                hadj = 0.8, font = 2, tick = FALSE)
            axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 2, 
                las = 2, cex.axis = 0.9)
            mtext(text = "Signature exposures", side = 2, font = 2, 
                cex = 1, line = 2.8)
            plot.new()
            par(mar = c(2, 3, 0, 0))
            legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)], 
                border = NA, bty = "n", pch = 15, xpd = TRUE, 
                ncol = 1, cex = 1.2, pt.cex = 1.5, horiz = TRUE)
        }
        else {
            lo = layout(mat = matrix(data = c(1, 2), nrow = 2), 
                heights = c(6, 2))
            par(mar = c(3, 4, 2, 1))
            b = barplot(contrib, axes = FALSE, horiz = FALSE, 
                col = cols, border = NA, names.arg = rep("", 
                  ncol(contrib)))
            axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 2, 
                las = 2, cex.axis = 0.9)
            mtext(text = "Signature exposure", side = 2, font = 2, 
                cex = 1, line = 2.8)
            plot.new()
            par(mar = c(2, 3, 0, 0))
            legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)], 
                border = NA, bty = "n", pch = 15, xpd = TRUE, 
                ncol = 1, cex = 1.2, pt.cex = 1.5, horiz = TRUE)
        }
    }
    else {
        aetiology = structure(list(aetiology = c("spontaneous deamination of 5-methylcytosine", 
            "APOBEC Cytidine Deaminase (C>T)", "defects in DNA-DSB repair by HR", 
            "exposure to tobacco (smoking) mutagens", "Unknown", 
            "defective DNA mismatch repair", "UV exposure", "Unknown", 
            "defects in polymerase-eta", "defects in polymerase POLE", 
            "exposure to alkylating agents", "Unknown", "APOBEC Cytidine Deaminase (C>G)", 
            "Unknown", "defective DNA mismatch repair", "Unknown", 
            "Unknown", "Unknown", "Unknown", "defective DNA mismatch repair", 
            "unknown", "exposure to aristolochic acid", "Unknown", 
            "exposures to aflatoxin", "Unknown", "defective DNA mismatch repair", 
            "Unknown", "Unknown", "exposure to tobacco (chewing) mutagens", 
            "Unknown")), .Names = "aetiology", row.names = c("Signature_1", 
            "Signature_2", "Signature_3", "Signature_4", "Signature_5", 
            "Signature_6", "Signature_7", "Signature_8", "Signature_9", 
            "Signature_10", "Signature_11", "Signature_12", "Signature_13", 
            "Signature_14", "Signature_15", "Signature_16", "Signature_17", 
            "Signature_18", "Signature_19", "Signature_20", "Signature_21", 
            "Signature_22", "Signature_23", "Signature_24", "Signature_25", 
            "Signature_26", "Signature_27", "Signature_28", "Signature_29", 
            "Signature_30"), class = "data.frame")
        plotData = as.data.frame(t(conv.mat.nmf.signatures))
        nsigs = nrow(plotData)
        if (is.null(color)) {
            color = c("coral4", "lightcyan4", "deeppink3", "lightsalmon1", 
                "forestgreen", "cornflowerblue")
        }
        colors = rep(color, each = 16)
        par(mfrow = c(nsigs, 1), oma = c(5, 4, 0, 0) + 0.1, mar = c(0, 
            0, 2.5, 0) + 0.1, las = 1, tcl = -0.25, font.main = 4, 
            xpd = NA)
        for (i in 1:nsigs) {
            ae.sig = names(which(coSineMat[i, ] == max(coSineMat[i, 
                ])))
            ae = as.character(aetiology[ae.sig, ])
            ae = paste0(ae.sig, " like; cosine-similarity: ", 
                round(max(coSineMat[i, ]), digits = 3), " \n Aetiology: ", 
                ae)
            d = as.matrix(plotData[i, ])
            if (is.na(yaxisLim)) {
                bh = ceiling(max(d, na.rm = TRUE) * 10)/10
            }
            else {
                bh = 0.3
            }
            barplot(d, xaxt = "n", yaxt = "n", col = colors, 
                beside = TRUE, ylim = c(-0.1, bh), cex.main = 1, 
                border = NA, font.axis = 2, font.lab = 2, adj = 0.25)
            if (show_title) {
                title(main = ae, cex.main = title_size, line = 0)
            }
            axis(side = 2, at = seq(0, bh, 0.1), pos = -2, las = 2, 
                lwd = axis_lwd, hadj = 1.1, font = 2, cex.axis = font_size)
            rect(xleft = seq(0, 192, 32), ybottom = -0.05, xright = 192, 
                ytop = -0.02, col = color, border = "gray70")
            if (i == nsigs) {
                text(labels = c("C>A", "C>G", "C>T", "T>A", "T>C", 
                  "T>G"), y = rep(-0.1, 6), x = seq(0, 192, 32)[2:7] - 
                  16, cex = font_size, font = 2, font.lab = 2, 
                  pos = 1.2)
            }
        }
    }
}


