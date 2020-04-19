# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 22:24:19 2016

@author: zza847
"""
direc=InputDir='/projects/p30007/Zexian/Alignment/Germline_37'
direc_normal='/projects/p30007/Zexian/Alignment/BBCAR/RAW_data'


from os import listdir
from os.path import isfile, join
import csv
import os
import os.path
import commands

picardtool='java -jar /projects/p30007/Zexian/tools/picard-tools-1.131/picard.jar'
GATKtool='java -jar /projects/p30007/Zexian/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
hg19Reference='/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta'
gold1000Indel ='/projects/p30007/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp='/projects/p30007/Zexian/reference/hg19/dbsnp_138.hg19.vcf'
annovar_call='perl /projects/p30007/Zexian/tools/annovar/table_annovar.pl'
snpeff_call='java -Xmx32g -jar /projects/p30007/Zexian/tools/snpEff/snpEff.jar'
HaloPlex_for_Haplotype='python /projects/p30007/Zexian/tools/DNAtools/HaloPlex_for_Haplotype.py'
takeExon='python /projects/p30007/Zexian/tools/DNAtools/Take_exon_splicing_from_annotated_VCF.py'
varScan_call='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/Varscan/VarScan.v2.3.9.jar'
varDict_call='/projects/p30007/Zexian/tools/VarDictJava/build/install/VarDict/bin/VarDict'
interval='/projects/p30007/Zexian/tools/DNAtools/S07604514_Padded.bed'
normal_pon='/projects/p30007/Zexian/Alignment/Germline_37/administrative/MuTect2_PON.vcf'

Raw_reads=direc+'/RAW_data'
RootFolder = direc+'/WES_Analysis'
ScriptFolder = InputDir+'/administrative/Step14codes/'
depthFolder=InputDir+'/administrative/Step13data/'
OutputFolder_normal = InputDir+'/administrative/Step2Data/'


BamFolder = RootFolder+'/BAM'
BamFolder = '/projects/b1042/lyglab/Zexian/Germline_37/BAM'

HaplotypeFolder=RootFolder+'/Haplotype'
MutecFolder=RootFolder+'/Mutect'
VarScanFolder=RootFolder+'/VarScan'
VarDcitFolder=RootFolder+'/VarDict'

Haplo_VCF=HaplotypeFolder+'/VCF'
Mutect_VCF=MutecFolder+'/VCF'
VarScan_VCF=VarScanFolder+'/VCF'
VarDict_VCF=VarDcitFolder+'/VCF'

Filter_Haplo_VCF=HaplotypeFolder+'/Filter_VCF'


Haplo_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/Haplo_annovar_AnnoVCF'
Mutect_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/Mutect_annovar_AnnoVCF'
VarScan_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/VarScan_annovar_AnnoVCF'
VarDict_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/VarDict_annovar_AnnoVCF'

Haplo_FinalVCF=HaplotypeFolder+'/annotated_VCF'
Mutect_FinalVCF=MutecFolder+'/annotated_VCF'
VarScan_FinalVCF=VarScanFolder+'/annotated_VCF'
VarDict_FinalVCF=VarDcitFolder+'/annotated_VCF'

Haplo_Exon=HaplotypeFolder+'/exon_VCF'
Haplo_Exon_5per=HaplotypeFolder+'/exon_5per_VCF'
Haplo_Exon_1per=HaplotypeFolder+'/exon_1per_VCF'

Mutec_Exon=MutecFolder+'/exon_VCF'
Mutec_Exon_5per=MutecFolder+'/exon_5per_VCF'
Mutec_Exon_1per=MutecFolder+'/exon_1per_VCF'

VarScan_Exon=VarScanFolder+'/exon_VCF'
VarScan_Exon_5per=VarScanFolder+'/exon_5per_VCF'
VarScan_Exon_1per=VarScanFolder+'/exon_1per_VCF'

VarDict_Exon=VarDcitFolder+'/exon_VCF'
VarDict_Exon_5per=VarDcitFolder+'/exon_5per_VCF'
VarDict_Exon_1per=VarDcitFolder+'/exon_1per_VCF'

Analysis = RootFolder+'/Analysis'
Variant_plot=Analysis+'/Variant_plot'

Mutect_preprocess=MutecFolder+'/ISOWN/Preprocess'
Mutect_midprocess=MutecFolder+'/ISOWN/Midprocess'
VarScan_preprocess=VarScanFolder+'/ISOWN/Preprocess'
VarScan_midprocess=VarScanFolder+'/ISOWN/Midprocess'
VarDict_preprocess=VarDcitFolder+'/ISOWN/Preprocess'
VarDict_midprocess=VarDcitFolder+'/ISOWN/Midprocess'

Varscan_CNV=direc+'/CNV/VarScan/'
Varscan_CNV_Process=direc+'/CNV/VarScan_Process/'
Varscan_CNV_mileup='/projects/b1042/lyglab/Zexian/Germline_37/mileup'

Validated_Samples=dict()
with open (InputDir+'/administrative/ValidaSamples.csv','r') as fin:
    title=next(fin)
    for line in fin:
        print(line.split(',')[0].strip())
        Validated_Samples[line.split(',')[0].strip()]=1

filename_dict=dict()
bamfile_dict=dict()
cancerfilename_dict=dict()
cancerbamfile_dict=dict()
listdirs = [d for d in os.listdir(Raw_reads) if os.path.isdir(os.path.join(Raw_reads, d))]
FileList = []
Seqinfors=dict()
pattern= '_R1.fastq.gz'

for sample_folder in listdirs:
    if sample_folder in Validated_Samples:
        #print(sample_folder)
        SampleName=sample_folder

        Vscan_Sample=Varscan_CNV+SampleName+'/'
        print(Vscan_Sample)
        if not os.path.exists(Vscan_Sample):
            os.makedirs(Vscan_Sample)
        script = open(ScriptFolder+'/'+SampleName+'_Varscan.sh', "w")
        #Print MSUB Header
        script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=39:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=4
module load picard/1.131
module load bwa/0.7.12
module load python
module load samtools/0.1.18 
module load java/jdk1.8.0_25
module load  R

''')    
        bam_string_normal=''
        bam_string_blood=''
        bamSampleFolder=BamFolder+'/'+sample_folder
        recalBam_normal=bamSampleFolder+'/'+SampleName+'_normal_recal_reads.bam'
        recalBam_blood=bamSampleFolder+'/'+SampleName+'_blood_recal_reads.bam'
        #call copy number variation using varscan
        with open (depthFolder+'/'+SampleName+'_blood_depth.txt') as fin:
            blood_depth=next(fin)
        with open (depthFolder+'/'+SampleName+'_normal_depth.txt') as fin:
            normal_depth=next(fin)
        data_ratio=float(float(blood_depth)/float(normal_depth))
        Copynumber_copynumber=Varscan_CNV+'/'+SampleName+'.copynumber'
        Copynumber_called=Varscan_CNV+'/'+SampleName+'.copynumber.called'
        Copynumber_called_homdel=Varscan_CNV+'/'+SampleName+'.copynumber.called.homdel'
        Copynumber_called_recenter=Varscan_CNV+'/'+SampleName+'.copynumber.called.recenter'
        Copynumber_seg=Varscan_CNV+'/'+SampleName+'.seg.pvalue'
        Copynumber_seg_merged=Varscan_CNV+'/'+SampleName+'.seg.pvalue.mergeSegment'
        markerPosition_file=Varscan_CNV+'/'+SampleName+'.markerPosition'
        print >>script, 'samtools mpileup -q 1 -f '+hg19Reference+' ' +recalBam_blood+' ' +recalBam_normal+' |awk \'{if($4 >= 6) print $0}\' | awk \'{if($7 != 0) print $0}\' | '+varScan_call+' copynumber -mpileup '+ Varscan_CNV+SampleName+' --mpileup 1 --p-value 0.01 --min-coverage 20 --min-map-qual 20 --min-base-qual 20 --data-ratio ' +str(data_ratio)+'\n'
        print >>script, varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called +' --output-homdel-file '+Copynumber_called_homdel+'\n'
        print >>script, 'delta=$(python /projects/p30007/Zexian/tools/varscan2/meanLogRatioByChromosome.py '+Copynumber_called+')'
        print >>script, 'cmp=$(awk -v delta=$delta \'END{if (delta < -0.01) {print "lt"} else {if (delta > 0.01) {print "gt"} else {print "eq"}}}\' < /dev/null)'
        print >>script, 'if [[ "$cmp" == "lt" ]]; then'
        print >>script, '   rd=$(echo $delta | sed \'s/-//\')'
        print >>script, '   '+varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called_recenter +' --output-homdel-file '+Copynumber_called_homdel+' --recenter-down $rd'
        print >>script, 'elif [[ "$cmp" == "gt" ]]; then'
        print >>script, '    '+varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called_recenter +' --output-homdel-file '+Copynumber_called_homdel+' --recenter-up $delta'
        print >>script, 'else'
        print >>script, '    cp '+Copynumber_copynumber+' '+Copynumber_called_recenter
        print >>script, 'fi'
        print >>script, 'Rscript /projects/p30007/Zexian/tools/DNAtools/CopyNumberR.R '+ Copynumber_called_recenter+' ' + Copynumber_seg+' '+SampleName+' '+markerPosition_file+'\n'
        print >>script, 'perl /projects/p30007/Zexian/tools/DNAtools/mergeSegment.pl '+Copynumber_seg+' --ref-arm-sizes /projects/p30007/Zexian/tools/DNAtools/hg19.len.bed --output-basename '+Copynumber_seg_merged+'\n'
        print >>script, 'python /projects/p30007/Zexian/tools/DNAtools/Clean_CBS_arm_file3.py '+ Varscan_CNV+'/'+SampleName+'.seg.pvalue.mergeSegment.events.tsv ' + Varscan_CNV+'/'+SampleName+ '.segmentedFile.clean.txt '+ SampleName+'\n'
        script.close()

#make the code for repetitive CNV for Varscan
patient_dict=dict()
with open (InputDir+'/administrative/ValidaSamples.csv','r') as fin:
    f_line=next(fin)
    for line in fin:
        items=line.split(',')
        patient_dict[items[0]]=items[1].strip()
case_seg=''
case_mark=''
control_seg=''
control_mark=''
all_seg=''
all_mark=''
for patient in patient_dict:
    disease=patient_dict[patient]
    if disease=='Case' :
        case_seg=case_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
        case_mark=case_mark+Varscan_CNV+'/'+patient+'.markerPosition '
    if disease=='Control' :
        control_seg=control_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
        control_mark=control_mark+Varscan_CNV+'/'+patient+'.markerPosition '
    if (disease=='Control' or disease =='Case' ) :
        all_seg=all_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
        all_mark=all_mark+Varscan_CNV+'/'+patient+'.markerPosition '
CNV_visulization = open(InputDir+'/administrative/'+'step14b_runGistic2.sh', "w")
CNV_visulization.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=39:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=4
module load picard/1.131
module load bwa/0.7.12
module load python
module load samtools/0.1.18 
module load java/jdk1.8.0_25
module load  R

''')    
print >> CNV_visulization, 'cd /projects/p30007/Zexian/tools/GISTIC_2_0_23/'
#print >> signature_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/Clean_vcf_for_signature_paired.py '+ InputDir+'\n'
print >> CNV_visulization, 'cat '+case_mark+' > ' +Varscan_CNV_Process+'/Case_normal.clean.markerPosition\n'
print >> CNV_visulization, 'cat '+case_seg+' > ' +Varscan_CNV_Process+'/Case_segmentedFile.clean.txt\n'
print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV_Process+'/Case/ -seg '+Varscan_CNV_Process+'/Case_segmentedFile.clean.txt -mk '+Varscan_CNV_Process+'/Case_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+ Varscan_CNV_Process+'/Case/log.txt \n'

print >> CNV_visulization, 'cat '+control_mark+' > ' +Varscan_CNV_Process+'/Control_normal.clean.markerPosition\n'
print >> CNV_visulization, 'cat '+control_seg+' > ' +Varscan_CNV_Process+'/Control_segmentedFile.clean.txt\n'
print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV_Process+'/Control/ -seg '+Varscan_CNV_Process+'/Control_segmentedFile.clean.txt -mk '+Varscan_CNV_Process+'/Control_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+Varscan_CNV_Process+'/Control/log.txt \n'

print >> CNV_visulization, 'cat '+all_mark+' > ' +Varscan_CNV_Process+'/All_normal.clean.markerPosition\n'
print >> CNV_visulization, 'cat '+all_seg+' > ' +Varscan_CNV_Process+'/All_segmentedFile.clean.txt\n'
print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV_Process+'/All/ -seg '+Varscan_CNV_Process+'/All_segmentedFile.clean.txt -mk '+Varscan_CNV_Process+'/All_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+Varscan_CNV_Process+'/All/log.txt\n'
CNV_visulization.close()
            
            
           












            
            
            
            
                