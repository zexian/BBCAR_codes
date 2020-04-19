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
ScriptFolder = InputDir+'/administrative/Step13codes/'
DataFolder = InputDir+'/administrative/Step13data/'

OutputFolder_normal = InputDir+'/administrative/Step2Data/'


BamFolder = RootFolder+'/BAM'
BamFolder = '/projects/b1042/lyglab/Zexian/Germline_37/BAM'
SampleDepthDirec=RootFolder+'/Analysis/Sample_Depth/'

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

Varscan_CNV=direc+'CNV/VarScan'
Varscan_CNV_mileup='/projects/b1042/lyglab/Zexian/Germline_37/mileup'

Validated_Samples=dict()
with open (InputDir+'/administrative/ValidaSamples.csv','r') as fin:
    title=next(fin)
    for line in fin:
        print(line.strip())
        Validated_Samples[line.strip()]=1

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
        Vscan_Sample=Varscan_CNV+'/'+SampleName
        #if not os.path.exists(Vscan_Sample):
        #    os.makedirs(Vscan_Sample)
        script = open(ScriptFolder+'/'+SampleName+'_GetReadDepth.sh', "w")
        #Print MSUB Header
        script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=39:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=2
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
        print >>script, 'samtools view -F 0x4 '+recalBam_normal+' | cut -f 1 | sort | uniq | wc -l > '+DataFolder+SampleName+'_normal_depth.txt\n'
        print >>script, 'samtools view -F 0x4 '+recalBam_blood+' | cut -f 1 | sort | uniq | wc -l > '+DataFolder+SampleName+'_blood_depth.txt\n'
        script.close()




                