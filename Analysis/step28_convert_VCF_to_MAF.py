# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:54:02 2016

@author: zexian
"""

import os
from os import listdir
import commands
#Might need rename files first: for f in *_1.fq; do mv $f $(echo ${f} |sed  's/_1.fq/_R1_.fq/'); done
direc=InputDir = '/projects/p30007/Zexian/Alignment/BBCAR_NEW'


picardtool='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/picard-tools-1.131/picard.jar'
GATKtool='java -jar /projects/p30007/Zexian/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
hg19Reference='/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta'
gold1000Indel ='/projects/p30007/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp='/projects/p30007/Zexian/reference/hg19/dbsnp_138.hg19.vcf'
annovar_call='perl /projects/p30007/Zexian/tools/annovar/table_annovar.pl'
snpeff_call='java -Xmx32g -jar /projects/p30007/Zexian/tools/snpEff/snpEff.jar'
takeExon='python /projects/p30007/Zexian/tools/DNAtools/Take_exon_splicing_from_annotated_VCF.py'
HaloPlex_for_Haplotype='python /projects/p30007/Zexian/tools/DNAtools/HaloPlex_for_Haplotype.py'
varScan_call='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/Varscan/VarScan.v2.3.9.jar'
varDict_call='/projects/p30007/Zexian/tools/VarDictJava/build/install/VarDict/bin/VarDict'
interval='/projects/p30007/Zexian/tools/DNAtools/S07604514_Padded.bed'
normal_pon='/projects/p30007/Zexian/Alignment/Germline_37/administrative/MuTect2_PON.vcf'


Raw_reads='/projects/p30007/Zexian/Alignment/BBCAR/RAW_data'
RootFolder = direc+'/WES_Analysis'
ScriptFolder = RootFolder+'/Scripts'
ScriptFolder=InputDir+'/administrative/Step2codes/'
BamFolder = RootFolder+'/BAM'
BamFolder = '/projects/b1042/lyglab/Zexian/BBCAR/BAM'
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

Haplo_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/Haplo_annovar_AnnoVCF'
Mutect_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/Mutect_annovar_AnnoVCF'
VarScan_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/VarScan_annovar_AnnoVCF'
VarDict_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/VarDict_annovar_AnnoVCF'

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

ScriptFolder=InputDir+'/administrative/Step28code/'
Raw_reads='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
OutF='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/MAF/'
listdirs = [d for d in os.listdir(Raw_reads)]
for file in listdirs:
    SampleName=file.replace('.hg19_multianno.vcf','')
    print(SampleName)
    script = open(ScriptFolder+'/'+SampleName+'.sh', "w")
    #Print MSUB Header
    script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomicsburst
#MSUB -l walltime=3:00:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=1
module load samtools
module load perl/5.22

''')    
    inputVCF=Raw_reads+file
    outputVCF=OutF+SampleName+'.maf'
    print >> script, 'perl /projects/p30007/Zexian/tools/vcf2maf/vcf2maf.pl --input-vcf '+inputVCF+' --output-maf '+outputVCF+' --tumor-id '+SampleName+'\n'
    script.close()

#for vep annotation 
#/usr/bin/perl /projects/p30007/Zexian/tools/VEP_data/vep/variant_effect_predictor.pl --species homo_sapiens --assembly GRCh37 --offline --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /projects/p30007/Zexian/tools/VEP_data/.vep --fasta /projects/p30007/Zexian/tools/VEP_data/.vep/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --format vcf --input_file /projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/1.hg19_multianno.vcf --output_file /projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/1.hg19_multianno.vep.vcf --fork 4 --check_allele --polyphen b --gmaf --maf_1kg --maf_esp --regulatory