# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 22:24:19 2016

@author: zza847
"""


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

MicroArray='/projects/p30007/Zexian/Alignment/MicroArray/ByPatient_intervals/'
Mutect_pass='/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Mutect/PASS_VCF/'
Mutect_BBCAR='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'

OutFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step24data/'
ScriptFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step24codes/'

listdirs = [d for d in os.listdir(MicroArray)]
for sample_micoarry in listdirs:
	print(sample_micoarry)
	sample=sample_micoarry.replace('_1','')
	sample=sample.replace('.txt','')
	sample_mutect=sample+'.vcf'
	sample_bbcar=sample+'.hg19_multianno.vcf'
	if not os.path.exists(Mutect_pass+sample_mutect):
		if os.path.exists(Mutect_BBCAR+sample_bbcar):
			print('in')
			script = open(ScriptFolder+'/'+sample+'.sh', "w")
			script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=47:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=4
module load picard/1.131
module load bwa/0.7.12
module load python
module load samtools/0.1.18 
module load java/jdk1.8.0_25
module load  R/3.4.3

''')
			print >> script, 'python /projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/step24_Format_VCF_to_Microarray_Mutect2.py '+MicroArray+sample_micoarry+' '+Mutect_BBCAR+sample_bbcar+' '+OutFolder+sample+'.txt'+'\n'
			script.close()

            
            
            
            
                