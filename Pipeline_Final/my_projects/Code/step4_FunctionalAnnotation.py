# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:54:02 2016

@author: zexian
"""

import os
from os import listdir
import commands
from os.path import isfile, join
import re
from os import listdir
import sys

#Might need rename files first: for f in *_1.fq; do mv $f $(echo ${f} |sed  's/_1.fq/_R1_.fq/'); done


annovar_call='perl /projects/p30007/Zexian/tools/annovar/table_annovar.pl'
snpeff_call='java -jar /projects/p30007/Zexian/tools/snpEff/snpEff.jar'
VEP_call ='/usr/bin/perl /projects/p30007/Zexian/tools/VEP_data/vep/variant_effect_predictor.pl'
ScriptFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step5codes/'
Tempfiles='/projects/p30007/Zexian/Alignment/BBCAR_NEW/Tempfiles/'
AnnoVCF='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno/'
FinalAnnoVCF='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_FinalAnno/'

InDirec='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/PASS_VCF/'
listdirs = [f for f in listdir(InDirec) if isfile(join(InDirec, f))]

	
for file in listdirs:
    SampleName=file.replace('.vcf','')
    script = open(ScriptFolder+SampleName+'_annovar.sh', "w")
    #Print MSUB Header
    script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=9:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=4
module load picard/1.131
module load bwa/0.7.12
module load samtools/1.2
module load python 
module load java/jdk1.8.0_25
module load  R/3.3.3

''')    

    print >> script, snpeff_call+' hg19 -ss 3 -canon '+InDirec+file+' > '+Tempfiles+'/'+SampleName+'.vcf\n'
    print >> script, annovar_call+' --thread 4 '+Tempfiles+'/'+SampleName+'.vcf'+' /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+AnnoVCF+SampleName+' -remove -protocol avsnp147,refGene,cytoBand,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp35a,revel,clinvar_20180603,snp138NonFlagged,cosmic80 -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
    print >> script, VEP_call+ ' --species homo_sapiens --assembly GRCh37 --offline --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /projects/p30007/Zexian/tools/VEP_data/.vep --fasta /projects/p30007/Zexian/tools/VEP_data/.vep/homo_sapiens/94_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --format vcf --input_file '+AnnoVCF+SampleName+' --output_file '+FinalAnnoVCF+SampleName+' --fork 4 --check_allele --polyphen b --gmaf --maf_1kg --maf_esp --regulatory'
    script.close()


