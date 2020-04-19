# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:34:16 2016

@author: zza847
"""
import csv
from os import listdir
from os.path import isfile, join
import re
from os import listdir
import sys
import commands
import glob
import os

Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step27data/'

onlyfiles = [f for f in listdir(Input_folder) if isfile(join(Input_folder, f))]

with open (Output_folder+'allMutationsInOneFile.txt','w') as fout:
    fout.write('Tumor_Sample_Barcode\tVariant_Classification\tsignificance\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlt\tALE\tDepth\n')    
    for file in onlyfiles:
        if file.endswith('.vcf'):
            sample=file.replace('.hg19_multianno.vcf','')
            print(file)
            with open(Input_folder+file,'r') as fin:
                titleline=next(fin)
                for line in fin:
                    items=line.split('\t')
                    chro=items[0]
                    chro_position=items[1]
                    End_Position=str(int(chro_position)+len(items[4]))
                    Ref=items[3]
                    Alt=items[4]
                    ALE=items[-1].split(':')[2]
                    Depth=str(int(items[-1].split(':')[1].split(',')[0])+int(items[-1].split(':')[1].split(',')[1]))
                    infor=items[7]
                    ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                    ANN_Info_Items=ANN_Info.split(',')
                    ANN_Info_Item= ANN_Info_Items[0]
                    Annotations=ANN_Info_Item.split('|')
                    Variant_Classification = re.search('ExonicFunc.refGene=(.+?);AA', line).group(1)
                    Varianttype2=re.search('AAChange.refGene=(.+?);cytoBand', line).group(1)
                    significance= Annotations[2]
                    Hugo_Symbol=Annotations[3]
                    if Variant_Classification =='frameshift_deletion':
                        Variant_Classification='Frame_Shift_Del'
                    if Variant_Classification =='frameshift_insertion':
                        Variant_Classification='Frame_Shift_Ins'
                    if Variant_Classification =='nonframeshift_deletion':
                        Variant_Classification='In_Frame_Del'
                    if Variant_Classification =='nonframeshift_insertion':
                        Variant_Classification='In_Frame_Ins'
                    if Variant_Classification =='nonsynonymous_SNV':
                        Variant_Classification='Missense_Mutation'
                    if Variant_Classification =='stopgain' or Variant_Classification =='stoploss':
                        Variant_Classification='Nonsense_Mutation'
                    if Variant_Classification =='synonymous_SNV':
                        Variant_Classification='Silent'
                    if Variant_Classification =='.':
                        Variant_Classification='NA'   
                    if Variant_Classification =='unknown':
                        Variant_Classification='NA'   
                    result=sample+'\t'+Variant_Classification+'\t'+Varianttype2+'\t'+Hugo_Symbol+'\t'+chro+'\t'+chro_position+'\t'+End_Position+'\t'+Ref+'\t'+Alt+'\t'+ALE+'\t'+Depth+'\n'
                    fout.write(result)              
                                   
                    
                    