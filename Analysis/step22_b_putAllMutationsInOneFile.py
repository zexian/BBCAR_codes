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
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step22data/'

onlyfiles = [f for f in listdir(Input_folder) if isfile(join(Input_folder, f))]

with open (Output_folder+'allMutationsInOneFile.txt','w') as fout:
    fout.write('sample\tvariant_type\tsignificance\tGene\tchro\tchro_position\tRef\tAlt\tALE\tDepth\n')    
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
                    Ref=items[3]
                    Alt=items[4]
                    ALE=items[-1].split(':')[2]
                    Depth=str(int(items[-1].split(':')[1].split(',')[0])+int(items[-1].split(':')[1].split(',')[1]))
                    infor=items[7]
                    ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                    ANN_Info_Items=ANN_Info.split(',')
                    ANN_Info_Item= ANN_Info_Items[0]
                    Annotations=ANN_Info_Item.split('|')
                    variant_type=Annotations[1]
                    significance= Annotations[2]
                    Gene=Annotations[3]
                    result=sample+'\t'+variant_type+'\t'+significance+'\t'+Gene+'\t'+chro+'\t'+chro_position+'\t'+Ref+'\t'+Alt+'\t'+ALE+'\t'+Depth+'\n'
                    fout.write(result)              
                                   
                    
                    