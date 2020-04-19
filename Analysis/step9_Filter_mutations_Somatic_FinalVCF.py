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
from shutil import copyfile
import os.path

Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
import os
import glob

files = glob.glob(Output_folder+'*')
for f in files:
    #print(f)
    os.remove(f)



Select_Dict=dict()
selected_file='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step8data/test.txt'
with open(selected_file,'r') as fin:
    title=next(fin)
    for line in fin:
        items=line.split('\t')
        key=items[0]+'_'+items[1]+'_'+items[2]
        Select_Dict[key.replace('"','')]=1


germline37='/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Mutect/PASS_VCF_Anno_Exonic/'
onlyfiles = [f for f in listdir(germline37) if isfile(join(germline37, f))]
for file in onlyfiles:
    if file.endswith('.vcf'):
        with open (germline37+file,'r') as fin, open (Output_folder+file,'w') as fout:
            for line in fin:
                fout.write('\t'.join(line.split('\t')[:-1])+'\n')






Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic/'
onlyfiles = [f for f in listdir(Input_folder) if isfile(join(Input_folder, f))]
for file in onlyfiles:
    if file.endswith('.vcf'):
        #print(file)
        if not os.path.isfile(Output_folder+file):
            #print(file)
            sample=file.replace('.hg19_multianno.vcf','')
            print(sample)
            with open(Input_folder+file,'r') as fin, open (Output_folder+file,'w') as fout:
                for line in fin:
                    if line[0:6]=='#CHROM':
                        fout.write(line)
                    if(line[0]!='#'):
                        items=line.split('\t')
                        key=sample+'_'+items[0]+'_'+items[1]
                        print(key)
                        if key in Select_Dict:
                            fout.write(line)







