# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:07:20 2016

@author: zza847
"""

import numpy
import copy
from os import listdir
from os.path import isfile, join
import re
from os import listdir
import os
import csv
import os
from os import listdir
import sys



Benign='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
Cancer='/projects/p30007/Zexian/Alignment/Cancer_10/WES_Analysis/Mutect/VCF_Anno_Exonic/'

Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step31data/'

allfiles = [f for f in listdir(Cancer) if isfile(join(Cancer, f))]

with open (Output_folder+'cross_all.txt','w') as fout:
    for file in allfiles:
        if file.endswith('.vcf'):
            sample=file.replace('.hg19_multianno.vcf','')
            benign_file=Benign+sample+'.hg19_multianno.vep.vcf'
            cancer_file=Cancer+sample+'.hg19_multianno.vcf'
            
            This_person=dict()
            with open (cancer_file,'r') as fin:
                for line in fin:
                    if(line[0]!='#'):
                        items=line.split('\t')
                        key=items[0]+items[1]+items[3]+items[4]
                        if key not in This_person:
                            This_person[key]=line
            with open (benign_file,'r') as fin:
                for line in fin:
                    if(line[0]!='#'):
                        items=line.split('\t')
                        key=items[0]+items[1]+items[3]+items[4]
                        if key in This_person:
                            fout.write(sample+'\t'+This_person[key].strip('\n')+'\t'+This_person[key].split('\t')[-1].split(':')[2]+'\t'+sample+'\t'+line.strip('\n')+'\t'+line.split('\t')[-1].split(':')[2]+'\n')


