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
import pandas as pd
import os
from os import listdir
import sys


Somatic_Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step34data/'


#take the germline matrix 
dirs=Somatic_Input_folder
patientMapping=dict()
allGenesWithLocationMap=dict()
allGenesMap=dict()
allVariants=dict()
genes=[]
books_number={}

allfiles = [f for f in listdir(dirs) if isfile(join(dirs, f))]


#load person file
for file in allfiles:
    if file.endswith('.hg19_multianno.vcf'):
        subject=file.split('.')[0]
        patientMapping[subject]=1
        #readFile(filename,subject)


#load dict for rnalist
geneDict=dict()
genelist=['POLB']
for gene in genelist:
    geneDict[gene]=0

for person in patientMapping:
    books_number[person]={}
    for geneName in geneDict:
        books_number[person][geneName]=0


with open (Output_folder+'mutation_somatic.txt','w') as fout:
    for sample in patientMapping:
        filedirec1=dirs+sample+'.hg19_multianno.vcf'
        for filedirec in [filedirec1]:
            with open (filedirec,'r') as fin:
                for line in fin:
                    if line[0]!='#':
                        items=line.split('\t')
                        GeneName=''
                        ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                        ANN_Info_Item=ANN_Info.split(',')[0]
                        GeneName=ANN_Info_Item.split('|')[3].split('&')[0]
                        if GeneName in geneDict:
                            genome1000 = re.search('1000g2015aug_all=(.+?);ExAC', line)
                            if genome1000 is not None:
                                genome1000=genome1000.group(1)
                            else:
                                genome1000=1
                            if genome1000=='.':
                                genome1000=0
                            #if float(genome1000)<0.05:
                            if True:
                                fout.write(sample+'\t'+GeneName+'\t'+line)
                                if GeneName in books_number[sample]:
                                    books_number[sample][GeneName]+=1
                                else:
                                    books_number[sample][GeneName]=1
                                
df = pd.DataFrame(books_number).T.fillna(0)
if len(df.columns)>1:
    df.to_csv(Output_folder+'number_Somatic.csv', line_terminator='\n')
    

Germline_Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Haplotype/VCF_Anno/'
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step34data/'

#take the germline matrix 
dirs=Germline_Input_folder
patientMapping=dict()
allGenesWithLocationMap=dict()
allGenesMap=dict()
allVariants=dict()
genes=[]
books_number={}

allfiles = [f for f in listdir(dirs) if isfile(join(dirs, f))]


#load person file
for file in allfiles:
    if file.endswith('_filtered_snps.hg19_multianno.vcf'):
        subject=file.split('_')[0]
        patientMapping[subject]=1
        #readFile(filename,subject)


#load dict for rnalist
geneDict=dict()
genelist=['POLB']
for gene in genelist:
    geneDict[gene]=0

for person in patientMapping:
    books_number[person]={}
    for geneName in geneDict:
        books_number[person][geneName]=0


with open (Output_folder+'mutation_germline.txt','w') as fout:
    for sample in patientMapping:
        filedirec1=dirs+sample+'_filtered_snps.hg19_multianno.vcf'
        filedirec2=dirs+sample+'_filtered_indel.hg19_multianno.vcf'
        for filedirec in [filedirec1,filedirec2]:
            with open (filedirec,'r') as fin:
                for line in fin:
                    if line[0]!='#':
                        items=line.split('\t')
                        GeneName=''
                        ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                        ANN_Info_Item=ANN_Info.split(',')[0]
                        GeneName=ANN_Info_Item.split('|')[3].split('&')[0]
                        if GeneName in geneDict:
                            genome1000 = re.search('1000g2015aug_all=(.+?);ExAC', line)
                            if genome1000 is not None:
                                genome1000=genome1000.group(1)
                            else:
                                genome1000=1
                            if genome1000=='.':
                                genome1000=0
                            #if float(genome1000)<0.05:
                            if True:
                                fout.write(sample+'\t'+GeneName+'\t'+line)
                                if GeneName in books_number[sample]:
                                    books_number[sample][GeneName]+=1
                                else:
                                    books_number[sample][GeneName]=1
                                
df = pd.DataFrame(books_number).T.fillna(0)
if len(df.columns)>1:
    df.to_csv(Output_folder+'number_Germline.csv', line_terminator='\n')
    


        