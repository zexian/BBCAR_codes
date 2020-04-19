# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import collections
import csv
from os import listdir
from os.path import isfile, join
import re
from os import listdir
import sys
import pandas as pd
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
from multiprocessing import Pool
import multiprocessing
from os.path import join
from os import listdir

patientGroup=dict()
patientMenaposal=dict()
patientAge=dict()


GeneDict=dict()
for gene in ['MHY11']:
    GeneDict[gene]=1    

with open ('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt','r') as fin:
    title=next(fin)
    for line in fin:
        items=line.split('\t')
        subject=items[0]
        group=''
        menaposal=''
        #print(items[3])
        age=int(items[3])
        if subject not in patientGroup:
            if items[1]=='Case':
                group=1
            if items[1]=='Control':
                group=0
            patientGroup[subject]=group

    

books={}
for person in patientGroup:
    books[person]={}
    books[person]['group']=int(patientGroup[person])
    #books[person]['menaposal']=int(patientMenaposal[person])
    #books[person]['age']=int(patientAge[person])


InDirec='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
for sample in patientGroup:
    filedirec=InDirec+sample+'.hg19_multianno.vcf'
    with open (filedirec,'r') as fin:
        person=sample
        for line in fin:
            if(line[0]!='#'):
                #if 'MODERATE' in line or 'HIGH' in line: 
                if True:
                    items=line.split('\t')
                    Pop_Fre = re.search('ExAC_NFE=(.+?);ExAC_OTH', line)
                    if Pop_Fre is not None:
                        Pop_Fre=Pop_Fre.group(1).split(',')[0]
                    if Pop_Fre == '.':
                        Pop_Fre=0
                    #print(Pop_Fre)
                    if True:
                    #if float(Pop_Fre)<=0.1:
                        GeneName = re.search('Gene.refGene=(.+?);GeneDetail', line)
                        if GeneName is not None:
                            GeneName=GeneName.group(1).split(',')[0]
                        genes=GeneName.split('\\x3b')
                        for gene in genes:
                            if gene in GeneDict:
                                print(sample)
                                key=items[0]+'_'+items[1]+items[3]+'_'+items[4]
                                key=gene
                                Freq=float(items[-1].split(':')[2])
                                print(key)
                                if key in books[person]:
                                    books[person][key]+=1
                                else:
                                    books[person][key]=1

df = pd.DataFrame(books).T.fillna(0)
df=df.loc[:, (df.sum() >= 1) ]
df.to_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step41data/MHY11.csv', line_terminator='\n')







