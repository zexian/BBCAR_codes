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
import statsmodels.api as sm

patientGroup=dict()
patientMenaposal=dict()
patientAge=dict()

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

            if items[2]=='Postmenopausal':
                menaposal=1
            if items[2]=='Premenopausal':
                menaposal=0
            patientMenaposal[subject]=menaposal

            patientAge[subject]=age

books={}
for person in patientGroup:
    books[person]={}
    books[person]['group']=int(patientGroup[person])
    #books[person]['menaposal']=int(patientMenaposal[person])
    #books[person]['age']=int(patientAge[person])


InDirec='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Vep/'
for sample in patientGroup:
    filedirec=InDirec+sample+'.hg19_multianno.vep.vcf'
    with open (filedirec,'r') as fin:
        person=sample
        for line in fin:
            if(line[0]!='#'):
                if 'MODERATE' in line or 'HIGH' in line: 
                #if True:
                    items=line.split('\t')
                    Pop_Fre = re.search('ExAC_NFE=(.+?);ExAC_OTH', line)
                    if Pop_Fre is not None:
                        Pop_Fre=Pop_Fre.group(1).split(',')[0]
                    if Pop_Fre == '.':
                        Pop_Fre=0
                    print(Pop_Fre)
                    if True:
                    #if float(Pop_Fre)<=0.1:
                        GeneName = re.search('Interpro_domain=(.+?);1000g2015aug_all', line)
                        if GeneName is not None:
                            GeneName=GeneName.group(1).split(',')[0]
                        keys=GeneName.split('\\x3b')
                        for key in keys:
                            print(key)
                            if key in books[person]:
                                books[person][key]+=1
                            else:
                                books[person][key]=1

df = pd.DataFrame(books).T.fillna(0)
df=df.loc[:, (df.sum() > 1) ]
df.to_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step37data/ProteinDomain_mutations_all.csv', line_terminator='\n')







