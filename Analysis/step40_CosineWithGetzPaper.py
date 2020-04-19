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

mutationType=dict()
with open ('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Hmatrix_Case.csv','r') as fin:
    title=next(fin)
    for line in fin:
        item=line.split(',')[0].strip('"')
        mutationType[item]=0

with open ('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step40data/Getz.csv','r') as fin:
    title=next(fin)
    for line in fin:
        items=line.split(',')
        firstone=items[0]
        number=float(items[3])
        left=firstone[3]
        right=firstone[4]
        fromM=firstone[0]
        toM=firstone[1]
        key=left+'['+fromM+'>'+toM+']'+right
        mutationType[key]=mutationType[key]+number

with open ('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step40data/Getz_summarized.csv','w') as fout:
    for key in mutationType:
        fout.write(key+','+str(mutationType[key])+'\n')


    

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
                    print(Pop_Fre)
                    if True:
                    #if float(Pop_Fre)<=0.1:
                        #GeneName = re.search('Interpro_domain=(.+?);1000g2015aug_all', line)
                        #if GeneName is not None:
                        #    GeneName=GeneName.group(1).split(',')[0]
                        #keys=GeneName.split('\\x3b')
                        #for key in keys:
                            #print(key)
                        if 'MutationNumber' in books[person]:
                            books[person]['MutationNumber']+=1
                        else:
                            books[person]['MutationNumber']=1

df = pd.DataFrame(books).T.fillna(0)
df=df.loc[:, (df.sum() > 1) ]
df.to_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step39data/MutationNumber.csv', line_terminator='\n')







