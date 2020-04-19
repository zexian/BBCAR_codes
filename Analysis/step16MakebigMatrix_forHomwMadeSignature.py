#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 09:41:53 2018

@author: zza847
"""


import csv
from os import listdir
from os.path import isfile, join
import re
from os import listdir
import sys
import os 
from os import listdir
from os.path import isfile, join

def Keyfunction(AA_Origi,Position,Direction,Ref,Alt_ori,five_AA,three_AA):
    if Direction=='reverse':
        if Alt_ori=='A':
            Alt_ori='T'
        elif Alt_ori=='T':
            Alt_ori='A'
        elif Alt_ori=='C':
            Alt_ori='G'
        elif Alt_ori=='G':
            Alt_ori='C'  
    if Position=='2':
        context=AA_Origi[0]+AA_Origi[2]
    elif Position=='1':
        context=five_AA[2]+AA_Origi[1]
    elif Position=='3':
        context=AA_Origi[1]+three_AA[0]
    makekey=Ref+Alt_ori+'_'+AA_Origi+'_'+Position+'_'+context
    return makekey

intputFo='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step15data/'
Out='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step16data/'
if not os.path.exists(Out):
    os.makedirs(Out)
                
index0Dict=dict()
index2Dict=dict()
index3Dict=dict()

onlyfiles = [f for f in listdir(intputFo) if isfile(join(intputFo, f))]
for file in onlyfiles:
    print(file)
    if file.endswith('hg19_multianno.csv'):
        inputfilename=intputFo+file

        with open (inputfilename,'r') as fin: 
            title=next(fin)
            for line in fin:
                items=line.split('\t')
                #makeKey=Keyfunction(items[9],items[11],items[7],items[8],items[6],items[14],items[16])
                makekey=items[8]+'_'+items[9]+'_'+items[10]+'_'+items[12]+'_'+items[13]+'_'+items[14].replace('\n','')
                makekey2=items[8]+'_'+items[9]+'_'+items[10]+'_'+items[13]+'_'+items[14].replace('\n','')
                makekey3=items[8]+'_'+items[9]+'_'+items[13]+'_'+items[14].replace('\n','')
                if makekey not in index0Dict:
                    index0Dict[makekey]=0
                if makekey2 not in index2Dict:
                    index2Dict[makekey2]=0
                if makekey3 not in index3Dict:
                    index3Dict[makekey3]=0                        
results=[]  
results2=[]                   
results3=[]
print(len(index0Dict))                   

onlyfiles = [f for f in listdir(intputFo) if isfile(join(intputFo, f))]
for file in onlyfiles:
    if file.endswith('.hg19_multianno.csv'):
        inputfilename=intputFo+file
        with open (inputfilename,'r') as fin: 
            currentDict=index0Dict.copy()
            currentDict2=index2Dict.copy()
            currentDict3=index3Dict.copy()
            title=next(fin)
            for line in fin:
                items=line.split('\t')
                #if items[18].strip('\n')=='no':
                makekey=items[8]+'_'+items[9]+'_'+items[10]+'_'+items[12]+'_'+items[13]+'_'+items[14].replace('\n','')
                makekey2=items[8]+'_'+items[9]+'_'+items[10]+'_'+items[13]+'_'+items[14].replace('\n','')
                makekey3=items[8]+'_'+items[9]+'_'+items[13]+'_'+items[14].replace('\n','')
                currentDict[makekey]=currentDict[makekey]+1
                currentDict2[makekey2]=currentDict2[makekey2]+1
                currentDict3[makekey3]=currentDict3[makekey3]+1

            result=[]
            result2=[]
            result3=[]
            result.append(file.replace('.hg19_multianno.csv',''))
            result2.append(file.replace('.hg19_multianno.csv',''))
            result3.append(file.replace('.hg19_multianno.csv',''))
            for key in sorted(currentDict):
                result.append(currentDict[key])
            results.append(result)
            for key in sorted(currentDict2):
                result2.append(currentDict2[key])
            results2.append(result2)
            for key in sorted(currentDict3):
                result3.append(currentDict3[key])
            results3.append(result3)                
with open (Out+'change_Codon_position_context.csv','w') as fout:
#with open (Out+'change_Codon_context.csv','w') as fout:
#with open (Out+'change_context.csv','w') as fout:
    keylist=''
    for key in sorted(currentDict):
        keylist=keylist+','+key
        print(repr(key))
    fout.write(keylist+'\n')
    for line in results:
        newline=str(line).replace('[','').replace(']','').replace('','')
        fout.write(newline+'\n')    
        
#with open (Out+'change_Codon_position_context.csv','w') as fout:
with open (Out+'change_Codon_context.csv','w') as fout:
#with open (Out+'change_context.csv','w') as fout:
    keylist=''
    for key in sorted(currentDict2):
        keylist=keylist+','+key
        print(repr(key))
    fout.write(keylist+'\n')
    for line in results2:
        newline=str(line).replace('[','').replace(']','').replace('','')
        fout.write(newline+'\n')      
    
#with open (Out+'change_Codon_position_context.csv','w') as fout:
#with open (Out+'change_Codon_context.csv','w') as fout:
with open (Out+'change_context.csv','w') as fout:
    keylist=''
    for key in sorted(currentDict3):
        keylist=keylist+','+key
        print(repr(key))
    fout.write(keylist+'\n')
    for line in results3:
        newline=str(line).replace('[','').replace(']','').replace('','')
        fout.write(newline+'\n')  
    