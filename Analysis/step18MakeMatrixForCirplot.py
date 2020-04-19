#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  8 21:36:30 2018

@author: zza847
"""

import itertools
from itertools import combinations
import os
from os import listdir
import commands
from os.path import isfile, join
import re
import sys

signumber={1:'SigA',2:'SigB',3:'SigC',4:'SigD',5:'SigE',6:'SigF'}
a = list(">".join(items) for items in itertools.permutations(['A','G','T','C'],2))
aa=['C>A', 'C>G', 'C>T', 'T>A', 'T>G','T>C']
b = list("".join(items) for items in itertools.product(['A','G','T','C'],repeat=2))
c= list("".join(items) for items in itertools.product(['A','G','T','C'],repeat=3))



InDirec='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step17data/'

Folder=InDirec+'/change_Codon_position_context/'
allSubSignatures = [f for f in listdir(Folder) if os.path.isdir(os.path.join(Folder, f))]
for subsig in allSubSignatures:
    SignatureFolder=Folder+subsig+'/'
    CirplotData=SignatureFolder+'Data/'
    CirplotPlot=SignatureFolder+'Plot/'
    if not os.path.exists(CirplotData):
        os.makedirs(CirplotData)
    if not os.path.exists(CirplotPlot):
        os.makedirs(CirplotPlot)        
    #take the actual number, loop through the signatures 
    lengthOfsig=int(re.search(r'\d+', subsig).group())
    for signum in range(1,(lengthOfsig+1)):
        #define the dictionaries
        side1st=dict()
        side2nd=dict()
        side3rd=dict()
        sideall=dict()
        books1={}
        books2={}
        books3={}
        booksall={}
        #populate the dictornary keys 
        for key1 in a:
            for key2 in b:
                side1st[key1+'_'+key2]=side2nd[key1+'_'+key2]=side3rd[key1+'_'+key2]=sideall[key1+'_'+key2]=0
        #populate the dictornary keys 
        for types in c:
            books1[types]={}
            books2[types]={}
            books3[types]={}
            booksall[types]={}
        #populate the dictornary keys 
        for key1 in books1:
            for key2 in side1st:
                books1[key1][key2]=0
        for key1 in books2:
            for key2 in side1st:
                books2[key1][key2]=0
        for key1 in books3:
            for key2 in side1st:
                books3[key1][key2]=0
        for key1 in booksall:
            for key2 in side1st:
                booksall[key1][key2]=0
        #read in the data 
        for class_type in ['all','Case','Control']:
            with open(SignatureFolder+class_type+'_'+'h_matrix.csv','r') as fin:
                title=next(fin)
                for line in fin:
                    line=line.strip('\n')
                    items=line.split(',')
                    item0=items[0].strip('\"')
                            #sig number 
                    Sig=signumber[signum]
                    number=float(items[signum])
                    #key 1 is A>C_AA
                    #print(disea + '  '+subsig)
                    key1_in=item0[0]+'>'+item0[2]+'_'+item0.split('_')[-2]+item0.split('_')[-1]
                    #key2 is the code ACT
                    key2_in=item0.split('_')[2]
                    #accumulate by position informaiton
                    position=item0.split('_')[3]
                    if position =='1':
                        side1st[key1_in]+=number
                        books1[key2_in][key1_in]+=number
                    elif position =='2':
                        side2nd[key1_in]+=number
                        books2[key2_in][key1_in]+=number
                    elif position =='3':
                        side3rd[key1_in]+=number
                        books3[key2_in][key1_in]+=number
                    sideall[key1_in]+=number
                    booksall[key2_in][key1_in]+=number
                    #now save the matrix 
                with open(CirplotData+class_type+'_'+Sig+'_y_1.csv','w') as fout:
                    for key in sorted(side1st):
                        fout.write(str(key)+','+str(side1st[key])+'\n')
                with open(CirplotData+class_type+'_'+Sig+'_y_2.csv','w') as fout:
                    for key in sorted(side2nd):
                        fout.write(str(key)+','+str(side2nd[key])+'\n')
                with open(CirplotData+class_type+'_'+Sig+'_y_3.csv','w') as fout:
                    for key in sorted(side3rd):
                        fout.write(str(key)+','+str(side3rd[key])+'\n')
                with open(CirplotData+class_type+'_'+Sig+'_y_all.csv','w') as fout:
                    for key in sorted(sideall):
                        fout.write(str(key)+','+str(sideall[key])+'\n')  
                #save he matrix data 
                with open(CirplotData+class_type+'_'+Sig+'_mat_1.csv','w') as fout:
                    key1 = sorted(books1)[1]
                    fout.write('title')
                    for key2 in sorted(books1[key1]):
                        fout.write(','+key2)
                    fout.write('\n')
                    for key1 in sorted(books1):
                        fout.write(str(key1))
                        for key2 in sorted(books1[key1]):
                            fout.write(','+str(books1[key1][key2]))
                        fout.write('\n')
                with open(CirplotData+class_type+'_'+Sig+'_mat_2.csv','w') as fout:
                    key1 = sorted(books2)[1]
                    fout.write('title')
                    for key2 in sorted(books2[key1]):
                        fout.write(','+key2)
                    fout.write('\n')
                    for key1 in sorted(books2):
                        fout.write(str(key1))
                        for key2 in sorted(books2[key1]):
                            fout.write(','+str(books2[key1][key2]))
                        fout.write('\n')
                with open(CirplotData+class_type+'_'+Sig+'_mat_3.csv','w') as fout:
                    key1 = sorted(books3)[1]
                    fout.write('title')
                    for key2 in sorted(books3[key1]):
                        fout.write(','+key2)
                    fout.write('\n')
                    for key1 in sorted(books3):
                        fout.write(str(key1))
                        for key2 in sorted(books3[key1]):
                            fout.write(','+str(books3[key1][key2]))
                        fout.write('\n')
                with open(CirplotData+class_type+'_'+Sig+'_mat_all.csv','w') as fout:
                    key1 = sorted(booksall)[1]
                    fout.write('title')
                    for key2 in sorted(booksall[key1]):
                        fout.write(','+key2)
                    fout.write('\n')
                    for key1 in sorted(booksall):
                        fout.write(str(key1))
                        for key2 in sorted(booksall[key1]):
                            fout.write(','+str(booksall[key1][key2]))
                        fout.write('\n')

Folder=InDirec+'/change_Codon_context/'
allSubSignatures = [f for f in listdir(Folder) if os.path.isdir(os.path.join(Folder, f))]
for subsig in allSubSignatures:
    SignatureFolder=Folder+subsig+'/'
    CirplotData=SignatureFolder+'Data/'
    CirplotPlot=SignatureFolder+'Plot/'
    if not os.path.exists(CirplotData):
        os.makedirs(CirplotData)
    if not os.path.exists(CirplotPlot):
        os.makedirs(CirplotPlot)        
    #take the actual number, loop through the signatures 
    lengthOfsig=int(re.search(r'\d+', subsig).group())
    for signum in range(1,(lengthOfsig+1)):
        #define the dictionaries
        sideall=dict()
        booksall={}
        #populate the dictornary keys 
        for key1 in a:
            for key2 in b:
                sideall[key1+'_'+key2]=0
        #populate the dictornary keys 
        for types in c:
            booksall[types]={}
        #populate the dictornary keys 
        for key1 in booksall:
            for key2 in side1st:
                booksall[key1][key2]=0
        #read in the data 
        for class_type in ['all','Case','Control']:
            with open(SignatureFolder+class_type+'_'+'h_matrix.csv','r') as fin:
                title=next(fin)
                for line in fin:
                    line=line.strip('\n')
                    items=line.split(',')
                    item0=items[0].strip('\"')
                            #sig number 
                    Sig=signumber[signum]
                    number=float(items[signum])
                    #key 1 is A>C_AA
                    key1_in=item0[0]+'>'+item0[2]+'_'+item0.split('_')[-2]+item0.split('_')[-1]
                    #key2 is the code ACT
                    key2_in=item0.split('_')[2]
                    #populate data 
                    sideall[key1_in]+=number
                    booksall[key2_in][key1_in]+=number
                    #now save the matrix 
                with open(CirplotData+class_type+'_'+Sig+'_y_all.csv','w') as fout:
                    for key in sorted(sideall):
                        fout.write(str(key)+','+str(sideall[key])+'\n')  
                #save he matrix data 
                with open(CirplotData+class_type+'_'+Sig+'_mat_all.csv','w') as fout:
                    key1 = sorted(booksall)[1]
                    fout.write('title')
                    for key2 in sorted(booksall[key1]):
                        fout.write(','+key2)
                    fout.write('\n')
                    for key1 in sorted(booksall):
                        fout.write(str(key1))
                        for key2 in sorted(booksall[key1]):
                            fout.write(','+str(booksall[key1][key2]))
                        fout.write('\n')


Folder=InDirec+'/change_context/'
allSubSignatures = [f for f in listdir(Folder) if os.path.isdir(os.path.join(Folder, f))]
for subsig in allSubSignatures:
    SignatureFolder=Folder+subsig+'/'
    CirplotData=SignatureFolder+'Data/'
    CirplotPlot=SignatureFolder+'Plot/'
    if not os.path.exists(CirplotData):
        os.makedirs(CirplotData)
    if not os.path.exists(CirplotPlot):
        os.makedirs(CirplotPlot)        
    #take the actual number, loop through the signatures 
    lengthOfsig=int(re.search(r'\d+', subsig).group())
    for signum in range(1,(lengthOfsig+1)):
        #define the dictionaries
        sideall=dict()
        #populate the dictornary keys 
        for key1 in a:
            for key2 in b:
                sideall[key1+'_'+key2]=0
        #read in the data 
        for class_type in ['all','Case','Control']:
            with open(SignatureFolder+class_type+'_'+'h_matrix.csv','r') as fin:
                title=next(fin)
                for line in fin:
                    line=line.strip('\n')
                    items=line.split(',')
                    item0=items[0].strip('\"')
                            #sig number 
                    Sig=signumber[signum]
                    number=float(items[signum])
                    #key 1 is A>C_AA
                    key1_in=item0[0]+'>'+item0[2]+'_'+item0.split('_')[-2]+item0.split('_')[-1]
                    #key2 is the code ACT
                    #populate data 
                    sideall[key1_in]+=number
                    #now save the matrix 
                with open(CirplotData+class_type+'_'+Sig+'_y_all.csv','w') as fout:
                    for key in sorted(sideall):
                        fout.write(str(key)+','+str(sideall[key])+'\n')  
