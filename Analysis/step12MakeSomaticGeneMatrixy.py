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



patientMapping=dict()
allGenesWithLocationMap=dict()
allGenesMap=dict()
allVariants=dict()
genes=[]

def getSift(Varianttype,Splicing,line):
    Sift = re.search('SIFT_score=(.+?);SIFT', line)
    #Sift = re.search('Polyphen2_HDIV_score=(.+?);Polyphen2_HDIV_rankscore', line)
    #Sift = re.search('CADD_raw=(.+?);CADD_raw_rankscore=', line)
    if Sift is not None:
        Sift=Sift.group(1)
    else:
        Sift=1
    if Sift=='.':
        if Splicing == 'splicing':
            Sift=0.3
        else:
            if Varianttype=='frameshift_insertion' or Varianttype=='frameshift_deletion' or  Varianttype=='stopgain' or  Varianttype=='stoploss' or Varianttype =='frameshift_block_substitution':
                Sift=0.5
            elif Varianttype=='nonframeshift_insertion' or Varianttype=='nonframeshift_deletion'  or Varianttype=='nonframeshift_block_substitution' or Varianttype=='nonsynonymous_SNV':
                Sift=0.2
            else:
                Sift=0.1
    else:
        Sift=1-float(Sift)
    return Sift

def getPP2(Varianttype,Splicing, line):
    #Sift = re.search('SIFT_score=(.+?);SIFT', line)
    Sift = re.search('Polyphen2_HDIV_score=(.+?);Polyphen2_HDIV_rankscore', line)
    #Sift = re.search('CADD_raw=(.+?);CADD_raw_rankscore=', line)
    if Sift is not None:
        Sift=Sift.group(1)
    else:
        Sift=1
    if Sift=='.':
        if Splicing == 'splicing':
            Sift=0.3
        else:
            if Varianttype=='frameshift_insertion' or Varianttype=='frameshift_deletion' or  Varianttype=='stopgain' or  Varianttype=='stoploss' or Varianttype =='frameshift_block_substitution':
                Sift=0.5
            elif Varianttype=='nonframeshift_insertion' or Varianttype=='nonframeshift_deletion'  or Varianttype=='nonframeshift_block_substitution' or Varianttype=='nonsynonymous_SNV':
                Sift=0.2
            else:
                Sift=0.1
    else:
        Sift=float(Sift)
    return Sift


def getCADD(Varianttype,Splicing, line):
    #Sift = re.search('SIFT_score=(.+?);SIFT', line)
    #Sift = re.search('Polyphen2_HDIV_score=(.+?);Polyphen2_HDIV_rankscore', line)
    Sift = re.search('CADD_raw=(.+?);CADD_raw_rankscore=', line)
    if Sift is not None:
        Sift=Sift.group(1)
    else:
        Sift=0
    if Sift=='.':
        if Splicing == 'splicing':
            Sift=3
        else:
            if Varianttype=='frameshift_insertion' or Varianttype=='frameshift_deletion' or  Varianttype=='stopgain' or  Varianttype=='stoploss' or Varianttype =='frameshift_block_substitution':
                Sift=7
            elif Varianttype=='nonframeshift_insertion' or Varianttype=='nonframeshift_deletion'  or Varianttype=='nonframeshift_block_substitution' or Varianttype=='nonsynonymous_SNV':
                Sift=2
            else:
                Sift=1
    else:
        Sift=float(Sift)
    return Sift


Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/Somatic_Matrix/'

dirs=Input_folder
allfiles = [f for f in listdir(dirs) if isfile(join(dirs, f))]

for file in allfiles:
    if file.endswith('hg19_multianno.vcf'):
        subject=file.split('.')[0]
        patientMapping[subject]=1
        #readFile(filename,subject)

geneDict=dict()
#obtain gene names first (to build dictionary)
for sample in patientMapping:
    filedirec1=dirs+sample+'.hg19_multianno.vcf'
    for filedirec in [filedirec1]:
        with open (filedirec,'r') as fin:
            for line in fin:
                items=line.split('\t')
                if(line[0]!='#'):
                    #if 'MODERATE' in line or 'HIGH' in line:
                    if True:
                        ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                        ANN_Info_Items=ANN_Info.split(',')
                        for ANN_Info_Item in ANN_Info_Items:
                            Annotations=ANN_Info_Item.split('|')
                            Gene=Annotations[3].split('&')[0]
                            if Gene not in geneDict:
                                geneDict[Gene]=0 

books_number={}
books_sift={}
books_pp2={}
books_cadd={}
books_ale={}

for geneName in geneDict:
    books_number[geneName]={}
    books_sift[geneName]={}
    books_pp2[geneName]={}
    books_cadd[geneName]={}
    books_ale[geneName]={}
    for person in patientMapping:
        books_number[geneName][person]={}
        books_sift[geneName][person]={}
        books_pp2[geneName][person]={}
        books_cadd[geneName][person]={}
        books_ale[geneName][person]={}


for sample in patientMapping:
    #sample='A1RD'
    #print(sample)
    filedirec1=dirs+sample+'.hg19_multianno.vcf'
    for filedirec in [filedirec1]:
        with open (filedirec,'r') as fin:
            for line in fin:
                if(line[0]!='#'):
                    items=line.split('\t')
                    #if(line[0]!='#') and ('MODERATE' in line or 'HIGH' in line):
                    #if True:
                    GeneName=''
                    ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                    ANN_Info_Item=ANN_Info.split(',')[0]
                    GeneName=ANN_Info_Item.split('|')[3].split('&')[0]
                    if GeneName in geneDict:
                        location=items[0]+'_'+items[1]+'_'+items[3]+'_'+items[4]
                        genome1000 = re.search('1000g2015aug_all=(.+?);ExAC', line)
                        if genome1000 is not None:
                            genome1000=genome1000.group(1)
                        else:
                            genome1000=1
                        if genome1000=='.':
                            genome1000=0
                        #if float(genome1000)<0.05:
                        if True:
                            Splicing = re.search('Func.refGene=(.+?);Gene.refG', line)
                            if Splicing is not None:
                                Splicing=Splicing.group(1)
                            Varianttype = re.search('ExonicFunc.refGene=(.+?);AA', line)
                            if Varianttype is not None:
                                Varianttype=Varianttype.group(1)
                            Sift=getSift(Varianttype,Splicing, line)
                            CADD=getCADD(Varianttype,Splicing ,line)
                            PP2=getPP2(Varianttype,Splicing ,line)
                            ALE=items[-1].split(':')[2]
                            if location in books_number[GeneName][sample]:
                                books_number[GeneName][sample][location]+=1
                            else:
                                books_number[GeneName][sample][location]=1
                            
                            if location in books_sift[GeneName][sample]:
                                books_sift[GeneName][sample][location]+=float(Sift)
                            else:
                                books_sift[GeneName][sample][location]=float(Sift)
                            
                            if location in books_pp2[GeneName][sample]:
                                books_pp2[GeneName][sample][location]+=float(PP2)
                            else:
                                books_pp2[GeneName][sample][location]=float(PP2)
                            
                            if location in books_cadd[GeneName][sample]:
                                books_cadd[GeneName][sample][location]+=float(CADD)
                            else:
                                books_cadd[GeneName][sample][location]=float(CADD)
                            if location in books_ale[GeneName][sample]:
                                books_ale[GeneName][sample][location]+=float(ALE)
                            else:
                                books_ale[GeneName][sample][location]=float(ALE)
for GeneName in books_number:
    df = pd.DataFrame(books_number[GeneName]).T.fillna(0)
    if len(df.columns)>1:
        RootFolder=Output_folder+GeneName+'/'
        if not os.path.exists(RootFolder):
            os.makedirs(RootFolder)
        df.to_csv(Output_folder+GeneName+'/number.csv', line_terminator='\n')
        
for GeneName in books_sift:
    df = pd.DataFrame(books_sift[GeneName]).T.fillna(0)
    if len(df.columns)>1:
        RootFolder=Output_folder+GeneName+'/'
        if not os.path.exists(RootFolder):
            os.makedirs(RootFolder)
        df.to_csv(Output_folder+GeneName+'/sift.csv', line_terminator='\n')
        
for GeneName in books_pp2:
    df = pd.DataFrame(books_pp2[GeneName]).T.fillna(0)
    if len(df.columns)>1:
        RootFolder=Output_folder+GeneName+'/'
        if not os.path.exists(RootFolder):
            os.makedirs(RootFolder)
        df.to_csv(Output_folder+GeneName+'/pp2.csv', line_terminator='\n')
        
for GeneName in books_cadd:
    df = pd.DataFrame(books_cadd[GeneName]).T.fillna(0)
    if len(df.columns)>1:
        RootFolder=Output_folder+GeneName+'/'
        if not os.path.exists(RootFolder):
            os.makedirs(RootFolder)
        df.to_csv(Output_folder+GeneName+'/cadd.csv', line_terminator='\n')
        
for GeneName in books_ale:
    df = pd.DataFrame(books_ale[GeneName]).T.fillna(0)
    if len(df.columns)>1:
        RootFolder=Output_folder+GeneName+'/'
        if not os.path.exists(RootFolder):
            os.makedirs(RootFolder)
        df.to_csv(Output_folder+GeneName+'/ale.csv', line_terminator='\n')

        