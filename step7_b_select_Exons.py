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

InputVCF='/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Mutect/PASS_VCF_Anno/'
OutputVCF='/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Mutect/PASS_VCF_Anno_Exonic/'

onlyfiles = [f for f in listdir(InputVCF) if isfile(join(InputVCF, f))]
for file in onlyfiles:
    if file.endswith('.vcf'):
        with open(InputVCF+file,'r') as fin, open (OutputVCF+file,'w') as fout:
            for line in fin:
                if line[0:6]=='#CHROM':
                    fout.write(line)
                if(line[0]!='#'):
                    items=line.split('\t')
                    infors=items[-2]
                    infoitems=infors.split(':')
                    if ',' in infoitems[1]:
                        ADs=infoitems[1].split(',')
                        AD1=int(ADs[0])
                        AD2=int(ADs[1])
                        DP_num=AD1+AD2
                        if DP_num>=20:
                            if 'Func.refGene=exonic' in line or 'Func.refGene=exonic;splicing' in line or 'Func.refGene=UTR3' in line or 'Func.refGene=UTR5' in line or 'Func.refGene=UTR5;UTR3' in line:
                                if 'snp138NonFlagged=r' in line  and 'cosmic80=.' in line:
                                    print('Germline')
                                else:
                                    fout.write(line)
onlyfiles = [f for f in listdir(InputVCF) if isfile(join(InputVCF, f))]
for file in onlyfiles:
    if file.endswith('.txt'):
        with open(InputVCF+file,'r') as fin, open (OutputVCF+file,'w') as fout:
            title=next(fin)
            fout.write(title)
            for line in fin:
                items=line.split('\t')
                infors=items[-2].replace('\n','')
                infoitems=infors.split(':')
                if ',' in infoitems[1]:
                    ADs=infoitems[1].split(',')
                    AD1=int(ADs[0])
                    AD2=int(ADs[1])
                    DP_num=AD1+AD2
                    if DP_num>=20:
                        reffun=items[6]
                        if 'exonic' == reffun or 'exonic;splicing' == reffun  or 'UTR3' == reffun or 'UTR5' == reffun or 'UTR5;UTR3' == reffun:
                            if items[99][0]=='r' and items[100]=='.':
                                print('Germline')
                            else:
                                fout.write(line)






