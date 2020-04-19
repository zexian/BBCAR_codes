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
import commands
import os

def convert(Snp_ori):
    SNP=''
    if Snp_ori=='A':
        SNP='T'
    elif Snp_ori=='C':
        SNP='G'
    elif Snp_ori=='T':
        SNP='A'
    elif Snp_ori=='G':
        SNP='C'
    return(SNP)

MicroArray='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step24data/'
OutFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step25data/'

listdirs = [d for d in os.listdir(MicroArray)]

with open (OutFolder+'Array_WES_Cross.txt','w') as fout:
    fout.write('sampleID'+'\t'+'SNPName'+'\t'+'chro'+'\t'+'pos'+'\t'+'Theta'+'\t'+'BAF'+'\t'+'GCScore'+'\t'+'WES_Ref'+'\t'+'WES_Ref_frequency'+'\t'+'WES_Alt'+'\t'+'WES_Alt_frequency'+'\t'+'WES_Depth'+'\t'+'Array_Ref'+'\t'+'Array_Ref_frequency'+'\t'+'Array_Ref_raw_frequency'+'\t'+'Array_Alt'+'\t'+'Array_Alt_frequency'+'\t'+'Array_Alt_raw_frequency'+'\n')
    for file in listdirs:
        sample=file.replace('.txt','')
        with open(MicroArray+file,'r') as fin:
            title=next(fin)
            print(title)
            for line in fin:
                if line.strip() !='':
                    items=line.strip().split('\t')
                    #print(line)
                    if items[0]=='YES':
                        print(line)
                        sampleID=items[2]
                        SNPName=items[1]
                        chro='chr'+str(items[3].replace('\"',''))
                        pos=str(items[4].replace('\"',''))
                        Theta=items[7]
                        BAF=items[8]
                        GCScore=items[5]
                        WES_Ref=items[20]
                        WES_Ref_frequency=str(1-float(items[22]))
                        WES_Alt=items[21]
                        WES_Alt_frequency=items[22]
                        WES_Depth=str(int(items[23].split(',')[0])+int(items[23].split(',')[1]))

                        SNP_front=items[14].replace('[','').replace(']','').split('/')[0]
                        SNP_back=items[14].replace('[','').replace(']','').split('/')[1]
                        
                        direction=items[15]
                        if direction=='-':
                            SNP_front=convert(SNP_front)
                            SNP_back=convert(SNP_back)
                        
                        snp_dict=dict()
                        snp_dict[SNP_back]=str(float(items[8]))
                        snp_dict[SNP_front]=str(1-float(items[8]))
                        
                        snp_dict2=dict()
                        snp_dict2[SNP_back]=str(float(items[7]))
                        snp_dict2[SNP_front]=str(1-float(items[7]))

                        Array_Ref=''
                        Array_Ref_frequency=''
                        Array_Ref_frequency_raw=''
                        Array_Alt=''
                        Array_Alt_frequency=''
                        Array_Alt_frequency_raw=''

                        if WES_Ref in snp_dict:
                            Array_Ref=WES_Ref
                            Array_Ref_frequency=snp_dict[WES_Ref]
                            Array_Ref_frequency_raw=snp_dict2[WES_Ref]
                            snp_dict.pop(WES_Ref, None)
                            snp_dict2.pop(WES_Ref, None)
                            Array_Alt=list(snp_dict.keys())[0]
                            Array_Alt_frequency=snp_dict[Array_Alt]
                            Array_Alt_frequency_raw=snp_dict2[Array_Alt]
                        else:
                            Array_Ref=list(snp_dict.keys())[0]
                            Array_Ref_frequency=snp_dict[Array_Ref]
                            Array_Ref_frequency_raw=snp_dict2[Array_Ref]
                            print(snp_dict)
                            Array_Alt=list(snp_dict.keys())[1]
                            Array_Alt_frequency=snp_dict[Array_Alt]
                            Array_Alt_frequency_raw=snp_dict2[Array_Alt]

                        fout.write(sampleID+'\t'+SNPName+'\t'+chro+'\t'+pos+'\t'+Theta+'\t'+BAF+'\t'+GCScore+'\t'+WES_Ref+'\t'+WES_Ref_frequency+'\t'+WES_Alt+'\t'+WES_Alt_frequency+'\t'+WES_Depth+'\t'+Array_Ref+'\t'+Array_Ref_frequency+'\t'+Array_Ref_frequency_raw+'\t'+Array_Alt+'\t'+Array_Alt_frequency+'\t'+Array_Alt_frequency_raw+'\n')











