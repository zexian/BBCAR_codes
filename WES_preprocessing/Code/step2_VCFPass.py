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

inputfile='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF/'
outputfile='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/PASS_VCF/'

onlyfiles = [f for f in listdir(inputfile) ]
for file in onlyfiles:
    if file.endswith('.vcf'):
        with open(inputfile+file,'r') as fin, open (outputfile+file,'w') as fout:
            for line in fin:
                if line[0:6]=='#CHROM':
                    fout.write(line)
                if(line[0]!='#'):
                    items=line.split('\t')
                    #if items[6]=='PASS' or items[6] == 'germline_risk':
                    if items[6]=='PASS':
                        if items[0] in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]:
                            fout.write(line)

