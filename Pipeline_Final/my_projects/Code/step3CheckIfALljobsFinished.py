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

inputfile='/projects/p30007/Zexian/Alignment/Cancer_10/WES_Analysis/Mutect/VCF/'

onlyfiles = [f for f in listdir(inputfile) ]

for file in onlyfiles:
    if file.endswith('.vcf'):
        indi='NO'
        with open(inputfile+file,'r') as fin:
            for line in fin:
                if(line[0]!='#'):
                    if 'chrX' in line:
                        indi='YES'
        if indi=='NO':
            print(file)
