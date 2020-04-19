
"""
Created on Thu Oct 20 15:41:08 2016

@author: zza847
"""

import csv
import re
import sys
import math
from pandas import *
import os
import pandas as pd
import scipy.stats as stats
import numpy as np
import sys


InputDir_VCF='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
InputDir_Clinical='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/clinical_final.txt'
OutputDir='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step11data/'
class1='Case'
class2='Control'
class1_size=0
class2_size=0


sampleDict=dict()

with open(InputDir_Clinical, 'r') as f:
  title = next(f)
  for line in f:
    items=line.split('\t')
    items[1]=items[1].replace('\n','')
    if items[1] ==class1 or items[1]== class2:
        if items[1]==class1:
            class1_size=class1_size+1
        if items[1]==class2:
            class2_size=class2_size+1
        sampleDict[items[0]]=items[1]
  
outdic_D=dict()
outdic_E=dict()
outdic_F=dict()
infodict=dict()


#for sample in sorted(sampleDict, key=sampleDict.get, reverse=True)[:10]:
for sample in sampleDict:
    print(sample)
    disease_type=sampleDict[sample]
    vcf=InputDir_VCF+'/'+sample+'.hg19_multianno.vcf'        
    with open(vcf, "r") as ins:
        title=next(ins)
        for line in ins:
            Varianttype = re.search('ExonicFunc.refGene=(.+?);AA', line).group(1)
            Varianttype2=re.search('AAChange.refGene=(.+?);cytoBand', line).group(1)
            Cosmic = re.search('cosmic80=(.+?);ALLELE_END', line).group(1)
            InterPro = re.search('Interpro_domain=(.+?);1000g', line).group(1)
            GERP = re.search('GERP(.+?);GERP', line).group(1)
            ConservationScore = re.search('phastCons20way_mammalian=(.+?);phastCons20way_mammalian_rankscore', line).group(1)
            GeneName = re.search('Gene.refGene=(.+?);GeneDetail', line).group(1)
            ExAC = re.search('ExAC_ALL=(.+?);ExAC', line).group(1)
            Sift = re.search('SIFT_score=(.+?);SIFT', line).group(1)
            CADD = re.search('CADD_phred=(.+?);DANN_score', line).group(1)
            FATHMM_score = re.search('FATHMM_score=(.+?);FATHMM_converted_rankscore', line).group(1)
            PP2 = re.search('Polyphen2_HDIV_score=(.+?);Polyphen2_HDIV_rankscore', line).group(1)
            genome1000 = re.search('1000g2015aug_all=(.+?);ExAC', line).group(1)
            esp = re.search('esp6500siv2_all=(.+?);SIFT', line).group(1)
            Splicing = re.search('Func.refGene=(.+?);Gene.refG', line).group(1)
            if Cosmic != '.':
                Cosmic = 'YES'
            else:
                Cosmic = 'NO'
            items=line.split('\t')
            keyf = GeneName+'^'+ str(items[0]) +'^'+str(items[1])+'^'+ str(items[3]) +'^'+items[4]
            if disease_type==class1:   
                if keyf in outdic_D:
                    outdic_D[keyf]=outdic_D[keyf]+1
                else:
                    outdic_D[keyf]=1
            if disease_type==class2: 
                if keyf in outdic_E:
                    outdic_E[keyf]=outdic_E[keyf]+1
                else:
                    outdic_E[keyf]=1
            valueaddition=(str(Varianttype)+'^'+str(Varianttype2)+'^'+str(items[2])+'^'+items[6]+'^'+str(Splicing)+'^'
                           +str(Cosmic)+'^'+str(InterPro)+'^'+ExAC+'^'+genome1000+'^' +str(esp)
                           +'^'+Sift+'^'+CADD+'^'+str(FATHMM_score)+'^'+str(PP2)+'^'+GERP+'^'+ConservationScore   )
            if keyf in infodict:
                continue
            else:
                infodict[keyf]=valueaddition

finalresult=[]

for key in outdic_D:
    if key in outdic_E:
        result= key+'^'+str(outdic_D[key])+'^'+str(outdic_E[key])
        outdic_E.pop(key, None)
    else:
        result= key+'^'+str(outdic_D[key])+'^'+str(0)
    finalresult.append(result)

for key in outdic_E:
    result= key+'^'+str(0)+'^'+str(outdic_E[key])
    finalresult.append(result)

Pval = 0
with open(OutputDir+'variant_list.txt', "w") as f:
    f.write('Gene Name'+"\t"+'Chr'+"\t"+'Position'+"\t"+'REF'+"\t"+'ALT'+"\t"+'Benign Percentage'+"\t"+'Control Percentage'+"\t"
            +'Variant Type'+"\t"+'Variant Type2'+"\t"+'RS_number'+"\t"+'filter'+"\t"+'splicing'+"\t"+'cosmic'+"\t"+'Protein Domain'+"\t"+'EXAC'
            +"\t"+'1000 genome'+"\t"+'ESP'+"\t"+'Sift'+"\t"+'CADD'+"\t"+'FATHMM_Score'+"\t"+'PP2'+"\t"+'GERP'+"\t"+'Conservation'+"\t"+'P_value'+'\n')
    for line in finalresult:
        Finalitems=line.split('^')
        B = float(Finalitems[5])/class1_size
        C = float(Finalitems[6])/class2_size
        Case= pd.DataFrame(["1"]*int(Finalitems[5]) + ["0"]*(class1_size-int(Finalitems[5])))
        Control=pd.DataFrame(["1"]*int(Finalitems[6]) + ["0"]*(class2_size-int(Finalitems[6])))
        case_table = pd.crosstab(index=Case[0], columns="count")
        control_table = pd.crosstab(index=Control[0], columns="count")
        observed = case_table
        sumtable=control_table.add(case_table, fill_value=0)
        national_ratios = sumtable/(len(Control)+len(Case))  # Get population ratios
        expected = national_ratios * len(Case)   # Get expected counts
        subtable=observed.subtract(expected,fill_value=0)
        chi_squared_stat = ((subtable*subtable)/expected).sum()
        observed2 = control_table
        expected2 = national_ratios * len(Control)   # Get expected counts
        subtable=observed2.subtract(expected2,fill_value=0)
        chi_squared_stat2 = ((subtable*subtable)/expected2).sum()
        #crit = stats.chi2.ppf(q = 0.95, df = 2)   # Df = number of variable categories - 1
        chi_squared_stat3=chi_squared_stat2+chi_squared_stat
        Pval = 1 - stats.chi2.cdf(x=chi_squared_stat3,  df=1)
        Pval=Pval[0]
        infokey=Finalitems[0]+'^'+ Finalitems[1] +'^'+Finalitems[2]+'^'+ Finalitems[3] +'^'+Finalitems[4]
        values='.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'+'^'+'.'
        if infokey in infodict:
            values=infodict[infokey]
        values=values.split('^')
        f.write(Finalitems[0]+"\t"+Finalitems[1]+"\t"+Finalitems[2]+"\t"+Finalitems[3]+"\t"+Finalitems[4]+"\t"+str(B)+"\t"+str(C)
        +"\t"+values[0]+"\t"+values[1]+"\t"+values[2]+"\t"+values[3]+"\t"+values[4]+"\t"+values[5]+"\t"+values[6]
        +"\t"+values[7]+"\t"+values[8]+"\t"+values[9]+"\t"+values[10]+"\t"+values[11]+"\t"+values[12]+"\t"+values[13]
        +"\t"+values[14]+"\t"+values[15]
            +"\t"+str(Pval)+'\n')



