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


MicroArray_File=sys.argv[1]
Mutect_File=sys.argv[2]
OutputFile=sys.argv[3]


Mutect_Dict=dict()
with open(Mutect_File,'r') as fin:
    title=next(fin)
    for line in fin:
        items=line.split('\t')
        chro=items[0]
        pos=items[1]
        ref=items[3]
        alt=items[4]
        tumor_alt=items[-1].split(':')[2]
        tumor_depth=items[-1].split(':')[1]
        key=chro+'_'+pos
        #print(key)
        value=chro+'\t'+pos+'\t'+ref+'\t'+alt+'\t'+tumor_alt+'\t'+tumor_depth
        Mutect_Dict[key]=value


with open(MicroArray_File,'r') as fin, open (OutputFile,'w') as fout:
    title=next(fin)
    fout.write('Found_WES'+'\t'+'SNP.Name'+'\t'+'Sample.ID'+'\t'+'Chr'+'\t'+'Position'+'\t'+'GC.Score'+'\t'+'R'+'\t'+'Theta'+'\t'+'B.Allele.Freq'+'\t'+'Log.R.Ratio' +'\t'+' CNV.Value'+'\t'+'CNV.Confidence'+'\t'+'Allele1...Plus'+'\t'+'Allele2...Plus'+'\t'+'SNP'+'\t'+'Plus.Minus.Strand'+'\t'+'Allele1...Top'+'\t'+'Allele2...Top'+'\t'+'chro'+'\t'+'pos'+'\t'+'ref'+'\t'+'alt'+'\t'+'tumor_alt'+'\t'+'tumor_depth'+'\n')
    for line in fin:
        if line.strip() !='':
            items=line.strip().split('\t')
            #print(line)
            chro_micro='chr'+str(items[2].replace('\"',''))
            pos_micro=str(items[3].replace('\"',''))
            key_micro=chro_micro+'_'+pos_micro
            #print(key_micro)
            if key_micro in Mutect_Dict:
                fout.write('YES'+'\t'+line.replace('\n','').replace('\"','')+'\t'+Mutect_Dict[key_micro]+'\n')
            else:
                ref=commands.getoutput('/software/samtools/1.6/bin/samtools faidx /projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta '+chro_micro+':'+str(pos_micro)+'-'+str(pos_micro))[-1:]
                value=chro_micro+'\t'+pos_micro+'\t'+ref+'\t'+'NA'+'\t'+'0'+'\t'+'0'
                fout.write('NO'+'\t'+line.replace('\n','').replace('\"','')+'\t'+value+'\n')



