# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:54:02 2016

@author: zexian
"""

import os
from os import listdir
import commands
from os.path import isfile, join
import re
from os import listdir
import sys

#Might need rename files first: for f in *_1.fq; do mv $f $(echo ${f} |sed  's/_1.fq/_R1_.fq/'); done


CirplotCode='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/step19Cirplot.R'
CirplotCode_V2='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/step19Cirplot_remove_position.R'
CirplotCode_V3='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/step19Cirplot_remove_position_context.R'

ScriptFolder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step19codes/'


InDirec='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step17data/'
	
filetype='change_Codon_position_context'
TypeFolder=InDirec+filetype+'/'
allSubTypeFolder = [f for f in listdir(TypeFolder) if os.path.isdir(os.path.join(TypeFolder, f))]
for folderss in allSubTypeFolder:
	SignatureFolder=TypeFolder+'/'+folderss+'/'
	script = open(ScriptFolder+filetype+'_'+folderss+'.sh', "w")
	    #Print MSUB Header
	script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=9:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=1
module load picard/1.131
module load bwa/0.7.12
module load samtools/1.2
module load python 
module load java/jdk1.8.0_25
module load  R/3.4.3

''')    
	print >> script, 'Rscript '+ CirplotCode+' '+SignatureFolder +'\n'
	script.close()

filetype='change_Codon_context'
TypeFolder=InDirec+filetype+'/'
allSubTypeFolder = [f for f in listdir(TypeFolder) if os.path.isdir(os.path.join(TypeFolder, f))]
for folderss in allSubTypeFolder:
	SignatureFolder=TypeFolder+'/'+folderss+'/'
	script = open(ScriptFolder+filetype+'_'+folderss+'.sh', "w")
	    #Print MSUB Header
	script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=9:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=1
module load picard/1.131
module load bwa/0.7.12
module load samtools/1.2
module load python 
module load java/jdk1.8.0_25
module load  R/3.4.3

''')    
	print >> script, 'Rscript '+ CirplotCode_V2+' '+SignatureFolder +'\n'
	script.close()





filetype='change_context'
TypeFolder=InDirec+filetype+'/'
allSubTypeFolder = [f for f in listdir(TypeFolder) if os.path.isdir(os.path.join(TypeFolder, f))]
for folderss in allSubTypeFolder:
	SignatureFolder=TypeFolder+'/'+folderss+'/'
	script = open(ScriptFolder+filetype+'_'+folderss+'.sh', "w")
	    #Print MSUB Header
	script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=9:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=1
module load picard/1.131
module load bwa/0.7.12
module load samtools/1.2
module load python 
module load java/jdk1.8.0_25
module load  R/3.4.3

''')    
	print >> script, 'Rscript '+ CirplotCode_V3+' '+SignatureFolder +'\n'
	script.close()




