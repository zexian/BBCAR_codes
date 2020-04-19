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
import glob
import os

Input_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/VCF_Anno_Exonic_Somatic/'
Output_folder='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step15data/'

removefiles = glob.glob(Output_folder+'/*')
for f in removefiles:
    os.remove(f)

AA_File='/projects/p30007/Zexian/Alignment/TCGA_BRCA/AA_Change.txt'
AA_Change_Dict=dict()
with open(AA_File,'r') as fin:
    for line in fin:
        items=line.split('\t')
        #print repr(items[1])
        AA_Change_Dict[items[0]]=items[1].strip('\r\n')

RefDict=dict()
with open ('/projects/p30007/Zexian/Alignment/TCGA_Signature/administrative/refSeqDataBase.fa','r') as fin:
    for line in fin:
        items=line.split('\t')
        ID=items[0]
        seq=items[1].replace('\n','')
        if ID not in RefDict:
            RefDict[ID]=seq
def ChangeStrand(Stringin):
    s=list(Stringin)
    for index,val in enumerate(s):
        sub_val='N'
        if val=='A':
            sub_val='T'
        if val=='T':
            sub_val='A'
        if val=='C':
            sub_val='G'
        if val=='G':
            sub_val='C'
        s[index]=sub_val
    StringOut=''.join(s)
    return StringOut

def parseString(Varaint_Change,AA_Change,variant_type):
    Ref_Secondary=Varaint_Change[-3]
    Alt_Secondary=Varaint_Change[-1]
    Variant_Location=re.search(r'\d+', Varaint_Change).group()
    AA_Ref=AA_Change[0:3]
    AA_ALt=AA_Change[-3:]
    AA_Location=AA_Change[3:-3]
    if 'stop_gained' in variant_type:
        AA_ALt=AA_Change[-1:]
        AA_Location=AA_Change[3:-1]
    if 'stop_lost' in variant_type or 'start_lost' in variant_type:
        AA_ALt=AA_Change[-1:]
        AA_Location=re.search(r'\d+', AA_Change).group()
    return Ref_Secondary,Alt_Secondary,Variant_Location,AA_Ref,AA_ALt,AA_Location

            # #(('missense_variant', 'Asn', 'Thr', 'C', 'NM_001114106.2', '556', '1667', '1777, 'AA_Change_Dict','forward')
def ReadCode(variant_type,AA_Ref,AA_ALt,Alt_Secondary,TID,AA_Location,Variant_Location,Transcript_position,AA_Change_Dict,direction):
    AA_3time=int(AA_Location)*3
    offset=AA_3time-int(Variant_Location)
    Mutation_positio=3-offset
    seq=RefDict[TID]
    startPosition=int(Transcript_position)-6+offset
    endPosition=int(Transcript_position)+3+offset
    if variant_type=='start_lost':
        startPosition=startPosition+3
    if endPosition > len(seq):
        seq=seq+'AAA'
    letter9=seq[startPosition:endPosition]
    if variant_type=='start_lost':
        letter9=letter9+'AAA'
    return Mutation_positio,letter9

def TakeFromDNA(direction,chro,chro_position,offset):
    startPosition=int(chro_position)-1
    endPosition=int(chro_position)+1
    AA_startPosition=int(chro_position)-2+offset
    AA_endPosition=int(chro_position)+offset
    letter3=commands.getoutput('/software/samtools/1.6/bin/samtools faidx /projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta '+chro+':'+str(startPosition)+'-'+str(endPosition))[-3:]
    AA_letter3=commands.getoutput('/software/samtools/1.6/bin/samtools faidx /projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta '+chro+':'+str(AA_startPosition)+'-'+str(AA_endPosition))[-3:]
    letter3=letter3.upper()
    AA_letter3=AA_letter3.upper()
    
    if direction=='reverse':
        letter3=letter3[::-1]
        letter3=ChangeStrand(letter3)
        
        AA_startPosition=int(chro_position)-offset
        AA_endPosition=int(chro_position)+2-offset
        AA_letter3=commands.getoutput('/software/samtools/1.6/bin/samtools faidx /projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta '+chro+':'+str(AA_startPosition)+'-'+str(AA_endPosition))[-3:]
        AA_letter3=AA_letter3[::-1]
        AA_letter3=ChangeStrand(AA_letter3)
    upstream=letter3[0]
    downstream=letter3[2]
    return(AA_letter3,upstream,downstream)

iiq=0
print('start')
onlyfiles = [f for f in listdir(Input_folder) if isfile(join(Input_folder, f))]
print('in total file number: ' + str(len(onlyfiles)))
for file in onlyfiles:
    if file.endswith('.vcf'):
        print(file)
        #iiq=iiq+1
        #print(iiq)
        with open(Input_folder+file,'r') as fin, open (Output_folder+file.replace('vcf','csv'),'w') as fout:
            fout.write('variant_type'+'\t'+'significance'+'\t'+'Gene'+'\t'+'chro'+'\t'+'chro_position'+'\t'+'Ref'+'\t'+'Alt'+'\t'+'Gene_Direction'+'\t'+'Ref_after_reverse'+'\t'+'Alt_after_reverse'+'\t'+'AA_Code'+'\t'+'AA_Name'+'\t'+'Mutated_AA_Position'+'\t'+'upstream'+'\t'+'downstream'+'\t'+'ALE'+'\t'+'Depth'+'\n')
            for line in fin:
                if line[0:6]=='#CHROM':
                    titleline=line
                if(line[0]!='#'):
                    if 'missense_variant' in line or 'synonymous_variant' in line or 'stop_gained' in line or 'stop_lost' in line or 'stop_retained_variant' in line or 'start_lost' in line:
                        #print(line)
                        result=''
                        cur_line=line
                        items=line.split('\t')
                        chro=items[0]
                        chro_position=items[1]
                        Ref=items[3]
                        Alt=items[4]
                        ALE=items[-1].split(':')[2]
                        #print(items[-1].split(':')[1])
                        Depth=str(int(items[-1].split(':')[1].split(',')[0])+int(items[-1].split(':')[1].split(',')[1]))
                        if items[3] in ['A','T','G','C']:
                            if items[4] in ['A','T','G','C']:
                                infor=items[7]
                                ANN_Info = re.search('ANN=(.+?);ANNOVAR_DATE', line).group(1)
                                ANN_Info_Items=ANN_Info.split(',')
                                for ANN_Info_Item in ANN_Info_Items:
                                    if 'missense_variant' in ANN_Info_Item or 'synonymous_variant' in  ANN_Info_Item or 'stop_gained' in ANN_Info_Item or 'stop_lost' in line or 'stop_retained_variant' in line or 'start_lost' in line:
                                        start_end='no'
                                        if 'splice_region_variant' in ANN_Info_Item:
                                            start_end='yes'
                                        #print(ANN_Info_Item)
                                        Annotations=ANN_Info_Item.split('|')
                                        #print(Annotations)
                                        variant_type=Annotations[1]
                                        significance= Annotations[2]
                                        Gene=Annotations[3]
                                        TID=Annotations[6]
                                        TID_items=TID.split('.')
                                        if len(TID_items)>2:
                                            TID=TID_items[0]+'.'+TID_items[1]
                                        Varaint_Change= Annotations[9].replace('c.','')   #Variant_Change: 1667A>C   #
                                        #print(Annotations[10])
                                        AA_Change=Annotations[10].replace('p.','')      ##AA_Change : Asn556Thr     #stop gained : Trp202*
                                        pos1=Annotations[11]         #1777/2397
                                        Transcript_position=pos1.split('/')[0]
                                        seq_length=pos1.split('/')[1]
                                        pos2=Annotations[12]    #1667/1962
                                        pos3=Annotations[13]    #556/653
                                        #('A', 'C', '1667', 'Asn', 'Thr', '556')
                                        Ref_Secondary,Alt_Secondary,Variant_Location,AA_Ref,AA_ALt,AA_Location = parseString(Varaint_Change,AA_Change,variant_type)
                                        AA_3time=int(AA_Location)*3
                                        offset=AA_3time-int(Variant_Location)
                                        Mutation_positio=3-offset
                                        direction = 'forward'
                                        if Ref!=Ref_Secondary:
                                            direction='reverse'
                                            #2,'TTTAACTAC'
                                        if TID in RefDict:
                                            seq=RefDict[TID]
                                            Mutation_positio,letter9=ReadCode(variant_type,AA_Ref,AA_ALt,Alt_Secondary,TID,AA_Location,Variant_Location,Transcript_position,AA_Change_Dict,direction)
                                            #print(letter9)                   #(('missense_variant', 'Asn', 'Thr', 'C', 'NM_001114106.2', '556', '1667', '1777, 'AA_Change_Dict','forward')
                                            letter_m3=letter9[3:-3]
                                            upstream=seq[int(Transcript_position)-2]
                                            if int(Transcript_position)>=len(seq):
                                                downstream=letter9[3+Mutation_positio]
                                            else:
                                                downstream=seq[int(Transcript_position)]    
                                            Code_Reference,check_upstream,check_downstream=TakeFromDNA(direction,chro,chro_position,offset)
                                            if letter_m3 in AA_Change_Dict:                                 
                                                AA_m3=AA_Change_Dict[letter_m3].strip('\n')
                                                s=list(letter_m3)
                                                s[Mutation_positio-1] = Alt_Secondary
                                                AA_ALt_Code=''.join(s)
                                                if 'stop_gained' in variant_type or 'stop_lost' in variant_type or 'start_lost' in variant_type:
                                                    AA_ALt=AA_Change_Dict[AA_ALt_Code]
                                                variant_type=variant_type.replace('&splice_region_variant','').replace('splice_region_variant&','')
                                                if (check_upstream==upstream) and (check_downstream==downstream) and (AA_Ref==AA_m3) and (AA_ALt==AA_Change_Dict[AA_ALt_Code]) and (Code_Reference==letter_m3): #quality check
                                                    result=variant_type+'\t'+significance+'\t'+Gene+'\t'+chro+'\t'+chro_position+'\t'+Ref+'\t'+Alt+'\t'+direction+'\t'+Ref_Secondary+'\t'+Alt_Secondary+'\t'+letter_m3+'\t'+AA_m3+'\t'+str(Mutation_positio)+'\t'+upstream+'\t'+downstream+'\t'+ALE+'\t'+Depth+'\n'
                                                    fout.write(result) 
                                                    break
                                                else:
                                                    print('in')
                                                    print(Code_Reference)
                                                    print(letter_m3)
                                                    Code_Reference,upstream,downstream=TakeFromDNA(direction,chro,chro_position,offset)
                                                    result=variant_type+'\t'+significance+'\t'+Gene+'\t'+chro+'\t'+chro_position+'\t'+Ref+'\t'+Alt+'\t'+direction+'\t'+Ref_Secondary+'\t'+Alt_Secondary+'\t'+Code_Reference+'\t'+AA_Ref+'\t'+str(Mutation_positio)+'\t'+upstream+'\t'+downstream+'\t'+ALE+'\t'+Depth+'\n'
                                                    fout.write(result)
                                                    break
                                        else:
                                            print('TID is not in dict')
                                            Code_Reference,upstream,downstream=TakeFromDNA(direction,chro,chro_position,offset)
                                            result=variant_type+'\t'+significance+'\t'+Gene+'\t'+chro+'\t'+chro_position+'\t'+Ref+'\t'+Alt+'\t'+direction+'\t'+Ref_Secondary+'\t'+Alt_Secondary+'\t'+Code_Reference+'\t'+AA_Ref+'\t'+str(Mutation_positio)+'\t'+upstream+'\t'+downstream+'\t'+ALE+'\t'+Depth+'\n'
                                            fout.write(result)
                                            break
                                       

                        
                        

