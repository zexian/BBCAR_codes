# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 22:24:19 2016

@author: zza847
"""
direc=InputDir='/projects/p30007/Zexian/Alignment/Germline_37'
direc_normal='/projects/p30007/Zexian/Alignment/BBCAR/RAW_data'


from os import listdir
from os.path import isfile, join
import csv
import os
import os.path
import commands

picardtool='java -jar /projects/p30007/Zexian/tools/picard-tools-1.131/picard.jar'
GATKtool='java -jar /projects/p30007/Zexian/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
hg19Reference='/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta'
gold1000Indel ='/projects/p30007/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp='/projects/p30007/Zexian/reference/hg19/dbsnp_138.hg19.vcf'
annovar_call='perl /projects/p30007/Zexian/tools/annovar/table_annovar.pl'
snpeff_call='java -Xmx32g -jar /projects/p30007/Zexian/tools/snpEff/snpEff.jar'
HaloPlex_for_Haplotype='python /projects/p30007/Zexian/tools/DNAtools/HaloPlex_for_Haplotype.py'
takeExon='python /projects/p30007/Zexian/tools/DNAtools/Take_exon_splicing_from_annotated_VCF.py'
varScan_call='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/Varscan/VarScan.v2.3.9.jar'
varDict_call='/projects/p30007/Zexian/tools/VarDictJava/build/install/VarDict/bin/VarDict'
interval='/projects/p30007/Zexian/tools/DNAtools/S07604514_Padded.bed'

Raw_reads=direc+'/RAW_data'
RootFolder = direc+'/WES_Analysis'
ScriptFolder = InputDir+'/administrative/Step3codes/'
OutputFolder_normal = InputDir+'/administrative/Step2Data/'
OutFolder = InputDir+'/administrative/'


BamFolder = RootFolder+'/BAM'
BamFolder = '/projects/b1042/lyglab/Zexian/Germline_37/BAM'
SampleDepthDirec=RootFolder+'/Analysis/Sample_Depth/'

HaplotypeFolder=RootFolder+'/Haplotype'
MutecFolder=RootFolder+'/Mutect'
VarScanFolder=RootFolder+'/VarScan'
VarDcitFolder=RootFolder+'/VarDict'

Haplo_VCF=HaplotypeFolder+'/VCF'
Mutect_VCF=MutecFolder+'/VCF'
VarScan_VCF=VarScanFolder+'/VCF'
VarDict_VCF=VarDcitFolder+'/VCF'

Filter_Haplo_VCF=HaplotypeFolder+'/Filter_VCF'


Haplo_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/Haplo_annovar_AnnoVCF'
Mutect_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/Mutect_annovar_AnnoVCF'
VarScan_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/VarScan_annovar_AnnoVCF'
VarDict_AnnoVCF='/projects/b1042/ClareLab/Zexian/Germline_37'+'/VarDict_annovar_AnnoVCF'

Haplo_FinalVCF=HaplotypeFolder+'/annotated_VCF'
Mutect_FinalVCF=MutecFolder+'/annotated_VCF'
VarScan_FinalVCF=VarScanFolder+'/annotated_VCF'
VarDict_FinalVCF=VarDcitFolder+'/annotated_VCF'

Haplo_Exon=HaplotypeFolder+'/exon_VCF'
Haplo_Exon_5per=HaplotypeFolder+'/exon_5per_VCF'
Haplo_Exon_1per=HaplotypeFolder+'/exon_1per_VCF'

Mutec_Exon=MutecFolder+'/exon_VCF'
Mutec_Exon_5per=MutecFolder+'/exon_5per_VCF'
Mutec_Exon_1per=MutecFolder+'/exon_1per_VCF'

VarScan_Exon=VarScanFolder+'/exon_VCF'
VarScan_Exon_5per=VarScanFolder+'/exon_5per_VCF'
VarScan_Exon_1per=VarScanFolder+'/exon_1per_VCF'

VarDict_Exon=VarDcitFolder+'/exon_VCF'
VarDict_Exon_5per=VarDcitFolder+'/exon_5per_VCF'
VarDict_Exon_1per=VarDcitFolder+'/exon_1per_VCF'

Analysis = RootFolder+'/Analysis'
Variant_plot=Analysis+'/Variant_plot'

Mutect_preprocess=MutecFolder+'/ISOWN/Preprocess'
Mutect_midprocess=MutecFolder+'/ISOWN/Midprocess'
VarScan_preprocess=VarScanFolder+'/ISOWN/Preprocess'
VarScan_midprocess=VarScanFolder+'/ISOWN/Midprocess'
VarDict_preprocess=VarDcitFolder+'/ISOWN/Preprocess'
VarDict_midprocess=VarDcitFolder+'/ISOWN/Midprocess'

Varscan_CNV=direc+'CNV/VarScan'
Varscan_CNV_mileup='/projects/b1042/lyglab/Zexian/Germline_37/mileup'

Validated_Samples=dict()
with open (InputDir+'/administrative/ValidaSamples.csv','r') as fin:
    title=next(fin)
    for line in fin:
        print(line.strip())
        Validated_Samples[line.strip()]=1

filename_dict=dict()
bamfile_dict=dict()
cancerfilename_dict=dict()
cancerbamfile_dict=dict()
listdirs = [d for d in os.listdir(OutputFolder_normal) if d.endswith('.vcf') ]
FileList = []
Seqinfors=dict()

string_normal=''
for file in listdirs:
    string_normal=string_normal+'-V '+OutputFolder_normal+file+' '

        #print(sample_folder)
#if not os.path.exists(Vscan_Sample):
#    os.makedirs(Vscan_Sample)
script = open(ScriptFolder+'/'+'step3generatorePON'+'.sh', "w")
#Print MSUB Header
script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=47:59:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=12
module load picard/1.131
module load bwa/0.7.12
module load python
module load samtools/0.1.18 
module load java/jdk1.8.0_25
module load  R

''')    

print >> script, GATKtool+' -T CombineVariants -R '+hg19Reference+' '+string_normal+'-minN 2  --setKey "null" --filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED -L '+interval+' -o '+ OutFolder+'MuTect2_PON.vcf'
script.close()


 #       bam_string_normal=''
 #       bam_string_blood=''
    #    bamSampleFolder=BamFolder+'/'+sample_folder
    #     if not os.path.exists(bamSampleFolder):
    #         os.makedirs(bamSampleFolder)
    #     for file in os.listdir(direc_normal+'/'+sample_folder+'/'):
    #         if file.endswith(pattern):
    #             SampleName=file.split('_')[0]
    #             readgroup=file.split('_')[1]+file.split('_')[2]
    #             #Make gene list to comparison CSV
    #             FileList.append(str(SampleName))
    #             R1Path = direc_normal+'/'+sample_folder+'/'+file
    #             R2Path = direc_normal+'/'+sample_folder+'/'+file.replace('R1','R2')
    #             subbam=bamSampleFolder+'/'+readgroup+'.bam'
    #             sortbam=bamSampleFolder+'/'+readgroup+'_sorted.bam'
    #             bam_string_normal=bam_string_normal+'I='+sortbam+' '
    #             print >> script, 'bwa mem -M -R \'@RG\\tID:'+readgroup+'\\tSM:'+SampleName+'normal\\tLB:library1\\tPL:ILLUMINA\\tPU:'+readgroup+'\' -t 12 '+hg19Reference+' ' +R1Path+ ' '+ R2Path+ ' | samtools view -bS - > '+subbam+'\n'
    #             print >> script, picardtool+' SortSam INPUT='+subbam+ ' OUTPUT='+sortbam+' SORT_ORDER=coordinate\n'
    #             print >> script, 'rm '+subbam +'\n'
    #     for file in os.listdir(Raw_reads+'/'+sample_folder+'/'):
    #         if file.endswith(pattern):
    #             SampleName=file.split('_')[0]
    #             readgroup=file.split('_')[1]+file.split('_')[2]
    #             #Make gene list to comparison CSV
    #             #print(SampleName)
    #             FileList.append(str(SampleName))
    #             R1Path = Raw_reads+'/'+sample_folder+'/'+file
    #             R2Path = Raw_reads+'/'+sample_folder+'/'+file.replace('R1','R2')
    #             subbam=bamSampleFolder+'/'+readgroup+'.bam'
    #             sortbam=bamSampleFolder+'/'+readgroup+'_sorted.bam'
    #             bam_string_blood=bam_string_blood+'I='+sortbam+' '
    #             print >> script, 'bwa mem -M -R \'@RG\\tID:'+readgroup+'\\tSM:'+SampleName+'blood\\tLB:library1\\tPL:ILLUMINA\\tPU:'+readgroup+'\' -t 12 '+hg19Reference+' ' +R1Path+ ' '+ R2Path+ ' | samtools view -bS - > '+subbam+'\n'
    #             print >> script, picardtool+' SortSam INPUT='+subbam+ ' OUTPUT='+sortbam+' SORT_ORDER=coordinate\n'
    #             print >> script, 'rm '+subbam +'\n'
    # #    #merge bam files 
    #     merged_bam_normal=bamSampleFolder+'/'+SampleName+'_normal_sorted_reads.bam'
    #     merged_bam_blood=bamSampleFolder+'/'+SampleName+'_blood_sorted_reads.bam'
    #     print >> script, picardtool+' MergeSamFiles ASSUME_SORTED=false CREATE_INDEX=true '+bam_string_normal+'MERGE_SEQUENCE_DICTIONARIES=false OUTPUT='+merged_bam_normal+' SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT \n'
    #     print >> script, picardtool+' MergeSamFiles ASSUME_SORTED=false CREATE_INDEX=true '+bam_string_blood+'MERGE_SEQUENCE_DICTIONARIES=false OUTPUT='+merged_bam_blood+' SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT \n'
    #     print >> script, 'rm '+bamSampleFolder+'/*_sorted.bam\n'
    # #    #mark duplicates and inex
    #     markDuplidateBam_normal=bamSampleFolder+'/'+SampleName+'_normal_sorted_reads.mdup.bam'
    #     markDuplidateBam_blood=bamSampleFolder+'/'+SampleName+'_blood_sorted_reads.mdup.bam'
    #     print >> script, picardtool+' MarkDuplicates CREATE_INDEX=true I='+merged_bam_normal+' O='+markDuplidateBam_normal+' M='+bamSampleFolder+'/'+SampleName+'_normal_sorted_reads.mdup_metrics.txt\n'
    #     print >> script, picardtool+' MarkDuplicates CREATE_INDEX=true I='+merged_bam_blood+' O='+markDuplidateBam_blood+' M='+bamSampleFolder+'/'+SampleName+'_blood_sorted_reads.mdup_metrics.txt\n'
    #     print >> script, 'rm '+merged_bam_normal +'\n'
    #     print >> script, 'rm '+merged_bam_blood +'\n'
    # #    #INdel-relign
    #     realignedBam_normal=bamSampleFolder+'/'+SampleName+'_normal_realigned.bam'
    #     realignedBam_blood=bamSampleFolder+'/'+SampleName+'_blood_realigned.bam'
    #     print >> script, GATKtool+' -nt 12 -T RealignerTargetCreator -R '+hg19Reference+' -I '+markDuplidateBam_normal+' -known '+gold1000Indel+' -o '+bamSampleFolder+'/'+SampleName+'_normal_realigner.intervals\n'
    #     print >> script, GATKtool+' -nt 12 -T RealignerTargetCreator -R '+hg19Reference+' -I '+markDuplidateBam_blood+' -known '+gold1000Indel+' -o '+bamSampleFolder+'/'+SampleName+'_blood_realigner.intervals\n'
    #     print >> script, GATKtool+' -T IndelRealigner -R '+hg19Reference+' -I '+markDuplidateBam_normal+' -known '+gold1000Indel+' -targetIntervals '+bamSampleFolder+'/'+SampleName+'_normal_realigner.intervals -o '+realignedBam_normal+'\n'
    #     print >> script, GATKtool+' -T IndelRealigner -R '+hg19Reference+' -I '+markDuplidateBam_blood+' -known '+gold1000Indel+' -targetIntervals '+bamSampleFolder+'/'+SampleName+'_blood_realigner.intervals -o '+realignedBam_blood+'\n'
    #     print >> script, 'rm '+markDuplidateBam_normal +'\n'
    #     print >> script, 'rm '+markDuplidateBam_blood +'\n'
    #     #base recalibrate
    #    recalBam_normal=bamSampleFolder+'/'+SampleName+'_normal_recal_reads.bam'
    #    recalBam_blood=bamSampleFolder+'/'+SampleName+'_blood_recal_reads.bam'
    #     print >> script, GATKtool+' -nct 12 -T BaseRecalibrator -R '+hg19Reference+' -I '+realignedBam_normal+' -knownSites '+dbsnp+' -knownSites '+gold1000Indel+' -L '+interval+' -o '+bamSampleFolder+'/'+SampleName+'_normal_recal_data.table\n'
    #     print >> script, GATKtool+' -nct 12 -T BaseRecalibrator -R '+hg19Reference+' -I '+realignedBam_blood+' -knownSites '+dbsnp+' -knownSites '+gold1000Indel+' -o '+bamSampleFolder+'/'+SampleName+'_blood_recal_data.table\n'
    #     print >> script, GATKtool+' -nct 12 -T PrintReads -R '+hg19Reference+' -I '+realignedBam_normal+' -BQSR '+bamSampleFolder+'/'+SampleName+'_normal_recal_data.table -o '+recalBam_normal+'\n'
    #     print >> script, GATKtool+' -nct 12 -T PrintReads -R '+hg19Reference+' -I '+realignedBam_blood+' -BQSR '+bamSampleFolder+'/'+SampleName+'_blood_recal_data.table -o '+recalBam_blood+'\n'
    #     print >> script, 'rm '+realignedBam_normal +'\n'
    #     print >> script, 'rm '+realignedBam_blood +'\n'
# #   #call haplogypecaller
#     Haplo_VCF_File=Haplo_VCF+'/'+SampleName+'.vcf'
#     print >> script, GATKtool+' -nct 12 -T HaplotypeCaller -R '+hg19Reference+' -I '+recalBam_blood+' --dbsnp '+dbsnp+' --genotyping_mode DISCOVERY -stand_emit_conf 20 -stand_call_conf 30 -o '+Haplo_VCF_File+'\n'
#     print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+' -V '+Haplo_VCF_File+' -selectType SNP -o '+Haplo_VCF+'/'+SampleName+'_raw_snps.vcf\n'
#     print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+Haplo_VCF+'/'+SampleName+'_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o '+Haplo_VCF+'/'+SampleName+'_filtered_snps.vcf\n'
#     print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+'  -V '+Haplo_VCF_File+' -selectType INDEL -o '+Haplo_VCF+'/'+SampleName+'_raw_indels.vcf\n'
#     print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+Haplo_VCF+'/'+SampleName+'_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter"  -o '+Haplo_VCF+'/'+SampleName+'_filtered_indel.vcf\n'
#     #HaloPlex for haplotype results /also pass 
#    print >> script, 'module load python/anaconda3.6'
#    print >> script, HaloPlex_for_Haplotype+' '+Haplo_VCF+'/'+SampleName+'_filtered_snps.vcf '+Filter_Haplo_VCF+'/'+SampleName +'_snp.vcf  \n'
#    print >> script, HaloPlex_for_Haplotype+' '+Haplo_VCF+'/'+SampleName+'_filtered_indel.vcf '+Filter_Haplo_VCF+'/'+SampleName +'_indel.vcf  \n'





    # Mutect_VCF_File=Mutect_VCF+'/'+SampleName+'.vcf'
    # print >> script, GATKtool+' -nct 12 -T MuTect2 -R '+hg19Reference+' -I:tumor '+recalBam_normal+' -I:normal ' +recalBam_blood+' --cosmic /projects/p30007/Zexian/reference/hg19/CosmicCodingMuts_chr_M_sorted.vcf --dbsnp '+dbsnp+' -o '+ Mutect_VCF_File+' --output_mode EMIT_VARIANTS_ONLY'
    # print >>script, 'samtools mpileup -B -q 1 -f '+hg19Reference+' ' +recalBam_blood+' '+recalBam_normal+' | '+varScan_call+' somatic -mpileup '+ VarScan_VCF+'/'+SampleName+ '_normal --min-coverage 20 --min-coverage-normal 20 --min-coverage-tumor 20 --min-var-freq 0.02 --min-freq-for-hom 0.75 --normal-purity 1.0 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf 1  \n'
    # print >>script, varScan_call+' processSomatic '+VarScan_VCF+'/'+SampleName+'_normal.snp.vcf --p-value 0.05'
    # print >>script, varScan_call+' processSomatic '+VarScan_VCF+'/'+SampleName+'_normal.indel.vcf --p-value 0.05'
    # print >> script, varDict_call +' -G '+hg19Reference+' -th 12 -f 0.02 -N '+SampleName+'normal -b \"'+recalBam_normal+'|'+recalBam_blood+'\" -c 1 -S 2 -E 3 /projects/p30007/Zexian/tools/DNAtools/self_defined.bed | /projects/p30007/Zexian/tools/VarDictJava/VarDict/testsomatic.R | /projects/p30007/Zexian/tools/VarDictJava/VarDict/var2vcf_paired.pl -N \"'+SampleName+'normal|'+SampleName+'blood\" -f 0.02 > '+ VarDict_VCF+ '/' + SampleName+'.vcf\n'
    
    #call ISOWN mutations
    #print >> script, 'python /projects/p30007/Zexian/tools/DNAtools/Preprocess_for_Varscan_to_Isown.py '+Varscan_VCF+ '/' + SampleName +'normal '+Mutect_VCF+ '/' + SampleName+'.vcf '+ Somatic_preprocess+'/'+SampleName+'_somatic.vcf \n'
    #print >> script, 'perl /projects/p30007/Zexian/tools/ISOWN_update//bin/database_annotation.pl '+Somatic_preprocess+'/'+SampleName+'_somatic.vcf  '+Somatic_midprocess+'/'+SampleName+'.vcf \n'

    #annotation and take exon
#    print >> script, 'module load python/anaconda3.6'
#    print >> script, annovar_call+' --thread 5 '+Haplo_VCF+'/'+SampleName+'_filtered_snps.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+Haplo_AnnoVCF+'/'+SampleName+'_snp -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, annovar_call+' --thread 5 '+Haplo_VCF+'/'+SampleName+'_filtered_indel.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+Haplo_AnnoVCF+'/'+SampleName+'_indel -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+Haplo_AnnoVCF+'/'+SampleName+'_snp.hg19_multianno.vcf > '+Haplo_FinalVCF+'/'+SampleName+'_annotated_snp.vcf\n'
#    print >> script, snpeff_call+' hg19 '+Haplo_AnnoVCF+'/'+SampleName+'_indel.hg19_multianno.vcf > '+Haplo_FinalVCF+'/'+SampleName+'_annotated_indel.vcf\n'
#    print >> script, takeExon+' '+Haplo_FinalVCF+'/'+SampleName+'_annotated_snp.vcf  '+Haplo_Exon+'/'+SampleName+'_exon_snp.vcf '+Haplo_Exon_5per+'/'+SampleName+'_exon_snp.vcf '+Haplo_Exon_1per+'/'+SampleName+'_exon_snp.vcf\n'
#    print >> script, takeExon+' '+Haplo_FinalVCF+'/'+SampleName+'_annotated_indel.vcf  '+Haplo_Exon+'/'+SampleName+'_exon_indel.vcf '+Haplo_Exon_5per+'/'+SampleName+'_exon_indel.vcf '+Haplo_Exon_1per+'/'+SampleName+'_exon_indel.vcf\n'
#
#    print >> script, annovar_call+' --thread 5 '+Mutect_VCF_File+' /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+Mutect_AnnoVCF+'/'+SampleName+' -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+Mutect_AnnoVCF+'/'+SampleName+'.hg19_multianno.vcf > '+Mutect_FinalVCF+'/'+SampleName+'_annotated.vcf\n'
#    print >> script, takeExon+' '+Mutect_FinalVCF+'/'+SampleName+'_annotated.vcf  '+Mutec_Exon+'/'+SampleName+'_exon.vcf '+Mutec_Exon_5per+'/'+SampleName+'_exon.vcf '+Mutec_Exon_1per+'/'+SampleName+'_exon.vcf\n'
#
#    all_annofile = [f for f in listdir(VarScan_VCF) if isfile(join(VarScan_VCF, f))]
#    for annofile in all_annofile:
#        gene_scr=SampleName+'_normal'
#        if gene_scr in annofile:
#            print >> script, annovar_call+' --thread 5 '+VarScan_VCF+'/'+annofile+' /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+SVARAnnoVCF+'/'+annofile[:-4]+' -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#            print >> script, snpeff_call+' hg19 '+SVARAnnoVCF+'/'+annofile[:-4]+'.hg19_multianno.vcf > '+VarScan_FinalVCF+'/'+annofile[:-4]+'_annotated.vcf\n'
#            print >> script, takeExon+' '+VarScan_FinalVCF+'/'+annofile[:-4]+'_annotated.vcf  '+VarScan_Exon+'/'+annofile[:-4]+'_exon.vcf '+VarScan_Exon_5per+'/'+annofile[:-4]+'_exon.vcf '+VarScan_Exon_1per+'/'+annofile[:-4]+'_exon.vcf\n'
#
#    print >> script, annovar_call+' --thread 4 '+VarDict_VCF+'/'+SampleName+'.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+VarDict_AnnoVCF+'/'+SampleName+' -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+VarDict_AnnoVCF+'/'+SampleName+'.hg19_multianno.vcf > '+VarDict_FinalVCF+'/'+SampleName+'_annotated.vcf\n'
#    print >> script, takeExon+' '+VarDict_FinalVCF+'/'+SampleName+'_annotated.vcf  '+VarDict_Exon+'/'+SampleName+'_exon.vcf '+VarDict_Exon_5per+'/'+SampleName+'_exon.vcf '+VarDict_Exon_1per+'/'+SampleName+'_exon.vcf\n'
#    
    #call copy number variation using varscan
#    with open (SampleDepthDirec+'/'+SampleName+'_blood_depth.txt') as fin:
#        blood_depth=next(fin)
#    with open (SampleDepthDirec+'/'+SampleName+'_normal_depth.txt') as fin:
#        normal_depth=next(fin)
#    data_ratio=float(float(blood_depth)/float(normal_depth))
#    Copynumber_copynumber=Varscan_CNV+'/'+SampleName+'.copynumber'
#    Copynumber_called=Varscan_CNV+'/'+SampleName+'.copynumber.called'
#    Copynumber_called_homdel=Varscan_CNV+'/'+SampleName+'.copynumber.called.homdel'
#    Copynumber_called_recenter=Varscan_CNV+'/'+SampleName+'.copynumber.called.recenter'
#    Copynumber_seg=Varscan_CNV+'/'+SampleName+'.seg.pvalue'
#    Copynumber_seg_merged=Varscan_CNV+'/'+SampleName+'.seg.pvalue.mergeSegment'
#    markerPosition_file=Varscan_CNV+'/'+SampleName+'.markerPosition'
#    print >>script, 'samtools mpileup -q 1 -f '+hg19Reference+' ' +recalBam_blood+' ' +recalBam_normal+' |awk \'{if($4 >= 6) print $0}\' | awk \'{if($7 != 0) print $0}\' | '+varScan_call+' copynumber -mpileup '+ Varscan_CNV+'/'+SampleName+' --mpileup 1 --p-value 0.01 --min-coverage 20 --min-map-qual 20 --min-base-qual 20 --data-ratio ' +str(data_ratio)+'\n'
#    print >>script, varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called +' --output-homdel-file '+Copynumber_called_homdel+'\n'
#    print >>script, 'delta=$(python /projects/p30007/Zexian/tools/varscan2/meanLogRatioByChromosome.py '+Copynumber_called+')'
#    print >>script, 'cmp=$(awk -v delta=$delta \'END{if (delta < -0.01) {print "lt"} else {if (delta > 0.01) {print "gt"} else {print "eq"}}}\' < /dev/null)'
#    print >>script, 'if [[ "$cmp" == "lt" ]]; then'
#    print >>script, '   rd=$(echo $delta | sed \'s/-//\')'
#    print >>script, '   '+varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called_recenter +' --output-homdel-file '+Copynumber_called_homdel+' --recenter-down $rd'
#    print >>script, 'elif [[ "$cmp" == "gt" ]]; then'
#    print >>script, '    '+varScan_call+' copyCaller '+Copynumber_copynumber+' --output-file '+Copynumber_called_recenter +' --output-homdel-file '+Copynumber_called_homdel+' --recenter-up $delta'
#    print >>script, 'else'
#    print >>script, '    cp '+Copynumber_copynumber+' '+Copynumber_called_recenter
#    print >>script, 'fi'
#    print >>script, 'Rscript /projects/p30007/Zexian/tools/DNAtools/CopyNumberR.R '+ Copynumber_called_recenter+' ' + Copynumber_seg+' '+SampleName+' '+markerPosition_file+'\n'
#    print >>script, 'perl /projects/p30007/Zexian/tools/DNAtools/mergeSegment.pl '+Copynumber_seg+' --ref-arm-sizes /projects/p30007/Zexian/tools/DNAtools/hg19.len.bed --output-basename '+Copynumber_seg_merged+'\n'
#    print >>script, 'python /projects/p30007/Zexian/tools/DNAtools/Clean_CBS_arm_file3.py '+ Varscan_CNV+'/'+SampleName+'.seg.pvalue.mergeSegment.events.tsv ' + Varscan_CNV+'/'+SampleName+ '.segmentedFile.clean.txt '+ SampleName+'\n'
##
#    ########summary seq information##########
#    fastq_read_number=Germ_snp=Germ_snp=Mutect_CALL=final_somatic=Germ_indel=mutec_somatic_exon=coverage=mutec_somatic=str(0)
#    #fastq_read_number=int( commands.getoutput('zcat '+R1Path+' | wc -l'))/4.0
#    #Germ_snp= commands.getoutput('grep -c PASS '+GVcfFolder+'/'+SampleName+'_filtered_snps.vcf')
#    #Germ_indel=commands.getoutput('grep -c PASS '+GVcfFolder+'/'+SampleName+'_filtered_indel.vcf')
#    #Mutect_CALL=commands.getoutput('grep -c PASS '+SVcffile)
#    #final_somatic=0
#    #if  os.path.exists(SFilteredVCF+'/'+SampleName+'_annotated.vcf'):
#    #mutec_somatic=commands.getoutput('grep -c PASS '+SMutecFinalVCF+'/'+SampleName+'_annotated.vcf')
#    #if  os.path.exists(Somatic_Exon+'/'+SampleName+'_annotated.vcf'):
#    #mutec_somatic_exon=commands.getoutput('grep -c PASS '+Mutec_Exon+'/'+SampleName+'_exon.vcf')
#    #coverage=commands.getoutput('python /projects/p30007/Zexian/tools/DNAtools/get_read_depth_paired.py '+InputDir+' '+ SampleName)
#    Seqinfors[str(SampleName)]=[fastq_read_number,Germ_snp,Germ_indel,Mutect_CALL,mutec_somatic,mutec_somatic_exon,coverage]



# #make the code for repetitive CNV for Varscan
# patient_dict=dict()
# with open ('/projects/p30007/Zexian/Alignment/Germline_37/WES_Analysis/Analysis/comparison.csv','r') as fin:
#     f_line=next(fin)
#     for line in fin:
#         items=line.split(',')
#         patient_dict[items[0]]=items[1]
# case_seg=''
# case_mark=''
# control_seg=''
# control_mark=''
# all_seg=''
# all_mark=''
# for patient in patient_dict:
#     disease=patient_dict[patient]
#     if disease=='Case' and patient != '1025':
#         case_seg=case_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
#         case_mark=case_mark+Varscan_CNV+'/'+patient+'.markerPosition '
#     if disease=='Control' and patient != '1401':
#         control_seg=control_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
#         control_mark=control_mark+Varscan_CNV+'/'+patient+'.markerPosition '
#     if (disease=='Control' or disease =='Case' ) and patient != '1401' and patient != '1025' and patient != '652':
#         all_seg=all_seg+Varscan_CNV+'/'+patient+'.segmentedFile.clean.txt '
#         all_mark=all_mark+Varscan_CNV+'/'+patient+'.markerPosition '
# CNV_visulization = open(Analysis+'/preprare_CNV.txt', "w")
# print >> CNV_visulization, 'cd /projects/p30007/Zexian/tools/GISTIC_2_0_23/'
# #print >> signature_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/Clean_vcf_for_signature_paired.py '+ InputDir+'\n'
# print >> CNV_visulization, 'cat '+case_mark+' > ' +Varscan_CNV+'/Case_normal.clean.markerPosition\n'
# print >> CNV_visulization, 'cat '+case_seg+' > ' +Varscan_CNV+'/Case_segmentedFile.clean.txt\n'
# print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV+'/Case/ -seg '+Varscan_CNV+'/Case_segmentedFile.clean.txt -mk '+Varscan_CNV+'/Case_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+ Varscan_CNV+'/Case/log.txt \n'

# print >> CNV_visulization, 'cat '+control_mark+' > ' +Varscan_CNV+'/Control_normal.clean.markerPosition\n'
# print >> CNV_visulization, 'cat '+control_seg+' > ' +Varscan_CNV+'/Control_segmentedFile.clean.txt\n'
# print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV+'/Control/ -seg '+Varscan_CNV+'/Control_segmentedFile.clean.txt -mk '+Varscan_CNV+'/Control_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+Varscan_CNV+'/Control/log.txt \n'

# print >> CNV_visulization, 'cat '+all_mark+' > ' +Varscan_CNV+'/All_normal.clean.markerPosition\n'
# print >> CNV_visulization, 'cat '+all_seg+' > ' +Varscan_CNV+'/All_segmentedFile.clean.txt\n'
# print >>CNV_visulization, './gistic2 gp_gistic2_from_seg -b '+ Varscan_CNV+'/All/ -seg '+Varscan_CNV+'/All_segmentedFile.clean.txt -mk '+Varscan_CNV+'/All_normal.clean.markerPosition -refgene /projects/p30007/Zexian/tools/GISTIC_2_0_23/refgenefiles/hg19.mat -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme > '+Varscan_CNV+'/All/log.txt\n'
# CNV_visulization.close()
            
            
           
# #############################m
# #make big comparison table 
# infor=dict()
# with open (Analysis+'/Cohort_Data3.csv','r') as fin:
#     firstrow=next(fin)
#     print(firstrow)
#     for line in fin:
#         items=line.split(',')
#         infor[str(items[1])]=[items[0],items[17],items[4], items[6], items[9], items[13],items[14], items[16], items[22] ]  #batch id ,case_control menopausal status, benign age, cancer age, interval days , race, ER, prolifiration class

# FileList.sort()
# print(len(infor))
# print(len(Seqinfors))
# with open(Analysis+'/comparison_temp.csv', "w") as ComparisonCSV:
#     ComparisonCSV.write('ID'+','+'case_control'+','+'batch'+','+'menopausal'+','+'Removed'+','+'benign_age'+','+'cancer_age'+','+'interval_days'+','+'race'+','+'ER'+','+'class'+','+'Fastq_read'+','+'Varscan_SNP'+','+'Varscan_INDEL'+','+'MUTECT_Call'+','+'Mutect_Somatic'+','+'Mutect_Somatic_Exon'+','+'Coverage'+'\n')
#     for file in FileList:
#         file =str(file)
#         ComparisonCSV.write(str(file)+','+str(infor[file][1])+','+str(infor[file][0])+','+str(infor[file][2])+',NO,'+str(infor[file][3])+','+str(infor[file][4])+','+str(infor[file][5])+','+str(infor[file][6])+','+str(infor[file][7])+','+str(infor[file][8])+','+str(Seqinfors[file][0])+','+str(Seqinfors[file][1])+','+str(Seqinfors[file][2])+','+str(Seqinfors[file][3])+','+str(Seqinfors[file][4])+','+str(Seqinfors[file][5])+','+str(Seqinfors[file][6])+'\n')
# ComparisonCSV.close()
# analysis_script = open(Analysis+'/make_big_table_germline_upload.sh', "w")
# analysis_script.write('''#!/bin/bash
# #MSUB -A b1042
# #MSUB -q genomics
# #MSUB -l walltime=24:00:00
# #MSUB -m a
# #MSUB -j oe
# #MOAB -W umask=0113
# #MSUB -l nodes=1:ppn=4
# #MSUB -N V1-6
# module load java/jdk1.8.0_25
# module load python/anaconda3
# module load R/3.3.1
# ''')    
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Alignment_summary.py '+  InputDir+'\n'
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Germ.py '+ InputDir+ ' Case ' + 'Control\n'
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Soma_before_ISOWN.py '+ InputDir+ ' Case ' + 'Control\n'
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Soma.py '+ InputDir+ ' Case ' + 'Control\n'
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/find_top_20_cancer.py '+ InputDir+'\n'
# #print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_gene_list_BBCAR_Soma.py '+ InputDir+ ' Case ' + 'Control\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/vcf_to_a_file_paired_mode.py '+ InputDir+'\n'
# print >> analysis_script, 'Rscript /projects/p30007/Zexian/tools/DNAtools/R_to_plot_variants.R '+ InputDir+'\n'
# analysis_script.close()
# #############################m#############################m
           

# #############################m
# #make the code for gene signature
# signature_script_upload = open(Analysis+'/signature_script_upload.sh', "w")
# signature_script_upload.write('''#!/bin/bash
# #MSUB -A b1042
# #MSUB -q genomics
# #MSUB -l walltime=24:00:00
# #MSUB -m a
# #MSUB -j oe
# #MOAB -W umask=0113
# #MSUB -l nodes=1:ppn=4
# module load java/jdk1.8.0_25
# module load python/anaconda3
# module load R/3.3.1
# module load gcc/5.1.0
# ''')    
# #print >> signature_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/Clean_vcf_for_signature_paired.py '+ InputDir+'\n'

# print >> signature_script_upload, 'Rscript /projects/p30007/Zexian/tools/DNAtools/Gene_Signature_Study_BBCAR_PAIR.R '+ InputDir+'\n'
# signature_script_upload.close()
            
            
            
            
            
            
                