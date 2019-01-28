# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:54:02 2016

@author: zexian
"""

import os
from os import listdir
import commands
#Might need rename files first: for f in *_1.fq; do mv $f $(echo ${f} |sed  's/_1.fq/_R1_.fq/'); done
direc=InputDir = '/projects/p30007/Zexian/Alignment/BBCAR_NEW'


picardtool='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/picard-tools-1.131/picard.jar'
GATKtool='java -jar /projects/p30007/Zexian/tools/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
hg19Reference='/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta'
gold1000Indel ='/projects/p30007/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
dbsnp='/projects/p30007/Zexian/reference/hg19/dbsnp_138.hg19.vcf'
annovar_call='perl /projects/p30007/Zexian/tools/annovar/table_annovar.pl'
snpeff_call='java -Xmx32g -jar /projects/p30007/Zexian/tools/snpEff/snpEff.jar'
takeExon='python /projects/p30007/Zexian/tools/DNAtools/Take_exon_splicing_from_annotated_VCF.py'
HaloPlex_for_Haplotype='python /projects/p30007/Zexian/tools/DNAtools/HaloPlex_for_Haplotype.py'
varScan_call='java -Xmx32g -d64 -jar /projects/p30007/Zexian/tools/Varscan/VarScan.v2.3.9.jar'
varDict_call='/projects/p30007/Zexian/tools/VarDictJava/build/install/VarDict/bin/VarDict'
interval='/projects/p30007/Zexian/tools/DNAtools/S07604514_Padded.bed'
normal_pon='/projects/p30007/Zexian/Alignment/Germline_37/administrative/MuTect2_PON.vcf'


Raw_reads='/projects/p30007/Zexian/Alignment/BBCAR/RAW_data'
RootFolder = direc+'/WES_Analysis'
ScriptFolder = RootFolder+'/Scripts'
ScriptFolder=InputDir+'/administrative/Step2codes/'
BamFolder = RootFolder+'/BAM'
BamFolder = '/projects/b1042/lyglab/Zexian/BBCAR/BAM'
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

Haplo_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/Haplo_annovar_AnnoVCF'
Mutect_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/Mutect_annovar_AnnoVCF'
VarScan_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/VarScan_annovar_AnnoVCF'
VarDict_AnnoVCF='/projects/b1042/ClareLab/Zexian/BBCAR'+'/VarDict_annovar_AnnoVCF'

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

listdirs = [d for d in os.listdir(Raw_reads) if os.path.isdir(os.path.join(Raw_reads, d))]
FileList = []
Seqinfors=dict()
pattern= '_R1.fastq.gz'

for sample_folder in listdirs:
    SampleName=sample_folder
    print(SampleName)
    script = open(ScriptFolder+'/'+SampleName+'_full.sh', "w")
    #Print MSUB Header
    script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomicsburst
#MSUB -l walltime=160:00:00
#MSUB -m a
#MSUB -j oe
#MOAB -W umask=0113
#MSUB -l nodes=1:ppn=24
module load picard/1.131
module load bwa/0.7.12
module load samtools/1.2
module load python 
module load java/jdk1.8.0_25
module load  R/3.4.3

''')    
    bam_string=''
    bamSampleFolder=BamFolder+'/'+sample_folder
    if not os.path.exists(bamSampleFolder):
        os.makedirs(bamSampleFolder)
    for file in os.listdir(Raw_reads+'/'+sample_folder+'/'):
        if file.endswith(pattern):
            SampleName=file.split('_')[0]
            readgroup=file.split('_')[1]+file.split('_')[2]
            #Make gene list to comparison CSV
            #print(SampleName)
            FileList.append(str(SampleName))
            #Makes Shell script with file name
            #bwa mem
            R1Path = Raw_reads+'/'+sample_folder+'/'+file
            R2Path = Raw_reads+'/'+sample_folder+'/'+file.replace('R1','R2')
            subbam=bamSampleFolder+'/'+readgroup+'.bam'
            sortbam=bamSampleFolder+'/'+readgroup+'_sorted.bam'
            bam_string=bam_string+'I='+sortbam+' '
            print >> script, 'bwa mem -M -R \'@RG\\tID:'+readgroup+'\\tSM:'+SampleName+'\\tLB:library1\\tPL:ILLUMINA\\tPU:'+readgroup+'\' -t 24 '+hg19Reference+' ' +R1Path+ ' '+ R2Path+ ' | samtools view -bS - > '+subbam+'\n'
            print >> script, 'java -jar /projects/p30007/Zexian/tools/picard-tools-1.131/picard.jar SortSam INPUT='+subbam+ ' OUTPUT='+sortbam+' SORT_ORDER=coordinate\n'
            print >> script, 'rm '+subbam +'\n'
    #merge bam files 
    merged_bam=bamSampleFolder+'/'+SampleName+'_sorted_reads.bam'
    print >> script, picardtool+' MergeSamFiles ASSUME_SORTED=false CREATE_INDEX=true '+bam_string+'MERGE_SEQUENCE_DICTIONARIES=false OUTPUT='+merged_bam+' SORT_ORDER=coordinate USE_THREADING=true VALIDATION_STRINGENCY=STRICT \n'
    print >> script, 'rm '+bamSampleFolder+'/*_sorted.bam\n'
    #mark duplicates and inex
    markDuplidateBam=bamSampleFolder+'/'+SampleName+'_sorted_reads.mdup.bam'
    print >> script, picardtool+' MarkDuplicates CREATE_INDEX=true I='+merged_bam+' O='+markDuplidateBam+' M='+bamSampleFolder+'/'+SampleName+'_sorted_reads.mdup_metrics.txt\n'
    print >> script, 'rm '+merged_bam +'\n'
    #INdel-relign
    realignedBam=bamSampleFolder+'/'+SampleName+'_realigned.bam'
    print >> script, GATKtool+' -nt 24 -T RealignerTargetCreator -R '+hg19Reference+' -I '+markDuplidateBam+' -known '+gold1000Indel+' -o '+bamSampleFolder+'/'+SampleName+'_realigner.intervals\n'
    print >> script, GATKtool+' -T IndelRealigner -R '+hg19Reference+' -I '+markDuplidateBam+' -known '+gold1000Indel+' -targetIntervals '+bamSampleFolder+'/'+SampleName+'_realigner.intervals -o '+realignedBam+'\n'
    print >> script, 'rm '+markDuplidateBam +'\n'
    #base recalibrate
    recalBam=BamFolder+'/'+SampleName+'_recal_reads.bam'
    print >> script, GATKtool+' -nct 24 -T BaseRecalibrator -R '+hg19Reference+' -I '+realignedBam+' -knownSites '+dbsnp+' -knownSites '+gold1000Indel+' -L '+interval+' -o '+bamSampleFolder+'/'+SampleName+'_recal_data.table\n'
    print >> script, GATKtool+' -nct 24 -T PrintReads -R '+hg19Reference+' -I '+realignedBam+' -BQSR '+bamSampleFolder+'/'+SampleName+'_recal_data.table -o '+recalBam+'\n'
    print >> script, 'rm '+realignedBam +'\n'
    #haplogypecaller 
    #Haplo_VCF_File=Haplo_VCF+'/'+SampleName+'.vcf'
    #print >> script, GATKtool+' -nct 8 -T HaplotypeCaller -R '+hg19Reference+' -I '+recalBam+' --dbsnp '+dbsnp+' --genotyping_mode DISCOVERY -stand_emit_conf 20 -stand_call_conf 30 -o '+Haplo_VCF_File+'\n'
    #print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+' -V '+Haplo_VCF_File+' -selectType SNP -o '+Haplo_VCF+'/'+SampleName+'_raw_snps.vcf\n'
    #print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+Haplo_VCF+'/'+SampleName+'_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o '+Haplo_VCF+'/'+SampleName+'_filtered_snps.vcf\n'
    #print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+'  -V '+Haplo_VCF_File+' -selectType INDEL -o '+Haplo_VCF+'/'+SampleName+'_raw_indels.vcf\n'
    #print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+Haplo_VCF+'/'+SampleName+'_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter"  -o '+Haplo_VCF+'/'+SampleName+'_filtered_indel.vcf\n'
    #HaloPlex for haplotype results /also pass 
    #print >> script, 'module load python/anaconda3.6'
    #print >> script, HaloPlex_for_Haplotype+' '+Haplo_VCF+'/'+SampleName+'_filtered_snps.vcf '+Filter_Haplo_VCF+'/'+SampleName +'_snp.vcf  \n'
    #print >> script, HaloPlex_for_Haplotype+' '+Haplo_VCF+'/'+SampleName+'_filtered_indel.vcf '+Filter_Haplo_VCF+'/'+SampleName +'_indel.vcf  \n'
    # #mutect
    Mutect_VCF_File=Mutect_VCF+'/'+SampleName+'.vcf'
    print >> script, GATKtool+' -nct 24 -T MuTect2 -R '+hg19Reference+' -I:tumor '+recalBam+' --cosmic /projects/p30007/Zexian/reference/hg19/CosmicCodingMuts_chr_M_sorted.vcf --dbsnp '+dbsnp+ ' -L '+interval+' --normal_panel '+normal_pon+' -o '+ Mutect_VCF_File+' --output_mode EMIT_VARIANTS_ONLY'
    script.close()



    ########germline take  exon###########
    #print >> script, 'module load python/anaconda3.6'
    #print >> script, takeExon+' '+GFinalVCF+'/'+SampleName+'_annotated_snp.vcf  '+GExonVCF+'/'+SampleName+'_exon_snp.vcf '+GExon_5per_VCF+'/'+SampleName+'_exon_snp.vcf '+GExon_1per_VCF+'/'+SampleName+'_exon_snp.vcf\n'
    #print >> script, takeExon+' '+GFinalVCF+'/'+SampleName+'_annotated_indel.vcf  '+GExonVCF+'/'+SampleName+'_exon_indel.vcf '+GExon_5per_VCF+'/'+SampleName+'_exon_indel.vcf '+GExon_1per_VCF+'/'+SampleName+'_exon_indel.vcf\n'

    # #mutect
    # Mutect_VCF_File=Mutect_VCF+'/'+SampleName+'.vcf'
    # print >> script, GATKtool+' -nct 8 -T MuTect2 -R '+hg19Reference+' -I:tumor '+recalBam+' --cosmic /projects/p30007/Zexian/reference/hg19/CosmicCodingMuts_chr_M_sorted.vcf --dbsnp '+dbsnp+' -o '+ Mutect_VCF_File+' --output_mode EMIT_VARIANTS_ONLY'
 
    #mutect2 call vcf
#    print >> script, GATKtool+' -nct 24 -T MuTect2 -R '+hg19Reference+' -I:tumor '+recalBam+' --cosmic /projects/p30007/Zexian/reference/hg19/CosmicCodingMuts_chr_M_sorted.vcf --dbsnp '+dbsnp+' -o '+ MutectVcffile+'\n'
#    print >> script, '/projects/p30007/Zexian/tools/VarDictJava/build/install/VarDict/bin/VarDict -G '+hg19Reference+' -th 16 -f 0.02 -N '+SampleName+' -b '+recalBam+' -c 1 -S 2 -E 3 -g 4 /projects/p30007/Zexian/tools/DNAtools/self_defined.bed | /projects/p30007/Zexian/tools/VarDictJava/VarDict/teststrandbias.R | /projects/p30007/Zexian/tools/VarDictJava/VarDict/var2vcf_valid.pl -N '+SampleName+' -f 0.02 > '+ VarDict_VCF+ '/' + SampleName+'.vcf\n'
#    
#    #test quality
#    print >> script, 'samtools flagstat '+ recalBam+' > '+Quality+'/'+SampleName+'_flagstat.txt' 
#    #germline varints gfilter
#    print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+' -V '+GVcffile+' -selectType SNP -o '+GVcfFolder+'/'+SampleName+'_raw_snps.vcf\n'
#    print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+GVcfFolder+'/'+SampleName+'_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o '+GVcfFolder+'/'+SampleName+'_filtered_snps.vcf\n'
#    print >> script, GATKtool+' -T SelectVariants -R '+hg19Reference+'  -V '+GVcffile+' -selectType INDEL -o '+GVcfFolder+'/'+SampleName+'_raw_indels.vcf\n'
#    print >> script, GATKtool+' -T VariantFiltration -R '+hg19Reference+'  -V '+GVcfFolder+'/'+SampleName+'_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter"  -o '+GVcfFolder+'/'+SampleName+'_filtered_indel.vcf\n'
#    #########Germline annotation###########
#    print >> script, annovar_call+' --thread 5 '+GVcfFolder+'/'+SampleName+'_filtered_snps.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+GAnnoVCF+'/'+SampleName+'_snp -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, annovar_call+' --thread 5 '+GVcfFolder+'/'+SampleName+'_filtered_indel.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+GAnnoVCF+'/'+SampleName+'_indel -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+GAnnoVCF+'/'+SampleName+'_snp.hg19_multianno.vcf > '+GFinalVCF+'/'+SampleName+'_annotated_snp.vcf\n'
#    print >> script, snpeff_call+' hg19 '+GAnnoVCF+'/'+SampleName+'_indel.hg19_multianno.vcf > '+GFinalVCF+'/'+SampleName+'_annotated_indel.vcf\n'
#    #########Somatic annotation###########
#    print >> script, annovar_call+' --thread 5 '+MutectVcffile+' /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+SAnnoVCF+'/'+SampleName+' -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+SAnnoVCF+'/'+SampleName+'.hg19_multianno.vcf > '+MutectFinalVCF+'/'+SampleName+'_annotated.vcf\n'
#    #vardict annoation
#    print >> script, annovar_call+' --thread 5 '+VarDict_VCF+ '/' + SampleName+'.vcf /projects/p30007/Zexian/tools/annovar/humandb/ -buildver hg19 -out '+VarDict_Anno_VCF+'/'+SampleName+' -remove -protocol avsnp147,refGene,cytoBand,cosmic80,dbnsfp31a_interpro,1000g2015aug_all,exac03,esp6500siv2_all,dbnsfp33a,revel,clinvar_20170130,snp138NonFlagged -operation f,g,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput\n' 
#    print >> script, snpeff_call+' hg19 '+VarDict_Anno_VCF+'/'+SampleName+'.hg19_multianno.vcf > '+Vardict_FinalVCF+'/'+SampleName+'_annotated.vcf\n'
#    print >> script, 'python /projects/p30007/Zexian/tools/DNAtools/Preprocess_for_Cosmic_vardict.py '+VarDict_VCF+ '/' + SampleName+'.vcf '+ Vardict_preprocess+'/'+SampleName+'_somatic.vcf \n'
#    print >> script, 'perl /projects/p30007/Zexian/tools/ISOWN_update/bin/database_annotation.pl '+Vardict_preprocess+'/'+SampleName+'_somatic.vcf  '+Vardict_midprocess+'/'+SampleName+'.vcf \n'
#    #########call somatic mutation###########
#    print >> script, 'python /projects/p30007/Zexian/tools/DNAtools/Preprocess_for_Cosmic_Mutect.py '+MutectVcffile +' '+ Mutect_preprocess+'/'+SampleName+'_somatic.vcf \n'
#    print >> script, 'perl /projects/p30007/Zexian/tools/ISOWN_update/bin/database_annotation.pl '+Mutect_preprocess+'/'+SampleName+'_somatic.vcf  '+Mutect_midprocess+'/'+SampleName+'.vcf \n'
#    print >> script, 'python /projects/p30007/Zexian/tools/DNAtools/Preprocess_for_Cosmic_Vardict.py '+VardictVcfFile +' '+ Vardict_preprocess+'/'+SampleName+'_somatic.vcf \n'
#    print >> script, 'perl /projects/p30007/Zexian/tools/ISOWN_update/bin/database_annotation.pl '+Vardict_preprocess+'/'+SampleName+'_somatic.vcf  '+Vardict_midprocess+'/'+SampleName+'.vcf \n'


    #########summary seq information##########
#    fastq_read_number=Germ_snp=Germ_snp=Mutect_CALL=final_somatic=Germ_indel=final_somatic_exon=coverage=str(0)
#    fastq_read_number=int( commands.getoutput('zcat '+R1Path+' | wc -l'))/4.0
#    Germ_snp= commands.getoutput('grep -c PASS '+GVcfFolder+'/'+SampleName+'_filtered_snps.vcf')
#    Germ_indel=commands.getoutput('grep -c PASS '+GVcfFolder+'/'+SampleName+'_filtered_indel.vcf')
#    Mutect_CALL=commands.getoutput('grep -c PASS '+MutectVcffile)
#    final_somatic=0
#    if  os.path.exists(MutectFilteredVCF+'/'+SampleName+'_annotated.vcf'):
#        final_somatic=commands.getoutput('grep -c PASS '+MutectFilteredVCF+'/'+SampleName+'_annotated.vcf')
#    if  os.path.exists(Somatic_Exon+'/'+SampleName+'_annotated.vcf'):
#        final_somatic_exon=commands.getoutput('grep -c PASS '+Somatic_Exon+'/'+SampleName+'_annotated.vcf')
#    coverage=commands.getoutput('python /projects/p30007/Zexian/tools/DNAtools/get_read_depth.py '+InputDir+' '+ SampleName)
#    coverage=commands.getoutput('samtools depth -b /projects/p30007/Zexian/tools/DNAtools/self_defined.bed '+ BamFolder+'/'+SampleName+'_recal_reads.bam  |  awk \'{sum+=$3; sumsq+=$3*$3} END { print \"Average = \",sum/NR; print \"Stdev = \",sqrt(sumsq/NR - (sum/NR)**2)}\'')
#    print(coverage)
#    Seqinfors[str(SampleName)]=[fastq_read_number,Germ_snp,Germ_indel,Mutect_CALL,final_somatic,final_somatic_exon,coverage]


#############################m
#make the code for somatic calling
# somatic_caller_script_upload = open(Analysis+'/somatic_caller_script_upload.sh', "w")
# somatic_caller_script_upload.write('''#!/bin/bash
# #MSUB -A b1042
# #MSUB -q genomics
# #MSUB -l walltime=47:00:00
# #MSUB -m a
# #MSUB -j oe
# #MOAB -W umask=0113
# #MSUB -l nodes=1:ppn=4
# #MSUB -N V1-6
# module load java/jdk1.8.0_25
# module load python/anaconda3
# module load R/3.3.1
# module load gcc/5.1.0
# ''')    
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/take_exon_for_Cosmic.py '+MutectFinalVCF+' '+Mutect_Exon_before_ISOWN+' '+Mutect_Exon_Rare_before_ISOWN+' '+Mutect_Exon_before_ISOWN_ALL+ '\n'
# print >> somatic_caller_script_upload, 'perl /projects/p30007/Zexian/tools/ISOWN_update/bin/run_isown.pl '+ Mutect_midprocess+' '+Somatic+'/Mutect_somatic.txt " -trainingSet /projects/p30007/Zexian/tools/ISWON/ISOWN/training_data/BRCA_100_TrainSet.arff -sanityCheck false -classifier nbc " \n'
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/take_exon_for_Cosmic.py '+Vardict_FinalVCF+' '+Vardict_Exon_before_ISOWN+' '+Vardict_Exon_Rare_before_ISOWN+ '\n'
# print >> somatic_caller_script_upload, 'perl /projects/p30007/Zexian/tools/ISOWN_update/bin/run_isown.pl '+ Vardict_midprocess+' '+Somatic+'/Vardict_somatic.txt " -trainingSet /projects/p30007/Zexian/tools/ISWON/ISOWN/training_data/BRCA_100_TrainSet.arff -sanityCheck false -classifier nbc " \n'

# print >> somatic_caller_script_upload, 'module load python/anaconda3.6'+ '\n'
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/after_process_for_Cosmic.py '+Somatic+'/Mutect_somatic.txt '+MutectFinalVCF +' ' +MutectFilteredVCF+ '\n'
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/take_exon_for_Cosmic.py '+MutectFilteredVCF+' '+Mutect_Exon+' '+Mutect_Exon_Rare+ '\n'
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/after_process_for_Cosmic.py '+Somatic+'/Vardict_somatic.txt '+Vardict_FinalVCF +' ' +Vardict_FilteredVCF+ '\n'
# print >> somatic_caller_script_upload, 'python /projects/p30007/Zexian/tools/DNAtools/take_exon_for_Cosmic.py '+Vardict_FilteredVCF+' '+Vardict_Exon+' '+Vardict_Exon_Rare+ '\n'
# somatic_caller_script_upload.close()

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
# with open(Analysis+'/comparison_temps.csv', "w") as ComparisonCSV:
#     ComparisonCSV.write('ID'+','+'case_control'+','+'batch'+','+'menopausal'+','+'benign_age'+','+'cancer_age'+','+'interval_days'+','+'race'+','+'ER'+','+'class'+','+'Fastq_read'+','+'Germline_SNP'+','+'Germline_INDEL'+','+'MUTECT_Call'+','+'Somatic'+','+'Somatic_Exon'+','+'Coverage'+'\n')
#     for file in FileList:
#         file =str(file)
#         ComparisonCSV.write(str(file)+','+str(infor[file][1])+','+str(infor[file][0])+','+str(infor[file][2])+','+str(infor[file][3])+','+str(infor[file][4])+','+str(infor[file][5])+','+str(infor[file][6])+','+str(infor[file][7])+','+str(infor[file][8])+','+str(Seqinfors[file][0])+','+str(Seqinfors[file][1])+','+str(Seqinfors[file][2])+','+str(Seqinfors[file][3])+','+str(Seqinfors[file][4])+','+str(Seqinfors[file][5])+','+str(Seqinfors[file][6])+'\n')
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
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Germ.py '+ InputDir+ ' Case ' + 'Control\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Soma_before_ISOWN.py '+ InputDir+ ' Case ' + 'Control\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_Varinats_list_BBCAR_Soma.py '+ InputDir+ ' Case ' + 'Control\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/find_top_20_cancer.py '+ InputDir+'\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/Make_All_gene_list_BBCAR_Soma.py '+ InputDir+ ' Case ' + 'Control\n'
# print >> analysis_script, 'python /projects/p30007/Zexian/tools/DNAtools/vcf_to_a_file.py '+ InputDir+'\n'
# print >> analysis_script, 'Rscript /projects/p30007/Zexian/tools/DNAtools/R_to_plot_variants.R '+ InputDir+'\n'


# #print >> analysis_script, 'Rscript /projects/p30007/Zexian/tools/DNAtools/Gene_Signature_Study.R '+ InputDir+'\n'
# analysis_script.close()


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
# #MSUB -N V1-6
# module load java/jdk1.8.0_25
# module load python/anaconda3
# module load R/3.3.1
# module load gcc/5.1.0
# ''')    
# print >> signature_script_upload, 'Rscript /projects/p30007/Zexian/tools/DNAtools/Gene_Signature_Study.R '+ InputDir+'\n'
# signature_script_upload.close()
