# BBCAR_codes 
Codes to handle the sequencing data for publication of "Somatic Genetic Aberrations in Benign Breast Disease and the Risk of Subsequent Breast Cancer"

Please follow the steps to process the WES data. Methods as descried in the publications. Start from the directory of ‘../my_projects/Code/step1_Fastq_to_VCF.py’ This is a python script to generate .sh files for read alignment and mutation calling, including haplotype caller, mutect2, Varscan2, and VarDict. Each sample will have a .sh file generated for server submission. Follow Step 2-7 to:
	1.	python scripts to select the ‘pass’ mutations
	2.	check job status
	3.	functional annotation
	4.	select exons
	5.	MICE imputation for missing functional annotations
	6.	Convert VCF to MAF, 
	7.	Prepare readdepth for CNV
	8.	Run Varscan2 CNV
 
In addition, step 0 is to generate PON file 

Codes in "Analysis" are for the summaries and statistical analysis. 
