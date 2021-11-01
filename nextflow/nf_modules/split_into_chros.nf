#!/usr/bin/env nextflow

/* 
 * Script to split a multi-chromosome VCF into single-chromosome VCFs
 *
 * @author
 * Likhitha Surapaneni <likhitha3996@gmail.com>
 *
 */
nextflow.enable.dsl=2
// params defaults

prefix = "out"
cpus = 1
params.outdir = ""

process splitVCF {
	/*
	Function to split a multi-chromosome VCF into single chromosome VCF

	Returns
	-------
	Returns 2 files per chromosome:
		1) A VCF format file for each splitted chromosome
		2) A tabix index for that VCF
	*/
	
	input:
        val(chr)
        path(vcf)


	output:
	tuple path("${prefix}.${chr}.vcf.gz"), path("${prefix}.${chr}.vcf.gz.tbi")

	"""
        bcftools index -t ${vcf}
	bcftools view -r ${chr} ${vcf} -o ${prefix}.${chr}.vcf.gz -O z
	bcftools index -t ${prefix}.${chr}.vcf.gz
	"""
}
