#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
params.outdir = ""
params.cpus = 1

process splitVCF {
  /*
  Split VCF file into multiple smaller VCFs based on a given number of bins

  Returns
  -------
  Returns 2 files:
      1) A VCF file
      2) A tabix index for that VCF
  */
  cpus params.cpus
  label 'bcftools'

  input:
  tuple path(vcf), path(vcf_index)
  val(bin_size)

  output:
  tuple val("${vcf}"), path("${prefix}*.vcf.gz"), path("${prefix}*.vcf.gz.tbi"), emit: files

  afterScript 'rm x*'

  script:
  """
  bcftools query -f'%CHROM\t%POS\n' ${vcf} | split -l ${bin_size}
  for file in x*; do
    bcftools view --no-version -T \${file} -Oz ${vcf} > ${prefix}.\${file}.vcf.gz
    bcftools index -t ${prefix}.\${file}.vcf.gz
  done
  """
}
