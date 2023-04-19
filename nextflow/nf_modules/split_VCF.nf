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
  tuple path(vcf), path(vcf_index), path(split_file), path(vep_config), val(output_dir)

  output:
  tuple val("${vcf}"), path("${prefix}*.vcf.gz"), path("${prefix}*.vcf.gz.tbi"), path(vep_config), val(output_dir)

  afterScript 'rm x*'

  script:
  """
  bcftools view --no-version -T ${split_file} -Oz ${vcf} > ${prefix}.${split_file}.vcf.gz
  bcftools index -t ${prefix}.${split_file}.vcf.gz
  """
}
