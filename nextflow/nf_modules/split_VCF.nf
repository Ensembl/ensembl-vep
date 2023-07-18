#!/usr/bin/env nextflow

/* 
 * Script to split a VCF file into multiple smaller VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"

process splitVCF {
  /*
  Split VCF file into multiple smaller VCFs based on a given number of bins

  Returns
  -------
  Tuple of original VCF, split VCF files, split VCF index files, vep config file, a output dir, and the index type of VCF file
  */
  
  cpus params.cpus
  label 'bcftools'

  input:
  tuple path(vcf), path(vcf_index), path(split_file), path(vep_config), val(index_type)

  output:
  tuple val("${vcf}"), path("${prefix}*.vcf.gz"), path("${prefix}*.vcf.gz.{tbi,csi}"), path(vep_config), val(index_type)

  afterScript 'rm x*'

  script:
  index_flag = index_type == "tbi" ? "-t" : "-c"
  
  """
  bcftools view --no-version -T ${split_file} -Oz ${vcf} > ${prefix}.${split_file}.vcf.gz
  bcftools index ${index_flag} ${prefix}.${split_file}.vcf.gz
  """
}
