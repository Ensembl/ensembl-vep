#!/usr/bin/env nextflow

/* 
 * Read chromosomes 
 */

nextflow.enable.dsl=2

prefix = ""

process generateSplits {
  /*
  Function to read variants from VCF file and write them to separate to split files

  Returns
  -------
  Tuple of VCF, VCF index, split files, vep config file, a output dir, and the index type of VCF file
  */
  
  cpus params.cpus
  label 'bcftools'

  input:
  tuple val(meta), path(vcf), path(vcf_index), path(vep_config)

  output:
  tuple val(meta), path(vcf), path(vcf_index), path("x*"), path(vep_config)

  shell:
  """
  bcftools query -f'%CHROM\t%POS\n' ${vcf} | uniq | split -a 3 -l ${params.bin_size}
  """
}
