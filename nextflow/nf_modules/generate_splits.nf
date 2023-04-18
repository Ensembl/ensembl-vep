#!/usr/bin/env nextflow

/* 
 * Read chromosomes 
 */

nextflow.enable.dsl=2

// defaults
prefix = ""
params.outdir = ""
params.cpus = 1


process generateSplits {
  /*
  Function to read variants from VCF file to split files

  Returns
  -------
  Returns:
      1) A VCF file
      2) A tabix index for that VCF
      3) A list of split files
  */
  cpus params.cpus
  label 'bcftools'

  input:
  tuple path(vcf), path(vcf_index), path(vep_config)
  val(bin_size)

  output:
  tuple path(vcf), path(vcf_index), path("x*"), path(vep_config)

  shell:
  """
  bcftools query -f'%CHROM\t%POS\n' ${vcf} | uniq | split -l ${bin_size}
  """
}
