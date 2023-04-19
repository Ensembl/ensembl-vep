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
  Function to read variants from VCF file and write them to separate to split files

  Returns
  -------
  Tuple of VCF, VCF index, split files, vep config file and a output dir
  */
  
  cpus params.cpus
  label 'bcftools'

  input:
  tuple path(vcf), path(vcf_index), path(vep_config), val(output_dir)
  val(bin_size)

  output:
  tuple path(vcf), path(vcf_index), path("x*"), path(vep_config), val(output_dir)

  shell:
  """
  bcftools query -f'%CHROM\t%POS\n' ${vcf} | uniq | split -l ${bin_size}
  """
}
