#!/usr/bin/env nextflow

/* 
 * Read chromosomes 
 */

nextflow.enable.dsl=2

// defaults
prefix = ""
params.outdir = ""
params.cpus = 1


process readChrVCF {
  /*
  Function to read chr from VCF file

  Returns
  -------
  Returns list of chromosomes
  */
  cpus params.cpus
  container "${params.singularity_dir}/bcftools.sif"

  input:
  path(vcf)
  path(vcf_index)

  output:
  path("scaffolds.txt")

  shell:
  """
  tabix ${vcf} -l > scaffolds.txt 
  """
}
