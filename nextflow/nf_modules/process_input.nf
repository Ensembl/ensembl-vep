#!/usr/bin/env nextflow

/* 
 * Script to process inputs channels
 */

nextflow.enable.dsl=2

// defaults
params.cpus = 1

process processInput {
  /*
  Generate input for subsequent jobs. It works like merge operator and it helps creating one-to-one or many-to-one or many-to-one relationship between vcf file and vep_config

  Returns
  -------
  Returns 2 files:
      1) A VCF file
      2) A tabix index for that VCF
      3) A VEP config file to run VEP on that VCF
  */
  cpus params.cpus

  input:
  path vcf
  path vep_config
  val output_dir

  output:
  tuple path(vcf), env(vcf_index), path(vep_config), val(output_dir)
  
  script:
  """
  vcf_index=`readlink -f ${vcf}`
  vcf_index=\${vcf_index}.tbi
  """
}
