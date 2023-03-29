#!/usr/bin/env nextflow

/* 
 * Script to check if the files are bgzipped and bgzip if not
 */

nextflow.enable.dsl=2

// defaults

params.cpus = 1

process checkVCF {
  /*
  Function to check input VCF files

  Returns
  -------
  Tuple of VCF, VCF index and GFF file
  */

  cpus params.cpus
  label 'vep'
  errorStrategy 'ignore'

  input:
  path input_vcf
  path input_vcf_index
  
  output:
  tuple path("*.gz", includeInputs: true), path ("*.gz.tbi", includeInputs: true)

  afterScript "rm *.vcf *.vcf.tbi"

  script:
  """
  bgzip -t ${input_vcf} || bgzip -c ${input_vcf} > ${input_vcf}.gz
  [ -f *gz.tbi ] || tabix -p vcf -f *.gz
  """

}
