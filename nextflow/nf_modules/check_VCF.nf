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
  path(input_vcf)
  
  output:
  tuple path("*.gz", includeInputs: true), path ("${input_vcf}*.tbi")

  afterScript "rm *.vcf"

  script:
  """
  bgzip -t ${input_vcf} || bgzip -c ${input_vcf} > ${input_vcf}.gz
  tabix -p vcf -f *.gz
  """

}
