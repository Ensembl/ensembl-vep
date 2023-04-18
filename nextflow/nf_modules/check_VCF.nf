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
  tuple path(input_vcf), path(input_vcf_index), path(vep_config)
  
  output:
  tuple path("*.gz", includeInputs: true), path ("*.gz.tbi", includeInputs: true), path(vep_config)

  afterScript "rm *.vcf *.vcf.tbi"

  script:
  """
  [ -f *gz ] || bgzip -c ${input_vcf} > ${input_vcf}.gz
  [ -f *gz.tbi ] || tabix -p vcf -f *.gz

  # quickly test tabix -- ensures both bgzip and tabix are okay
  chr=\$(tabix -l *.gz | head -n1)
  tabix *.gz \${chr}:1-10001
  """
}
