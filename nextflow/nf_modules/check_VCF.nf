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
  container "${params.singularity_dir}/vep.sif"
  errorStrategy 'ignore'

  input:
  path(input_vcf)
  path(vep_config)
  
  output:
  tuple (path("${input_vcf}"), path ("${input_vcf}.tbi"),path("${vep_config}") )

  script:
  """
  bgzip -t ${ input_vcf}
  tabix -p vcf -f ${ input_vcf}

  """

}
