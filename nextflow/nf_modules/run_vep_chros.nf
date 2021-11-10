#!/usr/bin/env nextflow

/* 
 * Script to run VEP on chromosome-wise split VCF files
 */

nextflow.enable.dsl=2

// defaults
prefix = "vep"
params.outdir = ""
params.cpus = 1

process chrosVEP {
  /*
  Function to run VEP on chromosome-wise split VCF files

  Returns
  -------
  Returns 2 files per chromosome:
      1) VEP output file for each chromosome-wise split VCF
      2) A tabix index for that VCF output file
  */
  cpus params.cpus

  input:
  tuple path(vcfFile), path(indexFile)
  path(vep_config)

  output:
  path("${prefix}-*.vcf.gz"), emit: vcfFile
  path("${prefix}-*.vcf.gz.tbi"), emit: indexFile

  script:
  """
  vep -i ${vcfFile} -o ${prefix}-${vcfFile} --vcf --compress_output bgzip --format vcf --config ${vep_config} 
  tabix -p vcf ${prefix}-${vcfFile}
  """	
}
