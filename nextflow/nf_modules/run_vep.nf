#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
 */

nextflow.enable.dsl=2

// defaults
prefix = "vep"
params.outdir = ""
params.cpus = 1

process runVEP {
  /*
  Run VEP on VCF files

  Returns
  -------
  Returns 2 files:
      1) VEP output file in VCF format
      2) A tabix index for that VCF output file
  */
  publishDir "${params.outdir}/vep-summary",
    pattern: "${prefix}-*.vcf.gz_summary.*",
    mode:'move'
  cpus params.cpus
  label 'vep'

  input:
  tuple val(original), path(vcfFile), path(indexFile)
  path(vep_config)
  
  output:
  tuple val(original), path("${prefix}-*.vcf.gz"), path("${prefix}-*.vcf.gz.tbi"), emit: vcf
  path("${prefix}-*.vcf.gz_summary.*")

  script:
  if( !vcfFile.exists() ) {
    exit 1, "VCF file is not generated: ${vcfFile}"
  }
  else if ( !indexFile.exists() ){
    exit 1, "VCF index file is not generated: ${indexFile}"
  }
  else {
    """
    vep -i ${vcfFile} -o ${prefix}-${vcfFile} --vcf --compress_output bgzip --format vcf --config ${vep_config}
    tabix -p vcf ${prefix}-${vcfFile}
    """	
  }
}
