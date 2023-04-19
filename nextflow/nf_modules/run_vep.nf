#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
 */

nextflow.enable.dsl=2

// defaults
params.outdir = ""
params.cpus = 1

process runVEP {
  /*
  Run VEP on VCF files

  Returns
  -------
  Tuple of original VCF, split VCF file after running VEP, VCF index, vep config file and a output dir
  */
  
  publishDir "${params.outdir}/vep-summary",
    pattern: "vep-${original}-${vep_config}-*.vcf.gz_summary.*",
    mode:'move'
  cpus params.cpus
  label 'vep'

  input:
  tuple val(original), path(vcfFile), path(indexFile), path(vep_config), val(output_dir)
  
  output:
  tuple val(original), path("*.vcf.gz"), path("*.vcf.gz.tbi"), val(output_dir), emit: vcf
  path("*.vcf.gz_summary.*")

  script:
  if( !vcfFile.exists() ) {
    exit 1, "VCF file is not generated: ${vcfFile}"
  }
  else if ( !indexFile.exists() ){
    exit 1, "VCF index file is not generated: ${indexFile}"
  }
  else {
    """
    out=vep-${original}-${vep_config}-${vcfFile}
    vep -i ${vcfFile} -o \${out} --vcf --compress_output bgzip --format vcf --config ${vep_config}
    tabix -p vcf \${out}
    """	
  }
}
