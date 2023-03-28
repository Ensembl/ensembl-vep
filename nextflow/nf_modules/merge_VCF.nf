#!/usr/bin/env nextflow

/* 
 * Script to merge VCF files into a single file
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
mergedVCF = "merged-file"
if ( params.output_prefix != "" ){
  mergedVCF = params.output_prefix
}
params.outdir = ""
params.cpus = 1

process mergeVCF {
  /*
  Merge VCF files into a single file

  Returns
  -------
  Returns 2 files:
      1) A VCF format file 
      2) A tabix index for that VCF
  */

  publishDir "${params.outdir}", 
    enabled: "${params.outdir}" as Boolean,
    mode:'move'
    
  cpus params.cpus
  label 'bcftools'
  cache 'lenient'
   
  input:
  tuple val(original), path(vcfFiles), path(indexFiles)

  output:
  path("${ original }_vep_${ mergedVCF }.vcf.gz*")

  script: 
  """
  mkdir -p temp
  bcftools concat --no-version -a ${ vcfFiles } -Oz -o temp-${ mergedVCF}.vcf.gz

  out=${ original }_vep_${ mergedVCF }.vcf.gz
  bcftools sort -T temp -Oz temp-${ mergedVCF }.vcf.gz -o \${out}
  bcftools index -t \${out}
  """
}
