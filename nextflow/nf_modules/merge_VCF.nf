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
  */
    
  cpus params.cpus
  label 'bcftools'
  cache 'lenient'
   
  input:
  tuple val(original), path(vcfFiles), path(indexFiles), val(output_dir)
  
  output:
  val("${output_dir}/${ original }_vep_${ mergedVCF }.vcf.gz")

  script: 
  """
  mkdir -p temp
  out=${ original }_vep_${ mergedVCF }.vcf.gz
  sorted_vcf=\$(echo ${vcfFiles} | xargs -n1 | sort | xargs)
  bcftools concat --no-version -a \${sorted_vcf} -Oz -o \${out}
  bcftools index -t \${out}
  
  mv \${out}* ${ output_dir }
  """
}
