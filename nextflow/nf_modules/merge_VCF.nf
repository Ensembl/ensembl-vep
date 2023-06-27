#!/usr/bin/env nextflow

/* 
 * Script to merge VCF files into a single file
 */

nextflow.enable.dsl=2

// defaults
merged_vcf = null
if ( params.output_prefix != "" ){
  merged_vcf = params.output_prefix + "_VEP.vcf.gz"
}

process mergeVCF {
  /*
  Merge VCF files into a single file
  */
    
  cpus params.cpus
  label 'bcftools'
  cache 'lenient'
   
  input:
  tuple val(original_vcf), path(vcf_files), path(index_files), val(output_dir), val(index_type)
  
  output:
  val("${output_dir}/${merged_vcf}")

  script:
  merged_vcf = merged_vcf ?: file(original_vcf).getName().replace(".vcf", "_VEP.vcf")
  index_flag = ${index_type} == "tbi" ? "-t" : "-c:
  
  """
  sorted_vcfs=\$(echo ${vcf_files} | xargs -n1 | sort | xargs)
  bcftools concat --no-version -a \${sorted_vcfs} -Oz -o ${merged_vcf}
  bcftools index ${index_flag} ${merged_vcf}
  
  # move the output file
  mkdir -p ${output_dir}
  mv ${merged_vcf}* ${output_dir}
  """
}
