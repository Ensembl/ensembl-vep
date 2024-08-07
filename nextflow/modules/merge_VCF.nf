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
  tuple val(meta), val(original_file), path(vcf_files), path(index_files), val(vep_config)
  
  output:
  val("${output_dir}/${merged_vcf}")

  script:
  index_type = meta.index_type
  one_to_many = meta.one_to_many
  output_dir = meta.output_dir
  
  merged_vcf = merged_vcf ?: file(original_file).getSimpleName() + "_VEP.vcf.gz"
  merged_vcf = one_to_many ? merged_vcf.replace(
    "_VEP.vcf", 
    "_" + file(vep_config).getName().replace(".ini", "") + "_VEP.vcf"
  ) : merged_vcf
  index_flag = index_type == "tbi" ? "-t" : "-c"
  
  """
  sorted_vcfs=\$(echo ${vcf_files} | xargs -n1 | sort | xargs)
  bcftools concat --no-version --naive \${sorted_vcfs} -Oz -o ${merged_vcf}
  bcftools index ${index_flag} ${merged_vcf}
  
  # move the output file
  mkdir -p ${output_dir}
  mv ${merged_vcf}* ${output_dir}
  """
}
