#!/usr/bin/env nextflow

/*
 * Script to merge VCF files into a single file
 */

nextflow.enable.dsl=2

process mergeVCF {
  /*
  Merge VCF files into a single file
  */

  cache 'lenient'
  cpus params.cpus
  label 'bcftools'
  cache 'lenient'

  input:
  tuple val(meta), val(original_file), path(vcf_files), path(index_files), val(vep_config)

  output:
  tuple val(file_base_name), val(vep_config), val("${output_dir}/${merged_vcf}")

  script:
  file_base_name = meta.file_base_name
  merged_vcf = null
  if ( params.output_prefix ){
    merged_vcf = params.output_prefix
  }
  else{
    merged_vcf = meta.file_base_name
  }

  if (meta.one_to_many) {
    merged_vcf += "_" + file(vep_config).getName().replace(".ini", "")
  }

  merged_vcf += '_VEP.vcf.gz'

  def index_type = meta.index_type
  output_dir = meta.output_dir

  def index_flag = index_type == "tbi" ? "-t" : "-c"

  """
  sorted_vcfs=\$(echo ${vcf_files} | xargs -n1 | sort | xargs)
  bcftools concat --no-version --naive \${sorted_vcfs} -Oz -o ${merged_vcf}
  bcftools index ${index_flag} ${merged_vcf}

  # move the output file
  mkdir -p ${output_dir}
  mv ${merged_vcf}* ${output_dir}
  """
}
