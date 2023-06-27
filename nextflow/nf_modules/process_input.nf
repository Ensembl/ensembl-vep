#!/usr/bin/env nextflow

/* 
 * Script to process inputs channels
 */
 
nextflow.enable.dsl=2

process processInput {
  /*
  Generate input for subsequent jobs. It works like merge operator and it helps creating one-to-one or many-to-one or many-to-one relationship between vcf file, vep_config, and output_dir

  Returns
  -------
  Returns:
      1) A VCF file
      2) A tabix index for that VCF
      3) A VEP config file
      4) A output dir
      5) The type of tabix index, either tbi or csi
  */
  cpus params.cpus

  input:
  path vcf
  path vep_config
  val output_dir

  output:
  tuple path(vcf), env(vcf_index), path(vep_config), val(output_dir), env(index_type)
  
  script:
  """
  vcf_filepath=`readlink -f ${vcf}`
  if [[ -f \${vcf_filepath}.csi ]]; then
    vcf_index=\${vcf_filepath}.csi
    index_type=csi
  else
    vcf_index=\${vcf_filepath}.tbi
    index_type=tbi
  fi
  """
}
