#!/usr/bin/env nextflow

/* 
 * Script to run VEP on split VCF files
 */

nextflow.enable.dsl=2

process runVEP {
  /*
  Run VEP on VCF files

  Returns
  -------
  Tuple of original VCF, split VCF file after running VEP, tabix index of that file, vep config file, a output dir, and the index type of VCF file
  */
  
  publishDir "${params.outdir}/vep-summary",
    pattern: "vep-${original}-${vep_config}-*.vcf.gz_summary.*",
    mode:'move'
  cpus params.cpus
  label 'vep'

  input:
  tuple val(original_vcf), path(vcf), path(vcf_index), path(vep_config), val(index_type)
  
  output:
  tuple val(original_vcf), path(out_vcf), path("${out_vcf}.{tbi,csi}"), val("${vep_config}"), val(index_type), emit: files
  path("*.vcf.gz_summary.*")

  script:
  out_vcf = "vep" + "-" + file(original_vcf).getSimpleName() + "-" + vep_config.getSimpleName() + "-" + vcf
  
  if( !vcf.exists() ) {
    exit 1, "VCF file is not generated: ${vcf}"
  }
  else if ( !vcf_index.exists() ){
    exit 1, "VCF index file is not generated: ${vcf_index}"
  }
  else {
    """
    vep -i ${vcf} -o ${out_vcf} --vcf --compress_output bgzip --format vcf --config ${vep_config}
    
    if [[ "${index_type}" == "tbi" ]]; then
      tabix -p vcf ${out_vcf}
    else
      tabix -C -p vcf ${out_vcf}
    fi
    """
  }
}
