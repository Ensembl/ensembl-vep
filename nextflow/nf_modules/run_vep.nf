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
  tuple val(meta), val(original_vcf), path(vcf), path(vcf_index), path(vep_config)
  
  output:
  tuple val(meta), val(original_vcf), path(out_vcf), path("${out_vcf}.{tbi,csi}"), val("${vep_config}"), emit: files

  script:
  index_type = meta.index_type
  out_vcf = "vep" + "-" + file(original_vcf).getSimpleName() + "-" + vep_config.getSimpleName() + "-" + vcf
  tabix_arg = index_type == 'tbi' ? '' : '-C'
  
  if( !vcf.exists() ) {
    exit 1, "VCF file is not generated: ${vcf}"
  }
  else if ( !vcf_index.exists() ){
    exit 1, "VCF index file is not generated: ${vcf_index}"
  }
  else if ( meta.filters != null ){
    def filters = meta.filters.split(",")
    def filter_arg = ""
    for (filter in filters) {
      filter_arg = filter_arg + "-filter \"" + filter + "\" "
    }
    """
    vep -i ${vcf} -o STDOUT --vcf --format vcf --config ${vep_config} | filter_vep -o filtered.vcf --only_matched --format vcf ${filter_arg}  
    
    bgzip filtered.vcf
    mv filtered.vcf.gz ${out_vcf}
    
    tabix ${tabix_arg} -p vcf ${out_vcf}
    """
  }
  else {
    """
    vep -i ${vcf} -o ${out_vcf} --vcf --compress_output bgzip --format vcf --config ${vep_config}
    
    tabix ${tabix_arg} -p vcf ${out_vcf}
    """
  }
}
