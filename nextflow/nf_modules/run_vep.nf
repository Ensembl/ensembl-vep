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
    pattern: "vep-${original}-${vep_config}-*.gz_summary.*",
    mode:'move'
  cpus params.cpus
  label 'vep'

  input:
  tuple val(meta), val(original_file), path(input), path(index), path(vep_config), val(format)
  
  output:
  tuple val(meta), val(original_file), path("${out}{.gz,}"), path("${out}{.gz,}.{tbi,csi}"), val("${vep_config}"), emit: files

  script:
  index_type = meta.index_type
  out = "vep" + "-" + file(original_file).getSimpleName() + "-" + vep_config.getSimpleName() + "-" + input.getName().replace(".gz", "")
  tabix_arg = index_type == 'tbi' ? '' : '-C'
  
  if( !input.exists() ) {
    exit 1, "VCF file is not generated: ${input}"
  }
  else if ( format == 'vcf' && !index.exists() ){
    exit 1, "VCF index file is not generated: ${index}"
  }
  else if ( meta.filters != null ){
    def filters = meta.filters.split(",")
    def filter_arg = ""
    for (filter in filters) {
      filter_arg = filter_arg + "-filter \"" + filter + "\" "
    }
    """
    vep -i ${input} -o STDOUT --vcf --config ${vep_config} | filter_vep -o filtered.vcf --only_matched ${filter_arg}

    # Sort, bgzip and tabix VCF
    (head -1000 filtered.vcf | grep "^#"; cat out.vcf | grep -v "^#" | sort -k1,1d -k2,2n) > sorted.vcf
    mv sorted.vcf ${out}
    bgzip ${out}
    tabix ${tabix_arg} ${out}.gz
    """
  }
  else {
    """
    vep -i ${input} -o out.vcf --vcf --config ${vep_config}
    
    # Sort, bgzip and tabix VCF
    (head -1000 out.vcf | grep "^#"; cat out.vcf | grep -v "^#" | sort -k1,1d -k2,2n) > sorted.vcf
    mv sorted.vcf ${out}
    bgzip ${out}
    tabix ${tabix_arg} ${out}.gz
    """
  }
}
