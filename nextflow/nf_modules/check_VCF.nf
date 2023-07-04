#!/usr/bin/env nextflow

/* 
 * Script to check if the files are bgzipped and bgzip if not
 */

nextflow.enable.dsl=2

process checkVCF {
  /*
  Function to check input VCF files

  Returns
  -------
  Tuple of VCF, VCF index, vep config file, a output dir, and the index type of VCF file
  */

  cpus params.cpus
  label 'vep'
  errorStrategy 'ignore'

  input:
  tuple path(vcf), path(vcf_index), path(vep_config), val(index_type)
  
  output:
  tuple path("*.gz", includeInputs: true), path ("*.gz.{tbi,csi}", includeInputs: true), path(vep_config), val(index_type)

  afterScript "rm *.vcf *.vcf.tbi *.vcf.csi"

  script:
  """
  [ -f *gz ] || bgzip -c ${vcf} > ${vcf}.gz
  
  if [[ "${index_type}" == "tbi" ]]; then
    [ -f *gz.tbi ] || tabix -p vcf -f *.gz
  else
    [ -f *gz.csi ] || tabix -C -p vcf -f *.gz
  fi
    
  # quickly test tabix -- ensures both bgzip and tabix are okay
  chr=\$(tabix -l *.gz | head -n1)
  tabix *.gz \${chr}:1-10001
  """
}
