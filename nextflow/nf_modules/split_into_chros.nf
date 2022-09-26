#!/usr/bin/env nextflow

/* 
 * Script to split a multi-chromosome VCF into single-chromosome VCFs
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
params.outdir = ""
params.cpus = 1

process splitVCF {
  /*
  Function to split a multi-chromosome VCF into single chromosome VCF

  Returns
  -------
  Returns 2 files per chromosome:
      1) A VCF format file for each splitted chromosome
      2) A tabix index for that VCF
  */
  cpus params.cpus
  container "${params.singularity_dir}/bcftools.sif"

  input:
  val(chr)
  path(vcf)
  path(vcf_index)

  output:
  tuple path("${prefix}.${chr}.vcf.gz"), path("${prefix}.${chr}.vcf.gz.tbi")

  script:
  """
  bcftools view -r ${chr} ${vcf} -o ${prefix}.${chr}.vcf.gz -O z
  bcftools index -t ${prefix}.${chr}.vcf.gz
  """
}
