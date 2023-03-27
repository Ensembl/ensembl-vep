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
  label 'bcftools'

  input:
  val(chr)
  path(vcf)
  path(vcf_index)
  val(bin_size)

  output:
  tuple path("${prefix}.${chr}.*vcf.gz"), path("${prefix}.${chr}.*vcf.gz.tbi"), emit: files

  script:
  """
  bcftools view --no-version -r ${chr} ${vcf} -o ${prefix}.${chr}.vcf.gz -O z
  bcftools index -t ${prefix}.${chr}.vcf.gz
  
  if [[ ${bin_size} ]]; then 
    bcftools query -f'%CHROM\t%POS\n' ${prefix}.${chr}.vcf.gz | split -l ${bin_size}
    
    for file in x*; do 
      bcftools view --no-version -T \${file} -Oz ${prefix}.${chr}.vcf.gz > ${prefix}.${chr}.\${file}.vcf.gz
      bcftools index -t ${prefix}.${chr}.\${file}.vcf.gz
    done
    
    rm ${prefix}.${chr}.vcf.gz
    rm ${prefix}.${chr}.vcf.gz.tbi
    rm x*
  fi
  """
}
