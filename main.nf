#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NF_VEP } from './nextflow/workflows/run_vep'

workflow {
  NF_VEP()
}
