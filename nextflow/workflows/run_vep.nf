/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

nextflow.enable.dsl=2

 // params default
params.cpus = 1
params.vep_config = "$PWD/vep_config/vep.ini"
params.outdir = "outdir"

params.vcf = null

params.output_prefix = "out"
params.bin_size = 100
params.skip_check = 0
params.help = false

// module imports
include { processInput } from '../nf_modules/process_input.nf'
include { checkVCF } from '../nf_modules/check_VCF.nf'
include { generateSplits } from '../nf_modules/generate_splits.nf'
include { splitVCF } from '../nf_modules/split_VCF.nf' 
include { mergeVCF } from '../nf_modules/merge_VCF.nf'  
include { runVEP } from '../nf_modules/run_vep.nf'

// print usage
if (params.help) {
  log.info """
Pipeline to run VEP
-------------------

Usage:
  nextflow run workflows/run_vep.nf --vcf <path-to-vcf> --vep_config vep_config/vep.ini

Options:
  --vcf VCF                 Sorted and bgzipped VCF. Alternatively, can also be a directory containing VCF files
  --bin_size INT            Number of variants used to split input VCF into multiple jobs. Default: 100
  --vep_config FILENAME     VEP config file. Alternatively, can also be a directory containing VEP INI files. Default: vep_config/vep.ini
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output directory. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <vcf>-<output_prefix>.vcf.gz
  --skip_check [0,1]        Skip check for tabix index file of input VCF. Enables use of cache with -resume. Default: 0
  """
  exit 1
}


def createChannels (input, pattern, is_vcf) {
  if (input instanceof String) {
    // if input is a String, process as file or directory
    files = file(input)
    if ( !files.exists() ) {
      exit 1, "The specified input does not exist: ${input}"
    }

    if (files.isDirectory()) {
      files = "${files}/${pattern}"
    }
    files = Channel.fromPath(files)
  } else {
    // if input is a Channel, just pass along
    files = input
  }
  
  return files;
}

workflow vep {
  take:
    vcf
    vep_config
  main:
    if (!vcf) {
      exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file"
    }

    if (!vep_config) {
      exit 1, "Undefined --vep_config parameter. Please provide a VEP config file"
    }
    
    vcf = createChannels(vcf, pattern="*.{vcf,gz}", true)
    vep_config = createChannels(vep_config, pattern="*.ini", false)
    
    vcf.count()
      .combine( vep_config.count() )
      .subscribe{ if ( it[0] != it[1] && it[0] != 1 && it[1] != 1 ) 
        exit 1, "Cannot map VCF and VEP config files to one-to-one, many-to-one, or, one-to-many scenario" 
      }
        
    // process input and create Channel
    // this works like 'merge' operator and thus might make the pipeline un-resumable
    // we might think of using 'toSortedList' and generate appropriate input from the 'processInput' module
    processInput(vcf, vep_config)
    
    // Prepare input VCF files (bgzip + tabix)
    checkVCF(processInput.out)
    
    // Generate split files that each contain bin_size number of variants from VCF
    generateSplits(checkVCF.out, params.bin_size)

    // Split VCF using split files
    splitVCF(generateSplits.out.transpose())

    // Run VEP for each split VCF file and for each VEP config
    runVEP(splitVCF.out.transpose())

    // Merge split VCF files (creates one output VCF for each input VCF)
    mergeVCF(runVEP.out.vcf.groupTuple())
  emit:
    mergeVCF.out
}

workflow {
  vep(params.vcf, params.vep_config)
}
