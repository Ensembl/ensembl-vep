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
include { checkVCF } from '../nf_modules/check_VCF.nf'
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
  --vcf VCF                 Sorted and bgzipped VCF. Alternatively, can also be a direcotry containing VCF files
  --bin_size INT            Input file is split into multiple files with a given number of variants. Enables faster run in expense of more jobs. Default: 100
  --vep_config FILENAME     VEP config file. Default: vep_config/vep.ini
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output dir. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>.vcf.gz
  --skip_check [0,1]        Skip checking for tabix index file of input VCF. Enables the first module to load from cache if -resume is used. Default: 0
  """
  exit 1
}


log.info 'Starting workflow.....'

workflow vep {
  take:
    vcf
    vep_config
  main:
    // VEP config required
    if ( params.vep_config ){
      vepFile = file(params.vep_config)
      if( !vepFile.exists() ){
        exit 1, "The specified VEP config does not exist: ${params.vep_config}"
      }
    }
    else
    {
      exit 1, "Undefined --vep_config parameter. Please provide a VEP config file"
    }

    // Input required
    if( !params.vcf) {
      exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file"
    }

    if(vcf instanceof String){
      vcfInput = file(vcf)
      if( !vcfInput.exists() ) {
        exit 1, "The specified VCF input does not exist: ${vcfInput}"
      }

      if (vcfInput.isDirectory()) {
        t=Channel.fromPath("${vcfInput}/{*.gz,*.vcf}").multiMap{ it ->
        vcf: it 
        vcf_index: "${it}.tbi"
        }
        checkVCF(t)
      } else {
        checkVCF(vcfInput,"${vcfInput}.tbi")
      }
    }
    else {
      //TODO 
      // checkVCF(vcf)
    }

    splitVCF(checkVCF.out, params.bin_size)
    vep_config = Channel.fromPath(vep_config, relative: true).first()
    runVEP(splitVCF.out.files.transpose(), vep_config)
    mergeVCF(runVEP.out.vcf.groupTuple())
  emit:
    mergeVCF.out
}

workflow {
  vep(params.vcf, params.vep_config)
}
