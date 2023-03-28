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
params.vcf_dir = null

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
  --vcf VCF                 VCF that will be split. Currently supports sorted and bgzipped file
  --vcf_dir                 VCF directory Directory containing input VCF files
  --bin_size INT            Input file is split into multiple files with a given number of variants. Enables faster run in expense of more jobs. Default: 100
  --vep_config FILENAME     VEP config file. Default: vep_config/vep.ini
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output dir. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>.vcf.gz
  --skip_check [0,1]        Skip checking for tabix index file of input VCF. Enables the first module to load from cache if -resume is used. Default: 0
  """
  exit 1
}

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
if( !params.vcf && !params.vcf_dir) {
  exit 1, "Undefined input parameters. Please provide the path to a VCF file or the directory containing VCF files"
}

// Cannot use both --vcf and --vcf_dir
else if ( params.vcf && params.vcf_dir) {
  exit 1, "Please specify one input, --vcf or --vcf_dir"
}

log.info 'Starting workflow.....'

workflow {
  input = params.vcf ? params.vcf : params.vcf_dir 
  vcfInput = file(input)
  if( !vcfInput.exists() ) {
    exit 1, "The specified VCF input does not exist: ${vcfInput}"
  }
  checkVCF(Channel.fromPath("${input}*.gz"), params.vep_config)
  splitVCF(checkVCF.out)
  runVEP(splitVCF.out.files.transpose(), params.vep_config)
  mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
}  
