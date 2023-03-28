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
params.output_prefix = "out"
params.bin_size = 100
params.skip_check = 0
params.help = false

// module imports
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
  --bin_size INT            Input file is split into multiple files with a given number of variants. Enables faster run in expense of more jobs. Default: 100
  --vep_config FILENAME     VEP config file. Default: vep_config/vep.ini
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output dir. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>.vcf.gz
  --skip_check [0,1]        Skip checking for tabix index file of input VCF. Enables the first module to load from cache if -resume is used. Default: 0
  """
  exit 1
}

// Input validation

if( !params.vcf) {
  exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file"
}

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified VCF file does not exist: ${params.vcf}"
}

check_bgzipped = "bgzip -t $params.vcf".execute()
check_bgzipped.waitFor()
if(check_bgzipped.exitValue()){
  exit 1, "The specified VCF file is not bgzipped: ${params.vcf}"
}

if ( !params.skip_check ){
  def sout = new StringBuilder(), serr = new StringBuilder()
  check_parsing = "tabix -p vcf -f $params.vcf".execute()
  check_parsing.consumeProcessOutput(sout, serr)
  check_parsing.waitFor()
  if( serr ){
    exit 1, "The specified VCF file has issues in parsing: $serr"
  }
}
vcf_index = "${params.vcf}.tbi"

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

log.info 'Starting workflow.....'

workflow {
  splitVCF(params.vcf, vcf_index, params.bin_size)
  runVEP(splitVCF.out.files.transpose(), params.vep_config)
  mergeVCF(runVEP.out.vcfFile.collect(), runVEP.out.indexFile.collect())
}  
