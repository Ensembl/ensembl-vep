/* 
 * Workflow to run VEP on chromosome based VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

nextflow.enable.dsl=2

 // params default
params.help = false
params.cpus = 1
params.outdir = "outdir"
params.singularity_dir=""
params.vep_config=""
params.chros=""
params.chros_file=""
params.vcf_dir = ""
params.vcf = ""

// module imports
include { splitVCF } from '../nf_modules/split_into_chros.nf' 
include { mergeVCF } from '../nf_modules/merge_chros_VCF.nf'  
include { chrosVEP } from '../nf_modules/run_vep_chros.nf'
include { readChrVCF } from '../nf_modules/read_chros_VCF.nf'
include { checkVCF } from '../nf_modules/check_VCF.nf'

 // print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP chromosome-wise'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C nf_config/nextflow.config run workflows/run_vep.nf --vcf <path-to-vcf> --chros 1,2 --vep_config'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf VCF                 VCF that will be split. Currently supports sorted and bgzipped file'
  log.info '  --vcf_dir                 VCF directory Directory containing input VCF files'
  log.info '  --outdir DIRNAME          Name of output dir. Default: outdir'
  log.info '  --vep_config FILENAME     VEP config file. Default: nf_config/vep.ini'
  log.info '  --chros LIST_OF_CHROS	Comma-separated list of chromosomes to generate. i.e. 1,2,... Default: 1,2,...,X,Y,MT'
  log.info '  --chros_file LIST_OF_CHROS_FILE Path to file containing list of chromosomes' 
  log.info '  --cpus INT	        Number of CPUs to use. Default 1.'
  log.info '  --output_prefix FILENAME_PREFIX   Output filename prefix. The generated output file will have name <output_prefix>.vcf.gz'
  exit 1
}
//Input validation

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
  chrosVEP(splitVCF.out)
  mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
}  
