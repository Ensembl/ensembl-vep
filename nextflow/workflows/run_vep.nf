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
params.vep_config = ""
params.chros = ""
params.chros_file = ""
params.skip_check = 0

// module imports
include { splitVCF } from '../nf_modules/split_VCF.nf' 
include { mergeVCF } from '../nf_modules/merge_VCF.nf'  
include { chrosVEP } from '../nf_modules/run_vep.nf'
include { readChrVCF } from '../nf_modules/read_VCF.nf'

 // print usage
if (params.help) {
  log.info """
Pipeline to run VEP
-------------------

Usage:
  nextflow run workflows/run_vep.nf --vcf <path-to-vcf> --vep_config vep_config/vep.ini

Options:
  --vcf VCF                 VCF that will be split. Currently supports sorted and bgzipped file
  --outdir DIRNAME          Name of output dir. Default: outdir
  --vep_config FILENAME     VEP config file. Default: vep_config/vep.ini
  --chros LIST_OF_CHROS     Comma-separated list of chromosomes to generate. i.e. 1,2,..., Default: 1,2,...X,Y,MT
  --chros_file FILE         Path to file containing list of chromosomes
  --cpus INT                Number of CPUs to use. Default 1.
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>.vcf.gz
  --skip_check [0,1]        Skip checking for tabix index file of input VCF. Enables the first module to load from cache if -resume is used. Default: 0
  --bin_size SIZE           Each chromosome file will contain at most SIZE number of variants. Enables faster run in expense of more jobs. By default, the input is only split by chromosome (SIZE=0). Default: 0
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
log.info params.chros
  if (params.chros){
    log.info 'Reading chromosome names from list'
    chr_str = params.chros.toString()
    chr = Channel.of(chr_str.split(','))
  }
  else if (params.chros_file) {
    log.info 'Reading chromosome names from file'
    chr = Channel.fromPath(params.chros_file).splitText().map{it -> it.trim()}
  }
  else {
    log.info 'Computing chromosome names from input'
    readChrVCF(params.vcf, vcf_index)
    chr = readChrVCF.out.splitText().map{it -> it.trim()}
  }
  splitVCF(chr, params.vcf, vcf_index, params.bin_size)
  chrosVEP(splitVCF.out.files.transpose(), params.vep_config)
  mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
}  
