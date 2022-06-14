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
params.vep_config=""
params.singularity_dir=""

// module imports
include { splitVCF } from '../nf_modules/split_into_chros.nf' 
include { mergeVCF } from '../nf_modules/merge_chros_VCF.nf'  
include { chrosVEP } from '../nf_modules/run_vep_chros.nf'

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
  log.info '  --outdir DIRNAME          Name of output dir. Default: outdir'
  log.info '  --vep_config FILENAME     VEP config file. Default: nf_config/vep.ini'
  log.info '  --chros LIST_OF_CHROS	Comma-separated list of chromosomes to generate. i.e. 1,2,... Default: 1,2,...,X,Y,MT'
  log.info '  --cpus INT	        Number of CPUs to use. Default 1.'
  exit 1
}

// Input validation
if( !params.chros) {
  exit 1, "Undefined --chros parameter. Please provide a comma-separated string with chromosomes: 1,2, ... Default: 1,2,...,X,Y,MT"
}
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

def sout = new StringBuilder(), serr = new StringBuilder()
check_parsing = "$params.singularity_dir/vep.sif tabix -p vcf -f $params.vcf".execute()
check_parsing.consumeProcessOutput(sout, serr)
check_parsing.waitFor()
if( serr ){
  exit 1, "The specified VCF file has issues in parsing: $serr"
}
vcf_index = "${params.vcf}.tbi"

vepFile = file(params.vep_config)
if( !vepFile.exists() ) {
  exit 1, "The specified VEP config does not exist: ${params.vep_config}"
}

log.info 'Starting workflow.....'

workflow {
  chr_str = params.chros.toString()
  chr = Channel.of(chr_str.split(','))
  splitVCF(chr, params.vcf, vcf_index)
  chrosVEP(splitVCF.out, params.vep_config)
  mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
}  
