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
  log.info '  nextflow -C nf_config/run_vep.config run workflows/run_vep.nf --vcf /mnt/d/work/repos/nf_genomics_pipelines/test-data/test.vcf.gz --chros 1,2 --vep_config /path/to/vep.ini'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf VCF                 VCF that will be split.'
  log.info '  --outdir DIRNAME          Name of output dir. Default: outdir'
  log.info '  --vep_config FILENAME     VEP config file. Default: nf_config/vep.ini'
  log.info '  --chros LIST_OF_CHROS	    Comma-separated list of chromosomes to generate. i.e. chr1,chr2,...'
  log.info '  --cpus INT	              Number of CPUs to use. Default 1.'
  exit 1
}

// Input validation
if( !params.chros) {
  exit 1, "Undefined --chros parameter. Please provide a comma-separated string with chromosomes: chr1,chr2, ..."
}
if( !params.vcf) {
  exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file"
}

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified VCF file does not exist: ${params.vcf}"
}

vepFile = file(params.vep_config)
if( !vepFile.exists() ) {
  exit 1, "The specified VEP config does not exist: ${params.vep_config}"
}

log.info 'Starting workflow.....'

workflow {
  chr = Channel.of(params.chros.split(','))
  splitVCF(chr, params.vcf)
  chrosVEP(splitVCF.out, params.vep_config)
  mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
}  
