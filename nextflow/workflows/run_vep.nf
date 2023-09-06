/* 
 * Workflow to run VEP on VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

nextflow.enable.dsl=2

 // params default
params.cpus = 1

params.vcf = null
params.vep_config = null
params.outdir = "outdir"

params.output_prefix = ""
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


def createInputChannels (input, pattern) {
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

def createOutputChannel (output) {
  if (output instanceof String){
    def dir = new File(output)
  
    // convert output dir to absolute path if necessary
    if (!dir.isAbsolute()) {
        output = "${launchDir}/${output}";
    }
    return Channel.fromPath(output)
  }
  
  return output;
}

workflow vep {
  take:
    vcf
    vep_config
    output_dir
  main:
    if (!vcf) {
      exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file."
    }

    if (!vep_config) {
      exit 1, "Undefined --vep_config parameter. Please provide a VEP config file."
    }
    
    vcf = createInputChannels(vcf, pattern="*.{vcf,gz}")
    vep_config = createInputChannels(vep_config, pattern="*.ini")

    vcf.count()
      .combine( vep_config.count() )
      .subscribe{ if ( it[0] != 1 && it[1] != 1 ) 
        exit 1, "Detected many-to-many scenario between VCF and VEP config files - currently not supported" 
      }
      
    output_dir = createOutputChannel(output_dir)
    
    // set if it is a one_to_many situation (single VCF and multiple ini file)
    // in this situation we produce output files with different names
    one_to_many = vcf.count()
      .combine( vep_config.count() )
      .map{ it[0] == 1 && it[1] != 1 }

    // process input and create Channel
    // this works like 'merge' operator and thus might make the pipeline un-resumable
    // we might think of using 'toSortedList' and generate appropriate input from the 'processInput' module
    processInput(vcf.combine(vep_config))
    
    // Prepare input VCF files (bgzip + tabix)
    checkVCF(processInput.out)
    
    // Generate split files that each contain bin_size number of variants from VCF
    generateSplits(checkVCF.out, params.bin_size)

    // Split VCF using split files
    splitVCF(generateSplits.out.transpose())

    // Run VEP for each split VCF file and for each VEP config
    runVEP(splitVCF.out.transpose())
    
    // Merge split VCF files (creates one output VCF for each input VCF)
    mergeVCF(runVEP.out.files.groupTuple(by: [0, 3, 4])
      .combine(one_to_many
      .combine(output_dir))
    )
  emit:
    mergeVCF.out
}

workflow {
  vep(params.vcf, params.vep_config, params.outdir)
}
