/* 
 * Workflow to run VEP on multiple input files
 *
 * Requires Nextflow: https://nextflow.io
 */

nextflow.enable.dsl=2

// module imports
include { checkVCF } from '../modules/check_VCF.nf'
include { generateSplits } from '../modules/generate_splits.nf'
include { splitVCF } from '../modules/split_VCF.nf' 
include { mergeVCF } from '../modules/merge_VCF.nf'  
include { runVEP as runVEPonVCF } from '../modules/run_vep.nf'
include { runVEP } from '../modules/run_vep.nf'

// print usage
if (params.help) {
  log.info """
Pipeline to run VEP
-------------------

Usage:
  nextflow run workflows/run_vep.nf --input <path-to-file> --vep_config vep_config/vep.ini

Options:
  --input FILE              Input file (if unsorted, use --sort to avoid errors). Alternatively, can also be a directory containing input files
  --bin_size INT            Number of lines to split input into multiple jobs. Default: 100
  --vep_config FILENAME     VEP config file. Alternatively, can also be a directory containing VEP INI files. Default: vep_config/vep.ini
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output directory. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <vcf>-<output_prefix>.vcf.gz

  --sort                    Sort VCF results from VEP (only required if input is unsorted; slower if enabled). Default: false
  --filter STRING           Comma-separated list of filter conditions to pass to filter_vep, such as "AF < 0.01,Feature is ENST00000377918".
                            Read more on how to write filters at https://ensembl.org/info/docs/tools/vep/script/vep_filter.html
                            Default: null (filter_vep is not run)
  """
  exit 1
}


def createInputChannels (input, pattern) {
  files = file(input)
  if ( !files.exists() ) {
    exit 1, "The specified input does not exist: ${input}"
  }

  if (files.isDirectory()) {
    files = "${files}/${pattern}"
  }
  files = Channel.fromPath(files)
  
  return files;
}

def createOutputChannel (output) {
  def dir = new File(output)

  // convert output dir to absolute path if necessary
  if (!dir.isAbsolute()) {
      output = "${launchDir}/${output}";
  }

  return Channel.fromPath(output)
}

workflow vep {
  take:
    inputs
  main:
    // Process input based on file extension
    inputs |
      branch {
        index: it.file =~ '\\.(tbi|csi)$'
        registry: it.file =~ '\\.registry$'
        ini: it.file =~ '\\.ini$'
        vcf: it.file =~ '\\.vcf(.gz)?$'
        other: true
      } |
      set { data }

    // Run VEP on VCF files
    data.vcf |
      checkVCF |
      // Generate split files that each contain bin_size number of variants
      generateSplits | transpose |
      // Split VCF using split files
      splitVCF | transpose |
      // Run VEP for each split VCF file and for each VEP config
      map { it + [format: 'vcf'] } | runVEPonVCF

    // Run VEP on non-VCF files
    data.other |
      map {
          // Split input by bin_size
          files = it.file.splitText(by: params.bin_size, file: true)
          res = []
          for (f : files) {
            // put it.file as index to avoid Nextflow errors
            res += [ meta: it.meta, original: it.file, file: f, index: it.file, vep_config: it.vep_config, format: 'other' ]
          }
          res
      } |
      flatten |
      runVEP

    // Merge split VCF files (creates one output VCF for each input VCF)
    out = runVEP.out.files
            .mix(runVEPonVCF.out.files)
            .groupTuple(by: [0, 1, 4])
    mergeVCF(out)
  emit:
    mergeVCF.out
}

workflow NF_VEP {
  if (!params.input) {
    exit 1, "Undefined --input parameter. Please provide the path to an input file."
  }

  if (params.vcf) {
    log.warn "The --vcf parameter is deprecated in Nextflow VEP. Please use --input instead."
  }

  if (!params.vep_config) {
    exit 1, "Undefined --vep_config parameter. Please provide a VEP config file."
  }

  input = createInputChannels(params.input, pattern="*")
  vep_config = createInputChannels(params.vep_config, pattern="*.ini")

  input.count()
    .combine( vep_config.count() )
    .subscribe{ if ( it[0] != 1 && it[1] != 1 ) 
      exit 1, "Detected many-to-many scenario between VCF and VEP config files - currently not supported" 
    }
    
  // set if it is a one-to-many situation (single VCF and multiple ini file)
  // in this situation we produce output files with different names
  one_to_many = input.count()
    .combine( vep_config.count() )
    .map{ it[0] == 1 && it[1] != 1 }

  output_dir = createOutputChannel(params.outdir)
  
  filters = Channel.of(params.filters)
  
  input
    .combine( vep_config )
    .combine( one_to_many )
    .combine( output_dir )
    .combine( filters )
    .map {
      data, vep_config, one_to_many, output_dir, filters ->
        meta = [:]
        meta.one_to_many = one_to_many
        meta.output_dir = output_dir
        meta.filters = filters
        
        // NOTE: csi is default unless a tbi index already exists
        meta.index_type = file(data + ".tbi").exists() ? "tbi" : "csi"

        index = data + ".${meta.index_type}"

        [ meta: meta, file: data, index: index, vep_config: vep_config ]
    }
    .set{ ch_input }
  
  vep(ch_input)
}

workflow {
  NF_VEP()
}
