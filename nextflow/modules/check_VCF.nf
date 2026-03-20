#!/usr/bin/env nextflow

/* 
 * Script to check if the files are bgzipped and bgzip if not
 */

def checkVCFheader (f) {
  // Check file extension
  if (!(f  =~ '\\.vcf$') && !(f =~ '\\.vcf\\.b?gz$')) {
    return false
  }

  // Check if file is compressed
  if (f =~ '\\.b?gz$') {
    InputStream fileStream = new FileInputStream(f.toString())
    InputStream gzip = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzip)
    BufferedReader data = new BufferedReader(decoder)
    lines = data.lines()
  } else {
    lines = f.readLines()
  }

  // Check file header
  boolean is_vcf_format = false
  boolean has_header = false
  for( line : lines ) {
    if (!line =~ '^#') {
      // stop inspecting file when reaching a line not starting with hash
      break
    } else if (line =~ '^##fileformat=') {
      is_vcf_format = true
    } else if (line =~ '^#CHROM') {
      has_header = true
    }
  }
  return is_vcf_format && has_header
}

process checkVCF {
  /*
  Function to check input VCF files

  Returns
  -------
  Tuple of VCF, VCF index, vep config file, a output dir, and the index type of VCF file
  */

  cache 'lenient'
  cpus params.cpus
  label 'vep'
  errorStrategy 'terminate'

  afterScript "rm *.vcf *.vcf.tbi *.vcf.csi"

  input:
  tuple val(meta), path(vcf), path(vcf_index), path(vep_config)
  
  output:
  tuple val(meta), val(vcf.simpleName), path("${vcf.simpleName}-checked.vcf.gz"), path("${vcf.simpleName}-checked.vcf.gz.${meta.index_type}"), path(vep_config)

  script:
  index_type = meta.index_type
  tabix_arg = index_type == 'tbi' ? '' : '-C'

  outfile_name = "${vcf.simpleName}-checked.vcf.gz"

  sorted_file = ""
  sort_cmd = ""
  if( params.sort ) {
    isGzipped = vcf.extension == 'gz'
    sorted_file = "${vcf.baseName}.sorted.vcf"
    sorted_file += isGzipped ? ".gz" : ""
    cat_cmd   = isGzipped ? "zcat ${vcf}" : "cat ${vcf}"
    sort_cmd += "(${cat_cmd} | head -1000 | grep '^#'; ${cat_cmd} | grep -v '^#' | sort -k1,1d -k2,2n)"  // Sort
    sort_cmd += isGzipped ? "| bgzip -c" : ""  // Conditionally compress
    sort_cmd += " > ${sorted_file}"  // Write to output file
  }
  else {
    sorted_file = vcf.toString()
  }

  """
  ${sort_cmd}
  [[ "${sorted_file}" == *.gz ]] && mv ${sorted_file} ${outfile_name} || bgzip -c ${sorted_file} > ${outfile_name}
  [ -f ${outfile_name}.${index_type} ] || tabix ${tabix_arg} -p vcf -f ${outfile_name}

  # quickly test tabix -- ensures both bgzip and tabix are okay
  chr=\$(tabix -l ${outfile_name} | head -n1)
  tabix ${outfile_name} \${chr}:1-10001
  """
}
