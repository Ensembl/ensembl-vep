## Nextflow VEP pipeline

The nextflow pipeline aims to run VEP faster utilising simple parallelisation. It is deployable on an individual Linux machine or on computing clusters running LSF or SLURM. The process can be summarised briefly by the following steps:

 * Splitting the input data into multiple files using a given number of bins (100 by default)
 * Running VEP on the split files in parallel
 * Merging VEP outputs into a single file

##### Table of contents

* [Installation and requirements](#install)
* [Pipeline setup](#setup)
* [Usage](#usage)
* [Example](#example)

---
<a id="install"></a>
### Installation and requirements

The nextflow pipeline requires **[Nextflow](https://www.nextflow.io)** (tested on 21.10.0).

<a id="setup"></a>
### Pipeline setup

The following config files are used and can be modified depending on user requirements:

* VEP config file
    ``` bash
    cp vep_config/vep.ini.template vep_config/vep.ini
    ```
* Nextflow config file
    `nextflow.config` has the default options for running the pipeline. The default options can be changed by directly modifying this file or by overriding them using command line options.
* VEP cache (optional, but faster than using `database`)
    VEP cache can be downloaded and extracted from the [Ensembl FTP][].
    The VEP config file needs to contain the absolute path to the cache directory.
    More detailed instructions are available in our [VEP cache documentation][].

[Ensembl FTP]: https://ftp.ensembl.org/pub/current_variation/indexed_vep_cache
[VEP cache documentation]: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

---
<a id="usage"></a>

### Usage

```bash
  nextflow run main.nf \
  -profile <profiles-to-use> \
  --input <path-to-file> \
  --vep_config <path-to-vep-config>
```

#### Options

```bash
[Mandatory]
  --input FILE              Input file (if unsorted, use --sort to avoid errors in indexing the output file). Alternatively, can also be a directory containing input files
  --vep_config FILENAME     VEP config file. Alternatively, can also be a directory containing VEP INI files. Default: 'vep_config/vep.ini'

[OPTIONAL]
  --bin_size INT            Number of lines to split input into multiple jobs. Default: 100
  --vep_version VERSION     VEP version to use from Docker Hub (such as 113.0); only required when using Docker or Singularity profile. Default: 'latest'
  --cpus INT                Number of CPUs to use. Default: 1
  --outdir DIRNAME          Name of output directory. Default: outdir
  --output_prefix PREFIX    Output filename prefix. The generated output file will have name <output_prefix>_VEP.vcf.gz.
                            NOTE: Do not use this parameter if you are expecting multiple output files.

  --sort                    Sort VCF results from VEP (only required if input is unsorted; slower if enabled). Default: false
  --filters STRING          Comma-separated list of filter conditions to pass to filter_vep, such as "AF < 0.01,Feature is ENST00000377918".
                            Read more on how to write filters at https://ensembl.org/info/docs/tools/vep/script/vep_filter.html
                            Default: null (filter_vep is not run)
```

---
<a id="example"></a>

### Example

You need to create a `vep.ini` file based on `vep_config/vep.ini.template`.
Always use absolute paths in `vep.ini` when indicating directories.

```bash
  nextflow run main.nf \
    --input ../examples/clinvar-testset/input.vcf \
    --vep_config ../vep_config/vep.ini
    -profile slurm,singularity
 ```
The above command starts the pipeline in SLURM using the Singularity image and
generates the output file upon completion.

#### Output validation

```bash
  work/singularity/quay.io-biocontainers-bcftools-*.img bcftools view \
    -H outdir/merged-file.vcf.gz \
    -r 1
```
Expected result:

```bash
1 925952  1019397 G A . . ALLELEID=1003021;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.11:g.925952G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;CSQ=A|upstream_gene_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000341065|protein_coding|||||||||||4360|1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000342066|protein_coding|2/14||||101|11|4|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000437963|protein_coding|2/5||||71|11|4|G/E|gGg/gAg|||1|cds_end_NF|HGNC|HGNC:28706,A|upstream_gene_variant|MODIFIER|LINC02593|ENSG00000223764|Transcript|ENST00000609207|retained_intron|||||||||||4936|-1||HGNC|HGNC:53933,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616016|protein_coding|2/14||||1057|548|183|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616125|protein_coding|1/11||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000617307|protein_coding|1/13||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618181|protein_coding|1/10||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618323|protein_coding|2/14||||1057|548|183|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618779|protein_coding|1/12||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000622503|protein_coding|1/13||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706
1 930139  1125147 C T . . ALLELEID=1110865;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.11:g.930139C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;CSQ=T|upstream_gene_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000341065|protein_coding|||||||||||173|1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000342066|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000437963|protein_coding||2/4||||||||||1|cds_end_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000616016|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000616125|protein_coding||1/10||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000617307|protein_coding||1/12||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618181|protein_coding||1/9||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618323|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618779|protein_coding||1/11||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000622503|protein_coding||1/12||||||||||1|cds_start_NF|HGNC|HGNC:28706
```
