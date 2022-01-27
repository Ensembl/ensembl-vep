## Nextflow VEP pipeline

The nextflow pipeline aims to run VEP faster utilising simple parallelisation. It is deployable on an individual Linux machine or on computing clusters running lsf or slurm (not tested). The process can be summarised briefly by the following steps:

 * Splitting the VCF chromosome-wise
 * Running VEP on chromosome-wise VCFs in parallel
 * Merging VEP outputs into a single file


##### Table of contents

* [Installation and requirements](#install)
* [Pipeline setup](#setup)
* [Usage](#usage)
* [Example](#example)

---
<a id="install"></a>
### Installation and requirements

The nextflow pipeline requires the following dependencies:

  * **[Nextflow](https://www.nextflow.io)** (tested on 21.10.0)
  * **[Singularity](https://sylabs.io/guides/3.7/admin-guide/installation.html)** (tested on 3.7)

<a id="setup"></a>
### Pipeline setup

#### Singularity images
Singularity images are required in order to run the following tools:

  * bcftools
  * VEP

The singularity images can be fetched by running:
```bash
   ./setup-images.sh
```

#### Config files

The following config files are used and can be modified depending on user requirements:

  * VEP config file
  ``` bash
      cp nf_config/vep.ini.template nf_config/vep.ini
  ```


  * Nextflow config file

    `nf_config/nextflow.config` has the default options for running the pipeline. The file can be modified to change the default options or override them using command line options

 Currently supported profiles for executors are standard (local), LSF and SLURM (untested!). As mentioned SLURM is untested at present, if you are running this pipeline on a slurm compute cluster and encounter problems, please contact us with details (raise a ticket on the github) and we can investigate.
 NB: If no profile is mentioned, the pipeline takes the standard profile.

---
<a id="usage"></a>
### Usage

```bash
  nextflow run workflows/run_vep.nf \
  -C nf_config/nextflow.config \
  --vcf <path-to-vcf> \
  --chros 1,2 \
  -profile <standard or lsf or slurm>
```

#### Options

```bash
  --vcf VCF               VCF that will be split. Currently supports sorted and bgzipped file
  --outdir DIRNAME        Name of output dir. Default: outdir
  --vep_config FILENAME   VEP config file. Default: nf_config/vep.ini
  --chros LIST_OF_CHROS   Comma-separated list of chromosomes to generate. i.e. 1,2,..., Default: 1,2,...X,Y,MT
  --cpus INT              Number of CPUs to use. Default 1.
```
NB: File paths are expected to be absolute paths.

---
<a id="example"></a>
### Example

```bash
  bgzip -c $PWD/examples/clinvar-testset/input.vcf > $PWD/examples/clinvar-testset/input.vcf.gz

  nextflow -C nf_config/nextflow.config \
    run workflows/run_vep.nf \
    --vcf $PWD/examples/clinvar-testset/input.vcf.gz \
    -profile lsf
 ```
The above commands start the pipeline and generate the output file upon completion.

#### Output validation

```bash
  singularity-images/bcftools.sif bcftools view \
    -H outdir/merged-file.vcf.gz \
    -r 1
```
Expected result

```bash
1 925952  1019397 G A . . ALLELEID=1003021;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.11:g.925952G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;CSQ=A|upstream_gene_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000341065|protein_coding|||||||||||4360|1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000342066|protein_coding|2/14||||101|11|4|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000437963|protein_coding|2/5||||71|11|4|G/E|gGg/gAg|||1|cds_end_NF|HGNC|HGNC:28706,A|upstream_gene_variant|MODIFIER|LINC02593|ENSG00000223764|Transcript|ENST00000609207|retained_intron|||||||||||4936|-1||HGNC|HGNC:53933,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616016|protein_coding|2/14||||1057|548|183|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616125|protein_coding|1/11||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000617307|protein_coding|1/13||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618181|protein_coding|1/10||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618323|protein_coding|2/14||||1057|548|183|G/E|gGg/gAg|||1||HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618779|protein_coding|1/12||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706,A|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000622503|protein_coding|1/13||||11|11|4|G/E|gGg/gAg|||1|cds_start_NF|HGNC|HGNC:28706
1 930139  1125147 C T . . ALLELEID=1110865;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.11:g.930139C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001627|intron_variant;ORIGIN=1;CSQ=T|upstream_gene_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000341065|protein_coding|||||||||||173|1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000342066|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000437963|protein_coding||2/4||||||||||1|cds_end_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000616016|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000616125|protein_coding||1/10||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000617307|protein_coding||1/12||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618181|protein_coding||1/9||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618323|protein_coding||2/13||||||||||1||HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000618779|protein_coding||1/11||||||||||1|cds_start_NF|HGNC|HGNC:28706,T|intron_variant|MODIFIER|SAMD11|ENSG00000187634|Transcript|ENST00000622503|protein_coding||1/12||||||||||1|cds_start_NF|HGNC|HGNC:28706
```
