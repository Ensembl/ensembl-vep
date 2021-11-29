## Nextflow VEP pipeline
The nextflow pipeline aims to run VEP faster. It contains the following steps:
 * Splitting VCF chromosome-wise
 * Running VEP on chromosome-wise VCFs in parallel
 * Merging VEP outputs into a single file


##### Table of contents
* [Installation and requirements](#install)
* [Pipeline setup](#setup)
* [Usage](#usage)

---
<a name="install"></a>
### Installation and requirements
The nextflow pipeline requires the following dependencies:
  * **Nextflow** (tested on 21.04.3.5560)
  * **Singularity** (tested on 3.7)

<a name="setup"></a>
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
The following config files are used:
  * VEP config file 
  ```bash
     cp nf_config/vep.ini.template nf_config/vep.ini
  ```
  Modify the vep config file as required

  * Nextflow config file
  `nf_config/nextflow.config` has the default options for running the pipeline. The file can be modified to change the default options or override them using command line options
  The default executor is LSF. This can be modified by following the nextflow [documentation](https://www.nextflow.io/docs/latest/executor.html)

---
<a name="usage"></a>
### Usage
```bash
  nextflow -C nf_config/nextflow.config run workflows/run_vep.nf --vcf <path-to-vcf> --chros 1,2 
```

#### Options
```bash
  --vcf VCF                VCF that will be split.
  --outdir DIRNAME         Name of output dir. Default: outdir
  --vep_config FILENAME    VEP config file. Default: nf_config/vep.ini
  --chros LIST_OF_CHROS    Comma-separated list of chromosomes to generate. i.e. 1,2,..., Default: All chromosomes
  --cpus INT               Number of CPUs to use. Default 1.
```

 
