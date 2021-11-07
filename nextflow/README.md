## Nextflow VEP pipeline
The nextflow pipeline aims to run VEP faster. It contains the following steps:
 * Splitting VCF chromosome-wise
 * Running VEP on chromosome-wise VCFs in parallel
 * Merging VEP outputs into a single file


##### Table of contents
* [Installation and requirements](#install)
  * [Fetching singularity images](#singularity)
* [Usage](#usage)

---
<a name="install"></a>
### Installation and requirements
The nextflow pipeline requires the following dependencies:
  * **Nextflow** (tested on 21.04.3.5560)
  * **Singularity** (tested on 3.7)

<a name="singularity"></a>
#### Fetching singularity images
By default all the images are expected under the `singularity-images` dir. 
  * bcftools
   ```bash
      singularity pull bcftools.sif docker://quay.io/biocontainers/bcftools:1.13--h3a49de5_0
   ```
  * [vep](https://www.ensembl.info/2021/05/24/cool-stuff-the-vep-can-do-singularity/)
  
In order to point to a different path to the image, the nextflow config file can be modified accordingly. 

---
<a name="usage"></a>
### Usage
```bash
  nextflow -C nf_config/run_vep.config run workflows/run_vep.nf --vcf /path/to/vcf --chros 1,2 --vep_config /path/to/vep.ini
```

#### Options
```bash
  --vcf VCF                VCF that will be split.
  --outdir DIRNAME         Name of output dir. Default: outdir
  --vep_config FILENAME    VEP config file. Default: nf_config/vep.ini
  --chros LIST_OF_CHROS    Comma-separated list of chromosomes to generate. i.e. chr1,chr2,...
  --cpus INT               Number of CPUs to use. Default 1.
```

