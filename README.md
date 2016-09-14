[![Coverage Status](https://coveralls.io/repos/github/willmclaren/ensembl-vep/badge.svg?branch=master)](https://coveralls.io/github/willmclaren/ensembl-vep?branch=master)
# ensembl-vep
* The Variant Effect Predictor (VEP) predicts the functional effects of genomic variants.
* [Haplosaurus](#haplo) uses phased genotype data to predict whole-transcript haplotype sequences.

> **!!! IMPORTANT !!!** This is pre-release code. Use at your own risk. Please continue to use the version of [VEP](http://www.ensembl.org/vep) in [ensembl-tools](https://github.com/Ensembl/ensembl-tools) if you are unsure.


### Installation and requirements
The VEP package requires Perl (>=5.10 recommended), the Ensembl API and a few other Perl modules. The Ensembl API in turn requires the [BioPerl](https://github.com/bioperl/bioperl-live), [DBI](http://search.cpan.org/~timb/DBI/DBI.pm) and [DBD::mysql](http://search.cpan.org/~michielb/DBD-mysql-4.036/lib/DBD/mysql.pm) packages to be installed.
* Ensembl API modules - use git to install ([Ensembl git install instructions](http://www.ensembl.org/info/docs/api/api_git.html))
  * ensembl-variation
  * ensembl
  * ensembl-funcgen
  * ensembl-io

> **IMPORTANT:** ensembl-variation and ensembl-io currently must be on the master or (when available) release/86 branch:

```bash
cd ensembl-variation
git checkout master
cd ../ensembl-io
git checkout master
cd ../
```

* CPAN modules - we recommend using [cpanminus](http://search.cpan.org/~miyagawa/Menlo-1.9003/script/cpanm-menlo) to install
  * [Set::IntervalTree](http://search.cpan.org/~benbooth/Set-IntervalTree/lib/Set/IntervalTree.pm)
  * [Bio::DB::HTS](http://search.cpan.org/dist/Bio-DB-HTS/) - requires compiled [htslib](https://github.com/samtools/htslib), set `$HTSLIB_DIR` to htslib path before installing Bio::DB::HTS
  * [JSON](http://search.cpan.org/dist/JSON/)

Additional requirements for non-core functionality
* [PerlIO::gzip](http://search.cpan.org/~nwclark/PerlIO-gzip-0.19/gzip.pm) - faster compressed file parsing
* [Bio::DB::BigFile](http://search.cpan.org/~lds/Bio-BigFile-1.07/lib/Bio/DB/BigFile.pm) - required for reading custom annotation data from BigWig files
 
---
## VEP

### Usage
```bash
perl vep.pl -i input.vcf -o out.txt -cache
```
vep.pl is compatible with the same downloadable caches as the ensembl-tools VEP. See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/index.html) for full command line instructions.

> Note that the documentation hosted on ensembl.org and linked to from here currently corresponds to the ensembl-tools version of VEP; this will be updated to reflect the new version soon. Almost all commands, flags and plugins should work on both versions.

### Differences to ensembl-tools VEP
This ensembl-vep repo is a complete rewrite of the VEP code intended to make the software faster, more robust and more easily extensible. Almost all functionality of the ensembl-tools version has been replicated, with the command line flags remaining largely unchanged. A summary of changes follows:

* **Script name:** For brevity and to distinguish the two versions, the new script is named vep.pl, with the version in ensembl-tools named variant_effect_predictor.pl.
* **Known/existing variants:** The alleles of your input variant are now compared to any known variants when using `--check_existing`. Previously this would require you to enable this functionality manually with `--check_alleles`. The old functionality can be restored using `--no_check_alleles`.
* **Allele frequencies:** Allele frequencies are now reported for the input allele only e.g. as `0.023` instead of `A:0.023,G:0.0005`. To reflect this change, the allele frequency fields are now named e.g. `AFR_AF` instead of `AFR_MAF`. The command line flags reflect this also, so `--maf` is now `--af` and `--maf_1kg` is now `--af_1kg`. Using the old flags will produce a deprecation message.
* **GFF and GTF files:** GFF and GTF files may now be used directly as a source of transcript annotation in place of, or even alongside, a cache or database source. Previously this involved [building a cache using gtf2vep.pl](http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gtf), which is now redundant. The files must first be bgzipped and tabix-indexed, and a FASTA file containing genomic sequence is required:
```bash
$ grep -v "#" data.gff | sort -k1,1 -k4,4n -k5,5n | bgzip -c > data.gff.gz
$ tabix -p gff data.gff.gz
$ perl vep.pl -i input.vcf -gff data.gff.gz -fasta genome.fa.gz
```
* **VCF custom annotations:** [VCF files used as a source of custom annotation](http://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html) will now have allele-specific data added from INFO fields; previously the whole content of each requested KEY=VALUE pair was reported.
* **New pick flags:** New flags added to aid [selecting amongst consequence output](http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick): `--pick_allele_gene`, `--flag_pick_allele_gene`
* **Runtime status:** vep.pl produces no runtime progress messages.
* **Deprecated:**
  * GVF output: `--gvf`
  * HTML output: `--html`
  * format conversion: `--convert`
  * pileup input: `--format pileup`
  * MAF flags (replaced by AF flags): `--gmaf` (`--af`), `--maf_1kg` (`--af`), `--maf_esp` (`--af_esp`), `--maf_exac` (`--af_exac`)

---
<a name="haplo"></a>
## Haplosaurus
haplo.pl is a local tool implementation of the same functionality that powers the [Ensembl transcript haplotypes view](http://www.ensembl.org/Homo_sapiens/Transcript/Haplotypes?t=ENST00000304748).

It shares much of the same command line functionality with vep.pl, and can use VEP caches, Ensembl databases, GFF and GTF files as sources of transcript data; all vep.pl command line flags relating to this functionality work the same with haplo.pl.

Input data must be a [VCF](http://samtools.github.io/hts-specs/VCFv4.3.pdf) containing phased genotype data for at least one individual; no other formats are currently supported.

When using a VEP cache as the source of transcript annotation, the first time you run haplo.pl with a particular cache it will spend some time scanning transcript locations in the cache.

```bash
perl haplo.pl -i input.vcf -o out.txt -cache
```

Output data is currently a simple tab-delimited file reporting all observed non-reference haplotypes. It has the following fields:

1. Transcript stable ID
2. CDS haplotype name
3. Comma-separated list of flags for CDS haplotype
4. Protein haplotype name
5. Comma-separated list of flags for protein haplotype
6. Comma-separated list of [frequency data](#haplofreq) for protein haplotype
7. Sample identifier
8. Number of copies of this haplotype observed in sample

<a name="haplofreq"></a>
### Frequency data
Haplotype frequencies may be loaded and assigned to observed haplotypes using `--haplotype_frequencies [file]`. The following files may be used:

* 1000 genomes frequencies (GRCh37): [protein_haplotype_freqs_1KG_e85_GRCh37.txt.gz](https://dl.dropboxusercontent.com/u/12936195/protein_haplotype_freqs_1KG_e85_GRCh37.txt.gz)
* 1000 genomes frequencies (GRCh38): [protein_haplotype_freqs_1KG_e85_GRCh38.txt.gz](https://dl.dropboxusercontent.com/u/12936195/protein_haplotype_freqs_1KG_e85_GRCh38.txt.gz)

> Note these files are temporarily hosted on 3rd party servers and may be subject to change or removal while the software remains in the development phase.


