> Used VEP? Please help us guide the future development of VEP by participating in our [user survey](https://docs.google.com/forms/d/e/1FAIpQLSe3PKvSviifr-zFjAc6xxQegOio17xTx_UlkDe4TMUoV0gbDg/viewform).

[![Coverage Status](https://coveralls.io/repos/github/Ensembl/ensembl-vep/badge.svg?branch=master)](https://coveralls.io/github/Ensembl/ensembl-vep?branch=master)
# ensembl-vep
* **VEP** (Variant Effect Predictor) predicts the functional effects of genomic variants.
* **Haplosaurus** uses phased genotype data to predict whole-transcript haplotype sequences.

##### Table of contents
* [Installation and requirements](#install)
* [VEP](#vep)
  * [Usage](#vepusage)
  * [Differences to ensembl-tools version](#vepdiffs)
* [Haplosaurus](#haplo)
  * [Usage](#haplousage)
  * [Output](#haplooutput)
  * [Flags](#haploflags)
  * [Frequency data](#haplofreq)

---
<a name="install"></a>

### Installation and requirements
The VEP package requires Perl (>=5.10 recommended, tested on 5.8, 5.10, 5.14, 5.18, 5.22) and the [DBI](http://search.cpan.org/~timb/DBI/DBI.pm) package installed.
The remaining dependencies can be installed using the included [INSTALL.pl](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer) script. Basic instructions:
```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```
The installer may also be used to check for updates to this and co-dependent packages, simply re-run INSTALL.pl.

See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) for full installation instructions.

#### Additional CPAN modules
The following modules are optional but most users will benefit from installing them. We recommend using [cpanminus](http://search.cpan.org/~miyagawa/Menlo-1.9003/script/cpanm-menlo) to install.
  * [DBD::mysql](http://search.cpan.org/~michielb/DBD-mysql/lib/DBD/mysql.pm) - required for database access (`--database` or `--cache` without `--offline`)
  * [Set::IntervalTree](http://search.cpan.org/~benbooth/Set-IntervalTree/lib/Set/IntervalTree.pm) - required for Haplosaurus, also confers speed updates to VEP
  * [JSON](http://search.cpan.org/dist/JSON/) - required for writing JSON output
  * [PerlIO::gzip](http://search.cpan.org/~nwclark/PerlIO-gzip-0.19/gzip.pm) - faster compressed file parsing
  * [Bio::DB::BigFile](http://search.cpan.org/~lds/Bio-BigFile-1.07/lib/Bio/DB/BigFile.pm) - required for reading custom annotation data from BigWig files
 
---

<a name="vep"></a>
## VEP

<a name="vepusage"></a>

### Usage
```bash
./vep -i input.vcf -o out.txt -offline
```
See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/index.html) for full command line instructions.

> Please report any bugs or issues by [contacting Ensembl](http://www.ensembl.org/info/about/contact/index.html) or creating a [GitHub issue](https://github.com/Ensembl/ensembl-vep/issues)

<a name="vepdiffs"></a>

### Differences to ensembl-tools VEP
This ensembl-vep repo is a complete rewrite of the VEP code intended to make the software faster, more robust and more easily extensible. Almost all functionality of the ensembl-tools version has been replicated, with the command line flags remaining largely unchanged. A summary of changes follows:

* **Tool name:** For brevity and to distinguish the two versions, the new command line tool is named `vep`, with the version in ensembl-tools named `variant_effect_predictor.pl`.
* **Speed:** A typical individual human genome of 4 million variants can now be processed in around 30 minutes on a quad-core machine using under 1GB of RAM.
* **Known/existing variants:** The alleles of your input variant are now compared to any known variants when using `--check_existing`. Previously this would require you to enable this functionality manually with `--check_alleles`. The old functionality can be restored using `--no_check_alleles`.
* **Allele frequencies:** Allele frequencies are now reported for the input allele only e.g. as `0.023` instead of `A:0.023,G:0.0005`. To reflect this change, the allele frequency fields are now named e.g. `AFR_AF` instead of `AFR_MAF`. The command line flags reflect this also, so `--gmaf` is now `--af` and `--maf_1kg` is now `--af_1kg`. Using the old flags will produce a deprecation message.
* **GFF and GTF files:** GFF and GTF files may now be used directly as a source of transcript annotation in place of, or even alongside, a cache or database source. Previously this involved [building a cache using gtf2vep](http://e87.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gtf), which is now redundant. The files must first be bgzipped and tabix-indexed, and a FASTA file containing genomic sequence is required:
```bash
grep -v "#" data.gff | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > data.gff.gz
tabix -p gff data.gff.gz
./vep -i input.vcf -gff data.gff.gz -fasta genome.fa.gz
```
* **VCF custom annotations:** [VCF files used as a source of custom annotation](http://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html) will now have allele-specific data added from INFO fields; previously the whole content of each requested KEY=VALUE pair was reported.
* **New pick flags:** New flags added to aid [selecting amongst consequence output](http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick): `--pick_allele_gene`, `--flag_pick_allele_gene`
* **Runtime status:** `vep` produces no runtime progress messages.
* **Deprecated:**
  * GVF output: `--gvf`
  * HTML output: `--html`
  * format conversion: `--convert`
  * pileup input: `--format pileup`
  * MAF flags (replaced by AF flags): `--gmaf` (`--af`), `--maf_1kg` (`--af_1kg`), `--maf_esp` (`--af_esp`), `--maf_exac` (`--af_exac`)
  * known variant allele checking (on by default, use `--no_check_alleles` to restore old behaviour):  `--check_alleles` 
  * cache building flags (replaced by internal Ensembl pipeline): `--build`, `--write_cache`

---
<a name="haplo"></a>
## Haplosaurus
`haplo` is a local tool implementation of the same functionality that powers the [Ensembl transcript haplotypes view](http://www.ensembl.org/Homo_sapiens/Transcript/Haplotypes?t=ENST00000304748). It takes phased genotypes from a VCF and constructs a pair of haplotype sequences for each overlapped transcript; these sequences are also translated into predicted protein haplotype sequences. Each variant haplotype sequence is aligned and compared to the reference, and an HGVS-like name is constructed representing its differences to the reference.

This approach offers an advantage over VEP's analysis, which treats each input variant independently. By considering the combined change contributed by all the variant alleles across a transcript, the compound effects the variants may have are correctly accounted for.

`haplo` shares much of the same command line functionality with `vep`, and can use VEP caches, Ensembl databases, GFF and GTF files as sources of transcript data; all `vep` command line flags relating to this functionality work the same with `haplo`.

<a name="haplousage"></a>
### Usage
Input data must be a [VCF](http://samtools.github.io/hts-specs/VCFv4.3.pdf) containing phased genotype data for at least one individual; no other formats are currently supported.

When using a VEP cache as the source of transcript annotation, the first time you run `haplo` with a particular cache it will spend some time scanning transcript locations in the cache.

```bash
./haplo -i input.vcf -o out.txt -cache
```

<a name="haplooutput"></a>
### Output
The default output format is a simple tab-delimited file reporting all observed non-reference haplotypes. It has the following fields:

1. Transcript stable ID
2. CDS haplotype name
3. Comma-separated list of [flags](#haploflags) for CDS haplotype
4. Protein haplotype name
5. Comma-separated list of [flags](#haploflags) for protein haplotype
6. Comma-separated list of [frequency data](#haplofreq) for protein haplotype
7. Comma-separated list of contributing variants
8. Comma-separated list of sample:count that exhibit this haplotype

Alternatively, JSON output matching the format of the [transcript haplotype REST endpoint](http://rest.ensembl.org/documentation/info/transcript_haplotypes_get) may be generated by using `--json`.

Each transcript analysed is summarised as a JSON object written to one line of the output file.

You may exclude fields in the JSON from being exported with `--dont_export field1,field2`. This may be used, for example, to exclude the full haplotype sequence and aligned sequences from the output with `--dont_export seq,aligned_sequences`.

> Note JSON output does not currently include side-loaded frequency data.

<a name="haploflags"></a>
### Flags
Haplotypes may be flagged with one or more of the following:
* **indel**: haplotype contains an insertion or deletion (indel) relative to the reference.
* **frameshift:** haplotype contains at least one indel that disrupts the reading frame of the transcript.
* **resolved_frameshift:** haplotype contains two or more indels whose combined effect restores the reading frame of the transcript.
* **stop_changed:** indicates either a STOP codon is gained (protein truncating variant, PTV) or the existing reference STOP codon is lost.
* **deleterious_sift_or_polyphen:** haplotype contains at least one single amino acid substitution event flagged as deleterious (SIFT) or probably damaging (PolyPhen2).

<a name="haplofreq"></a>
### Frequency data
Haplotype frequencies may be loaded and assigned to observed haplotypes using `--haplotype_frequencies [file]`. The following files may be used:

* 1000 genomes frequencies (GRCh37): [protein_haplotype_freqs_1KG_e85_GRCh37.txt.gz](https://dl.dropboxusercontent.com/u/12936195/protein_haplotype_freqs_1KG_e85_GRCh37.txt.gz)
* 1000 genomes frequencies (GRCh38): [protein_haplotype_freqs_1KG_e85_GRCh38.txt.gz](https://dl.dropboxusercontent.com/u/12936195/protein_haplotype_freqs_1KG_e85_GRCh38.txt.gz)

> Note these files are temporarily hosted on 3rd party servers and may be subject to change or removal while the software remains in the development phase.


