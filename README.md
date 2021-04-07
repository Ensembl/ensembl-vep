# ensembl-vep

[![GitHub](https://img.shields.io/github/license/Ensembl/ensembl-vep.svg)](https://github.com/Ensembl/ensembl-vep/blob/release/104/LICENSE)
[![Coverage Status](https://coveralls.io/repos/github/Ensembl/ensembl-vep/badge.svg?branch=release/104)](https://coveralls.io/github/Ensembl/ensembl-vep?branch=release/104)
[![Docker Build Status](https://img.shields.io/docker/build/ensemblorg/ensembl-vep.svg)](https://hub.docker.com/r/ensemblorg/ensembl-vep)
[![Docker Hub Pulls](https://img.shields.io/docker/pulls/ensemblorg/ensembl-vep.svg)](https://hub.docker.com/r/ensemblorg/ensembl-vep)

* **VEP** (Variant Effect Predictor) predicts the functional effects of genomic variants.
* **Haplosaurus** uses phased genotype data to predict whole-transcript haplotype sequences.
* **Variant Recoder** translates between different variant encodings.

##### Table of contents
* [Installation and requirements](#install)
* [VEP](#vep)
  * [Usage](#vepusage)
* [Haplosaurus](#haplo)
  * [Usage](#haplousage)
  * [Output](#haplooutput)
  * [REST](#haploREST)
  * [Flags](#haploflags)
* [Variant Recoder](#recoder)
  * [Usage](#recoderusage)
  * [Output](#recoderoutput) 

---
<a name="install"></a>

### Installation and requirements
The VEP package requires:

  * **gcc**, **g++** and **make**
  * **Perl** (>=5.10 recommended, tested on 5.10, 5.14, 5.18, 5.22, 5.26)
  * Perl libraries [Archive::Zip](https://metacpan.org/pod/Archive::Zip) and [DBI](https://metacpan.org/pod/DBI)

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
 
#### Docker
A docker image for VEP is available from [DockerHub](https://hub.docker.com/r/ensemblorg/ensembl-vep).

See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#docker) for the Docker installation instructions.
 
---

<a name="vep"></a>
## VEP

<a name="vepusage"></a>

### Usage
```bash
./vep -i input.vcf -o out.txt -offline
```
See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic) for full command line instructions.

> Please report any bugs or issues by [contacting Ensembl](http://www.ensembl.org/info/about/contact/index.html) or creating a [GitHub issue](https://github.com/Ensembl/ensembl-vep/issues)

---
<a name="haplo"></a>
## Haplosaurus
`haplo` is a local tool implementation of the same functionality that powers the [Ensembl transcript haplotypes view](http://www.ensembl.org/Homo_sapiens/Transcript/Haplotypes?t=ENST00000304748). It takes phased genotypes from a VCF and constructs a pair of haplotype sequences for each overlapped transcript; these sequences are also translated into predicted protein haplotype sequences. Each variant haplotype sequence is aligned and compared to the reference, and an HGVS-like name is constructed representing its differences to the reference.

This approach offers an advantage over VEP's analysis, which treats each input variant independently. By considering the combined change contributed by all the variant alleles across a transcript, the compound effects the variants may have are correctly accounted for.

`haplo` shares much of the same command line functionality with `vep`, and can use VEP caches, Ensembl databases, GFF and GTF files as sources of transcript data; all `vep` command line flags relating to this functionality work the same with `haplo`.

<a name="haplousage"></a>
### Usage
Input data must be a [VCF](http://samtools.github.io/hts-specs/VCFv4.3.pdf) containing phased genotype data for at least one individual and file must be sorted by chromosome and genomic position; no other formats are currently supported.

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
6. Comma-separated list of frequency data for protein haplotype
7. Comma-separated list of contributing variants
8. Comma-separated list of sample:count that exhibit this haplotype

The altered haplotype sequences can be obtained by switching to JSON output using `--json` which will display them by default.
Each transcript analysed is summarised as a JSON object written to one line of the output file.

The [JSON output](#haploREST) structure matches the format of the [transcript haplotype REST endpoint](https://rest.ensembl.org/documentation/info/transcript_haplotypes_get).

You may exclude fields in the JSON from being exported with `--dont_export field1,field2`. This may be used, for example, to exclude the full haplotype sequence and aligned sequences from the output with `--dont_export seq,aligned_sequences`.

> Note JSON output does not currently include side-loaded frequency data.



<a name="haploREST"></a>
### REST service
The [transcript haplotype REST endpoint](https://rest.ensembl.org/documentation/info/transcript_haplotypes_get).
returns arrays of protein_haplotypes and cds_haplotypes for a given transcript. The default haplotype record includes:

* **population_counts**: the number of times the haplotype is seen in each population
* **population_frequencies**: the frequency of the haplotype  in each population
* **contributing_variants**:  variants contributing to the haplotype
* **diffs**: differences between the reference and this haplotype
* **hex**: the md5 hex of this haplotype sequence
* **other_hexes**: the md5 hex of other related haplotype sequences (
        CDSHaplotypes that translate to this ProteinHaplotype or
        ProteinHaplotype representing the translation of this CDSHaplotype)
* **has_indel**: does the haplotype contain insertions or deletions
* **type**: the type of haplotype - cds, protein
* **name**: a human readable name for the haplotype (sequence id + REF or a change description)
* **flags**: [flags](#haploflags) for the haplotype
* **frequency**: haplotype frequency in full sample set
* **count**: haplotype count in full sample set

The REST service does not return raw sequences, sample-haplotype assignments and the aligned sequences used to generate
differences by default.



<a name="haploflags"></a>
### Flags
Haplotypes may be flagged with one or more of the following:
* **indel**: haplotype contains an insertion or deletion (indel) relative to the reference.
* **frameshift:** haplotype contains at least one indel that disrupts the reading frame of the transcript.
* **resolved_frameshift:** haplotype contains two or more indels whose combined effect restores the reading frame of the transcript.
* **stop_changed:** indicates either a STOP codon is gained (protein truncating variant, PTV) or the existing reference STOP codon is lost.
* **deleterious_sift_or_polyphen:** haplotype contains at least one single amino acid substitution event flagged as deleterious (SIFT) or probably damaging (PolyPhen2).

<a name="bioperl-ext"></a>
### bioperl-ext
`haplo` can make use of a fast compiled alignment algorithm from the [bioperl-ext](https://github.com/bioperl/bioperl-ext) package; this can speed up analysis, particularly in longer transcripts where insertions and/or deletions are introduced. The bioperl-ext package is no longer maintained and requires some tweaking to install. The following instructions install the package in `$HOME/perl5`; edit `PREFIX=[path]` to change this. You may also need to edit the `export` command to point to the path created for the architecture on your machine.

```bash
git clone https://github.com/bioperl/bioperl-ext.git
cd bioperl-ext/Bio/Ext/Align/
perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
perl Makefile.PL PREFIX=~/perl5
make
make install
cd -
export PERL5LIB=${PERL5LIB}:${HOME}/perl5/lib/x86_64-linux-gnu/perl/5.22.1/
```

If successful the following should print `OK`:

```bash
perl -MBio::Tools::dpAlign -e"print qq{OK\n}"
```

---
<a name="recoder"></a>
## Variant Recoder
`variant_recoder` is a tool for translating between different variant encodings. It accepts as input any format supported by VEP (VCF, variant ID, HGVS), with extensions to allow for parsing of potentially ambiguous HGVS notations. For each input variant, `variant_recoder` reports all possible encodings including variant IDs from [all sources imported into the Ensembl database](http://www.ensembl.org/info/genome/variation/species/sources_documentation.html) and HGVS (genomic, transcript and protein), reported on Ensembl, RefSeq and LRG sequences.

<a name="recoderusage"></a>
### Usage
`variant_recoder` depends on database access for identifier lookup, and cannot be used in offline mode as per VEP. The output format is JSON and the [JSON perl module](http://search.cpan.org/dist/JSON/) is required.

```bash
./variant_recoder --id [input_data_string]
./variant_recoder -i [input_file] --species [species]
```

<a name="recoderoutput"></a>
### Output
Output is a JSON array of objects, one per input variant, with the following keys:
* **input**: input string
* **id**: variant identifiers
* **hgvsg**: HGVS genomic nomenclature
* **hgvsc**: HGVS transcript nomenclature
* **hgvsp**: HGVS protein nomenclature
* **spdi**: Genomic SPDI notation
* **vcf_string**: VCF format (optional)
* **warnings**: Warnings generated e.g. for invalid HGVS

Use `--pretty` to pre-format and indent JSON output.

Example output:
```bash
./variant_recoder --id "AGT:p.Met259Thr" --pretty
[
   {
     "warnings" : [
         "Possible invalid use of gene or protein identifier 'AGT' as HGVS reference; AGT:p.Met259Thr may resolve to multiple genomic locations"
      ],
     "C" : {
        "input" : "AGT:p.Met259Thr",
        "id" : [
           "rs699",
           "CM920010",
           "COSV64184214"
        ],
        "hgvsg" : [
           "NC_000001.11:g.230710048A>G"
        ],
        "hgvsc" : [
           "ENST00000366667.6:c.776T>C",
           "ENST00000679684.1:c.776T>C",
           "ENST00000679738.1:c.776T>C",
           "ENST00000679802.1:c.776T>C",
           "ENST00000679854.1:n.1287T>C",
           "ENST00000679957.1:c.776T>C",
           "ENST00000680041.1:c.776T>C",
           "ENST00000680783.1:c.776T>C",
           "ENST00000681269.1:c.776T>C",
           "ENST00000681347.1:n.1287T>C",
           "ENST00000681514.1:c.776T>C",
           "ENST00000681772.1:c.776T>C",
           "NM_001382817.3:c.776T>C",
           "NM_001384479.1:c.776T>C"
        ],
        "hgvsp" : [
           "ENSP00000355627.5:p.Met259Thr",
           "ENSP00000505981.1:p.Met259Thr",
           "ENSP00000505063.1:p.Met259Thr",
           "ENSP00000505184.1:p.Met259Thr",
           "ENSP00000506646.1:p.Met259Thr",
           "ENSP00000504866.1:p.Met259Thr",
           "ENSP00000506329.1:p.Met259Thr",
           "ENSP00000505985.1:p.Met259Thr",
           "ENSP00000505963.1:p.Met259Thr",
           "ENSP00000505829.1:p.Met259Thr",
           "NP_001369746.2:p.Met259Thr",
           "NP_001371408.1:p.Met259Thr"
        ],
        "spdi" : [
           "NC_000001.11:230710047:A:G"
        ]
     }
   }
]
```

<a name="recoderopts"></a>
### Options
`variant_recoder` shares many of the same command line flags as VEP. Others are unique to `variant_recoder`.

* `-id|--input_data [input_string]`: a single variant as a string.
* `-i|--input_file [input_file]`: input file containing one or more variants, one per line. Mixed formats disallowed.
* `--species`: species to use (default: homo_sapiens).
* `--grch37`: use GRCh37 assembly instead of GRCh38.
* `--genomes`: set database parameters for [Ensembl Genomes](http://ensemblgenomes.org/) species.
* `--pretty`: write pre-formatted indented JSON.
* `--fields [field1,field2]`: limit output fields. Comma-separated list, one or more of: `id`, `hgvsg`, `hgvsc`, `hgvsp`, `spdi`, `vcf_string`.
* `--host [db_host]`: change database host from default `ensembldb.ensembl.org` (UK); geographic mirrors are `useastdb.ensembl.org` (US East Coast) and `asiadb.ensembl.org` (Asia). `--user`, `--port` and `--pass` may also be set.
* `--pick`, `--per_gene`, `--pick_allele`, `--pick_allele_gene`, `--pick_order`: set and customise transcript selection process, see [VEP documentation](http://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
