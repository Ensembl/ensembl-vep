[![Coverage Status](https://coveralls.io/repos/github/willmclaren/ensembl-vep/badge.svg?branch=master)](https://coveralls.io/github/willmclaren/ensembl-vep?branch=master)
# ensembl-vep
* The Variant Effect Predictor (VEP) predicts the functional effects of genomic variants.
* [Haplosaurus][haplo] uses phased genotype data to predict whole-transcript haplotype sequences.

### Installation and requirements
The VEP package requires Perl (>=5.10 recommended), the Ensembl API and a few other Perl modules. The Ensembl API in turn requires the [BioPerl](https://github.com/bioperl/bioperl-live), [DBI](http://search.cpan.org/~timb/DBI/DBI.pm) and [DBD::mysql](http://search.cpan.org/~michielb/DBD-mysql-4.036/lib/DBD/mysql.pm) packages to be installed.
* Ensembl API modules ([installation guide](http://www.ensembl.org/info/docs/api/api_installation.html), [git install](http://www.ensembl.org/info/docs/api/api_git.html))
  * ensembl-variation
  * ensembl
  * ensembl-funcgen
  * ensembl-io
* CPAN modules - we recommend using [cpanminus](http://search.cpan.org/~miyagawa/Menlo-1.9003/script/cpanm-menlo) to install
  * [Set::IntervalTree](http://search.cpan.org/~benbooth/Set-IntervalTree/lib/Set/IntervalTree.pm)
  * [Bio::DB::HTS](http://search.cpan.org/dist/Bio-DB-HTS/)
  * [JSON](http://search.cpan.org/dist/JSON/)

Additional requirements for non-core functionality
* [PerlIO::gzip](http://search.cpan.org/~nwclark/PerlIO-gzip-0.19/gzip.pm) - faster compressed file parsing
* [Bio::DB::BigFile](http://search.cpan.org/~lds/Bio-BigFile-1.07/lib/Bio/DB/BigFile.pm) - required for reading custom annotation data from BigWig files
 
### Usage
```bash
perl vep.pl -i input.vcf -o out.txt -cache
```
See [documentation](http://www.ensembl.org/info/docs/tools/vep/script/index.html) for full command line instructions.

> Note that this documentation currently corresponds to the ensembl-tools version of VEP; this will be updated to reflect the new version soon. Almost all commands and flags should work on both versions.

---
## Haplosaurus
[haplo]:yes
