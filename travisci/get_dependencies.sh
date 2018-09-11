#!/bin/bash

echo 'Getting BioPerl'
if [ ! -f release-1-6-924.zip ]; then
  wget https://github.com/bioperl/bioperl-live/archive/release-1-6-924.zip
  unzip -q release-1-6-924.zip
fi

echo 'Getting HTSlib'
if [ ! -d htslib ]; then
  wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
  tar xvjf htslib-1.9.tar.bz2
  mv htslib-1.9 htslib
fi

echo 'Getting Bio::DB::HTS'
if [ ! -d Bio-HTS ]; then
  git clone --branch release/v2.11 --depth 1 https://github.com/Ensembl/Bio-HTS.git
fi

echo 'Getting jksrc'
if [ ! -f v335_base.tar.gz ]; then
  wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
  tar xzf v335_base.tar.gz
fi
