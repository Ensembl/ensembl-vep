#!/bin/bash

echo 'Getting BioPerl'
if [ ! -f bioperl-release-1-2-3.zip ]; then
  wget https://github.com/bioperl/bioperl-live/archive/bioperl-release-1-2-3.zip
  unzip -q bioperl-release-1-2-3.zip
fi

echo 'Getting HTSlib'
if [ ! -d htslib ]; then
  git clone --branch master --depth 1 https://github.com/samtools/htslib.git
fi

echo 'Getting Bio::DB::HTS'
if [ ! -d Bio-HTS ]; then
  git clone --branch master --depth 1 https://github.com/Ensembl/Bio-HTS.git
fi

echo 'Getting jksrc'
if [ ! -f jksrc.zip ]; then
  wget http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
  unzip -q jksrc.zip
fi
