#!/bin/bash

if [[ "$COVERALLS" = 'true' && ! "${TRAVIS_PERL_VERSION}" =~ '5.10' ]]
then
  awk '!/Devel::Cover/' ensembl/cpanfile > tmp_cpanfile && mv tmp_cpanfile ensembl/cpanfile
  awk '!/Devel::Cover/' cpanfile > tmp_cpanfile && mv tmp_cpanfile cpanfile
fi

cpanm --quiet --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
cpanm --quiet --installdeps --with-recommends --notest .

if [[ "$COVERALLS" = 'true' && ! "${TRAVIS_PERL_VERSION}" =~ '5.10' ]]
then 
  cpanm --quiet -n Devel::Cover::Report::Coveralls
fi