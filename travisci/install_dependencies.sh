#!/bin/bash

if [[ "$COVERALLS" = 'true' && ! "${TRAVIS_PERL_VERSION}" =~ '5.10' ]]
then
  sed -i '/Devel::Cover/d' ensembl/cpanfile
  sed -i '/Devel::Cover/d' cpanfile
fi

cpanm --quiet --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
cpanm --quiet --installdeps --with-recommends --notest .

if [[ "$COVERALLS" = 'true' && ! "${TRAVIS_PERL_VERSION}" =~ '5.10' ]]
then 
  cpanm --quiet -n Devel::Cover::Report::Coveralls
fi
