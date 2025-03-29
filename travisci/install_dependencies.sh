#!/bin/bash

if [[ "$COVERALLS" = 'true' ]]
then
  sed -i '/Devel::Cover/d' ensembl/cpanfile
  sed -i '/Devel::Cover/d' cpanfile
fi

cpanm --quiet --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
cpanm --quiet --installdeps --with-recommends --notest .

if [[ "$COVERALLS" = 'true' ]]
then 
  cpanm --quiet -n Devel::Cover::Report::Coveralls
fi
