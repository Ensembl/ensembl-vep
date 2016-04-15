#!/bin/bash

export PERL5LIB=$PWD/bioperl-live:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/modules:$PWD/ensembl-io/modules:$PWD/ensembl-funcgen/modules:$PWD/ensembl-variation/modules

export PATH=$PATH:$PWD/htslib

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,ensembl-io,+ignore,ensembl-funcgen,+ignore,ensembl-variation' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/t $SKIP_TESTS
else
  perl $PWD/ensembl-test/scripts/runtests.pl $PWD/t $SKIP_TESTS
fi

rt=$?
if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit $rt
fi
