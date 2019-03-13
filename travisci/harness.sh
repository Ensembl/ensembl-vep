#!/bin/bash

# export DEPS=$HOME/dependencies

export PERL5LIB=$DEPS/bioperl-live:$PWD/ensembl-test/modules:$PWD/ensembl/modules:$PWD/modules:$PWD/ensembl-io/modules:$PWD/ensembl-funcgen/modules:$PWD/ensembl-variation/modules:$DEPS/Bio-HTS/blib/lib:$DEPS/Bio-HTS/blib/arch:$PERL5LIB

# export HTSLIB_DIR=$DEPS/htslib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DEPS/htslib

export PATH=$PATH:$DEPS/htslib

echo "Running test suite"
echo "Using $PERL5LIB"
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,modules/Bio/EnsEMBL/VEP/Pipeline,+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl,+ignore,ensembl-io,+ignore,ensembl-funcgen,+ignore,ensembl-variation,+ignore,Bio-HTS' perl $PWD/ensembl-test/scripts/runtests.pl -verbose $PWD/t $SKIP_TESTS
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
