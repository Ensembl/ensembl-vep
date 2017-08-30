### docker container for ensembl-vep
FROM ubuntu

# update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install apache2 build-essential cpanminus curl git libmysqlclient-dev libpng12-dev libssl-dev manpages mysql-client openssl perl perl-base unzip vim wget
# install ensembl dependencies
RUN cpanm DBI DBD::mysql

# create vep user
RUN useradd -r -m -U -d /home/vep -s /bin/bash -c "VEP User" -p '' vep
RUN usermod -a -G sudo vep
USER vep
ENV HOME /home/vep
WORKDIR $HOME

# clone git repositories
RUN mkdir -p src
WORKDIR $HOME/src
RUN git clone https://github.com/Ensembl/ensembl.git
RUN git clone https://github.com/Ensembl/ensembl-vep.git

# get VEP dependencies
WORKDIR $HOME/src
RUN ensembl-vep/travisci/get_dependencies.sh
ENV PERL5LIB $PERL5LIB:$HOME/src/bioperl-live-release-1-6-924
ENV KENT_SRC $HOME/src/kent-335_base/src
ENV HTSLIB_DIR $HOME/src/htslib
ENV MACHTYPE x86_64
ENV CFLAGS "-fPIC"
ENV DEPS $HOME/src

# and run the complilation/install as root
USER root
RUN ensembl-vep/travisci/build_c.sh

# install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# install bioperl-ext, faster alignments for haplo
WORKDIR $HOME/src
RUN git clone https://github.com/bioperl/bioperl-ext.git
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
RUN perl Makefile.PL
RUN make
RUN make install

# install perl dependencies
WORKDIR $HOME/src
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile .

# switch back to vep user
USER vep

# update bash profile
RUN echo >> $HOME/.profile && \
echo PATH=$HOME/src/ensembl-vep:\$PATH >> $HOME/.profile && \
echo export PATH >> $HOME/.profile

# setup environment
ENV PATH $HOME/src/ensembl-vep:$PATH

# run INSTALL.pl
WORKDIR $HOME/src/ensembl-vep
RUN chmod u+x *.pl
RUN ./INSTALL.pl -a a -l
