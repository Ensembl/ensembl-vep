### docker container for ensembl-vep
FROM ubuntu:16.04

# update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    apache2 \
    build-essential \
    cpanminus \
    curl \
    git \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    locales \
    manpages \
    mysql-client \
    openssl \
    perl \
    perl-base \
    unzip \
    vim \
    wget

# install ensembl dependencies
RUN cpanm DBI DBD::mysql

# create vep user
RUN useradd -r -m -U -d /opt/vep -s /bin/bash -c "VEP User" -p '' vep
RUN usermod -a -G sudo vep
USER vep
ENV OPT /opt/vep
WORKDIR $OPT

# clone git repositories
RUN mkdir -p src
WORKDIR $OPT/src
RUN git clone https://github.com/Ensembl/ensembl.git
RUN git clone https://github.com/Ensembl/ensembl-vep.git

# get VEP dependencies
WORKDIR $OPT/src
RUN ensembl-vep/travisci/get_dependencies.sh
ENV PERL5LIB $PERL5LIB:$OPT/src/bioperl-live-release-1-6-924
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT/src/htslib
ENV MACHTYPE x86_64
ENV CFLAGS "-fPIC"
ENV DEPS $OPT/src

# and run the complilation/install as root
USER root
RUN ensembl-vep/travisci/build_c.sh

# install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# install bioperl-ext, faster alignments for haplo
WORKDIR $OPT/src
RUN git clone https://github.com/bioperl/bioperl-ext.git
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
RUN perl Makefile.PL
RUN make
RUN make install

# install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT/src
RUN wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip
RUN unzip -q ensembl-xs.zip
RUN mv ensembl-xs-2.3.2 ensembl-xs
RUN rm ensembl-xs.zip
WORKDIR ensembl-xs
RUN perl Makefile.PL
RUN make
RUN make install

# install perl dependencies
WORKDIR $OPT/src
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile .

# configure locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# switch back to vep user
USER vep

# update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$OPT/src/ensembl-vep:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile

# setup environment
ENV PATH $OPT/src/ensembl-vep:$PATH

# run INSTALL.pl
WORKDIR $OPT/src/ensembl-vep
RUN chmod u+x *.pl
RUN ./INSTALL.pl -a a -l
