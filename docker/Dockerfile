#### docker container to build ensembl-vep - first stage ####
FROM alpine:3.8 as builder

# Add missing library to build first stage
RUN apk --no-cache add \
    perl \
    bash \
    git \
    unzip \
    wget

# Working directory
WORKDIR /opt/vep/src

# Clone ensembl git repository and extract useful ensemb core file
RUN git clone --depth 1 https://github.com/Ensembl/ensembl.git && \
    cp ensembl/cpanfile ensembl_cpanfile && \
    rm -rf ensembl

# Clone ensembl-vep git repository
RUN git clone --depth 1 https://github.com/Ensembl/ensembl-vep.git && chmod u+x ensembl-vep/*.pl

# Download bioperl-ext
RUN git clone --depth 1 https://github.com/bioperl/bioperl-ext.git
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL

# Download ensembl-xs
WORKDIR /opt/vep/src
RUN wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip && \
    unzip -q ensembl-xs.zip && mv ensembl-xs-2.3.2 ensembl-xs && rm -rf ensembl-xs.zip

# Clone/Download other repositories
# bioperl-live is needed so the cpanm dependencies installation from the ensembl-vep/cpanfile file takes less disk space
RUN ensembl-vep/travisci/get_dependencies.sh && mv bioperl-live-release-1-6-924 bioperl-live

# A lot of cleanup on the imported libraries, in order to reduce the docker image
RUN rm -rf Bio-HTS/.git Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci \
           bioperl-ext/.git bioperl-ext/README* bioperl-ext/t bioperl-ext/examples \
           bioperl-live/Changes bioperl-live/DEPENDENCIES bioperl-live/doc bioperl-live/examples bioperl-live/ide \
           bioperl-live/maintenance bioperl-live/scripts bioperl-live/t bioperl-live/travis_scripts \
           ensembl-vep/.git ensembl-vep/docker \
           ensembl-xs/.git ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/README ensembl-xs/t ensembl-xs/travisci \
           htslib/.git htslib/INSTALL htslib/NEWS htslib/README* htslib/test \
           kent-335_base/README kent-335_base/python kent-335_base/java



#### docker container to run ensembl-vep  - second stage ####
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
    mysql-client \
    openssl \
    perl \
    perl-base \
    unzip \
    vim && \
    rm -rf /var/lib/apt/lists/*

# create vep user
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
RUN useradd -r -m -U -d $OPT -s /bin/bash -c "VEP User" -p '' vep && usermod -a -G sudo vep && mkdir -p $OPT_SRC
USER vep

# Copy downloaded libraries (first stage) to this image (second stage)
COPY --chown=vep:vep --from=builder $OPT_SRC $OPT_SRC

# Setup VEP environment
WORKDIR $OPT_SRC
ENV PERL5LIB $PERL5LIB:$OPT_SRC/bioperl-live
ENV KENT_SRC $OPT_SRC/kent-335_base/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV MACHTYPE x86_64
ENV CFLAGS "-fPIC"
ENV DEPS $OPT/src
ENV PATH $OPT_SRC/ensembl-vep:$PATH

# update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$OPT_SRC/ensembl-vep:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile

# run the complilation/install as root
USER root
RUN ensembl-vep/travisci/build_c.sh

# install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# install bioperl-ext, faster alignments for haplo
WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
RUN perl Makefile.PL && make && make install

# install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT_SRC/ensembl-xs
RUN perl Makefile.PL && make && make install

# install ensembl perl dependencies
WORKDIR $OPT_SRC
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl_cpanfile . && \
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile .

# configure locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# switch back to vep user
USER vep

# run INSTALL.pl and remove the ensemb-vep tests and travis
WORKDIR $OPT_SRC/ensembl-vep
RUN ./INSTALL.pl -a a -l && rm -rf t travisci .travis.yml

