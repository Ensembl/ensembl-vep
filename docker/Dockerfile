###################################################
# Stage 1 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:16.04 as builder

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    git \
    libpng-dev \
    perl \
    perl-base \
    unzip \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV HTSLIB_DIR $OPT_SRC/htslib

# Working directory
WORKDIR $OPT_SRC

# Clone/download repositories/libraries
# Clone ensembl git repository and extract useful ensemb core file
RUN git clone --depth 1 https://github.com/Ensembl/ensembl.git && \
    cp ensembl/cpanfile ensembl_cpanfile && \
    rm -rf ensembl && \
    # Clone ensembl-vep git repository
    git clone --depth 1 https://github.com/Ensembl/ensembl-vep.git && chmod u+x ensembl-vep/*.pl && \
    # Clone ensembl-variation git repository and compile C code
    git clone --depth 1 https://github.com/Ensembl/ensembl-variation.git && \
    mkdir var_c_code && \
    cp ensembl-variation/C_code/*.c ensembl-variation/C_code/Makefile var_c_code/ && \
    rm -rf ensembl-variation && \
    chmod u+x var_c_code/* && \
    # Clone bioperl-ext git repository
    git clone --depth 1 https://github.com/bioperl/bioperl-ext.git && \
    # Download ensembl-xs
    wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip && \
    unzip -q ensembl-xs.zip && mv ensembl-xs-2.3.2 ensembl-xs && rm -rf ensembl-xs.zip && \
    # Clone/Download other repositories:  bioperl-live is needed so the cpanm dependencies installation from the ensembl-vep/cpanfile file takes less disk space
    ensembl-vep/travisci/get_dependencies.sh && mv bioperl-live-release-1-6-924 bioperl-live

# A lot of cleanup on the imported libraries, in order to reduce the docker image
RUN rm -rf Bio-HTS/.git Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci \
           bioperl-ext/.git bioperl-ext/README* bioperl-ext/t bioperl-ext/examples \
           bioperl-live/Changes bioperl-live/DEPENDENCIES bioperl-live/doc bioperl-live/examples bioperl-live/ide \
           bioperl-live/maintenance bioperl-live/scripts bioperl-live/t bioperl-live/travis_scripts \
           ensembl-vep/.git ensembl-vep/docker \
           ensembl-xs/.git ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/README ensembl-xs/t ensembl-xs/travisci \
           htslib/.git htslib/INSTALL htslib/NEWS htslib/README* htslib/test && \
    # Only keep needed kent-335_base libraries for VEP
    mv kent-335_base kent-335_base_bak && mkdir -p kent-335_base/src && \
    cp -R kent-335_base_bak/build kent-335_base_bak/confs kent-335_base/ && \
    cp -R kent-335_base_bak/src/lib kent-335_base_bak/src/inc kent-335_base_bak/src/jkOwnLib kent-335_base/src/ && \
    cp kent-335_base_bak/src/*.sh kent-335_base/src/ && \
    rm -rf kent-335_base_bak

# Setup bioperl-ext
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL

# Install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# Compile Variation LD C scripts
WORKDIR $OPT_SRC/var_c_code
RUN make && rm -f Makefile *.c


###################################################
# Stage 2 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:16.04

# Update aptitude and install some required packages
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

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PERL5LIB $PERL5LIB:$OPT_SRC/bioperl-live
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV DEPS $OPT_SRC
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Create vep user
RUN useradd -r -m -U -d $OPT -s /bin/bash -c "VEP User" -p '' vep && usermod -a -G sudo vep && mkdir -p $OPT_SRC
USER vep

# Copy downloaded libraries (stage 1) to this image (stage 2)
COPY --chown=vep:vep --from=builder $OPT_SRC $OPT_SRC
#############################################################

# Change user to root for the following complilations/installations
USER root


WORKDIR $OPT_SRC
# Install/compile some C libraries
RUN ensembl-vep/travisci/build_c.sh && \
    # Install ensembl perl dependencies (cpanm)
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl_cpanfile . && \
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile . && \
    # Configure "locale", see https://github.com/rocker-org/rocker/issues/19
    echo "$LANG UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=$LANG

# Install bioperl-ext, faster alignments for haplo
WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
RUN perl Makefile.PL && make && make install && rm -f Makefile.PL

# Install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT_SRC/ensembl-xs
RUN perl Makefile.PL && make && make install && rm -f Makefile.PL

# Switch back to vep user
USER vep

# Final steps
WORKDIR $OPT_SRC/ensembl-vep
# Update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile && \
    # Run INSTALL.pl and remove the ensemb-vep tests and travis
    ./INSTALL.pl -a a -l && rm -rf t travisci .travis.yml

