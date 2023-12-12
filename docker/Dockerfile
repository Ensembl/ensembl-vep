ARG BRANCH=main

###################################################
# Stage 1 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:23.04 as builder
ARG DEBIAN_FRONTEND=noninteractive

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
# many dependencies are needed to ensure that that vulnerabilities are squashed
RUN apt-get update && apt-get -y install --no-install-recommends \
    build-essential \
    git \
    libpng-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    perl \
    perl-base \
    python3 \
    python3-pip \
    unzip \
    wget && \
    mv /usr/lib/python3.11/EXTERNALLY-MANAGED /usr/lib/python3.11/EXTERNALLY-MANAGED.old && \
    python3 -m pip install --no-cache-dir --upgrade pip && \
    python3 -m pip uninstall cryptography && \
    python3 -m pip install --no-cache-dir cryptography==41.0.4 setuptools==65.5.1 && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV BRANCH

# Working directory
WORKDIR $OPT_SRC

# Clone/download repositories/libraries
RUN if [ "$BRANCH" = "main" ]; \
    then export BRANCH_OPT=""; \
    else export BRANCH_OPT="-b $BRANCH"; \
    fi && \
    # Get ensembl cpanfile in order to get the list of the required Perl libraries
    wget -q "https://raw.githubusercontent.com/Ensembl/ensembl/$BRANCH/cpanfile" -O "ensembl_cpanfile" && \
    # Clone ensembl-vep git repository
    git clone $BRANCH_OPT --depth 1 https://github.com/Ensembl/ensembl-vep.git && chmod u+x ensembl-vep/*.pl && \
    # Clone ensembl-variation git repository and compile C code
    git clone $BRANCH_OPT --depth 1 https://github.com/Ensembl/ensembl-variation.git && \
    mkdir var_c_code && \
    cp ensembl-variation/C_code/*.c ensembl-variation/C_code/Makefile var_c_code/ && \
    rm -rf ensembl-variation && \
    chmod u+x var_c_code/* && \
    # Clone bioperl-ext git repository - used by Haplosaurus
    git clone --depth 1 https://github.com/bioperl/bioperl-ext.git && \
    # Download ensembl-xs - it contains compiled versions of certain key subroutines used in VEP
    wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip && \
    unzip -q ensembl-xs.zip && mv ensembl-xs-2.3.2 ensembl-xs && rm -rf ensembl-xs.zip && \
    # Clone/Download other repositories: bioperl-live is needed so the cpanm dependencies installation from the ensembl-vep/cpanfile file takes less disk space
    ensembl-vep/travisci/get_dependencies.sh && \
    # Only keep the bioperl-live"Bio" library
    mv bioperl-live bioperl-live_bak && mkdir bioperl-live && mv bioperl-live_bak/Bio bioperl-live/ && rm -rf bioperl-live_bak && \
    ## A lot of cleanup on the imported libraries, in order to reduce the docker image ##
    rm -rf Bio-HTS/.??* Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/MANIFEST* Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci \
           bioperl-ext/.??* bioperl-ext/Bio/SeqIO bioperl-ext/Bio/Tools bioperl-ext/Makefile.PL bioperl-ext/README* bioperl-ext/t bioperl-ext/examples \
           ensembl-vep/.??* ensembl-vep/docker \
           ensembl-xs/.??* ensembl-xs/TODO ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/MANIFEST ensembl-xs/README ensembl-xs/t ensembl-xs/travisci \
           htslib/.??* htslib/INSTALL htslib/NEWS htslib/README* htslib/test && \
    # Only keep needed kent-335_base libraries for VEP - used by Bio::DB::BigFile (bigWig parsing)
    mv kent-335_base kent-335_base_bak && mkdir -p kent-335_base/src && \
    cp -R kent-335_base_bak/src/lib kent-335_base_bak/src/inc kent-335_base_bak/src/jkOwnLib kent-335_base/src/ && \
    cp kent-335_base_bak/src/*.sh kent-335_base/src/ && \
    rm -rf kent-335_base_bak

# # Build LOFTEE for GRCh37 & GRRCh38
# # This section is included for completion but commented out as LOFTEE is not a cononical plugin
# RUN mkdir /loftee-GRCh37 /loftee-GRCh38 && \
#     cd /loftee-GRCh37 && \
#     git clone https://github.com/konradjk/loftee.git && \
#     cd /loftee-GRCh38 && \
#     git clone -b ${LOFTEE_GRCh38_BRANCH} https://github.com/konradjk/loftee.git

# Setup bioperl-ext
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL

# Install htslib binaries (for 'bgzip' and 'tabix')
# htslib requires the packages 'zlib1g-dev', 'libbz2-dev' and 'liblzma-dev'
WORKDIR $HTSLIB_DIR
RUN make install && rm -f Makefile *.c

# Compile Variation LD C scripts
WORKDIR $OPT_SRC/var_c_code
RUN make && rm -f Makefile *.c


###################################################
# Stage 2 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:23.04
ARG DEBIAN_FRONTEND=noninteractive

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
# some packages have been moved here form the ubuntu-packages.txt
# due to the need to clear vulnerabilities
RUN apt-get update && apt-get -y install --no-install-recommends \
    software-properties-common \
    build-essential \
    cpanminus \
    curl \
    libgd-dev \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    locales \
    openssl \
    perl \
    perl-base \
    pkg-config \
    python3 \
    python3-pip \
    sqlite3 \
    uuid-dev \
    unzip \
    vim && \
    # following lines are needed to install updates to compromiosed python packages
    mv /usr/lib/python3.11/EXTERNALLY-MANAGED /usr/lib/python3.11/EXTERNALLY-MANAGED.old && \
    python3 -m pip install --no-cache-dir --upgrade pip && \
    python3 -m pip uninstall cryptography && \
    python3 -m pip install --no-cache-dir cryptography==41.0.4 setuptools==65.5.1 && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PERL5LIB_TMP $PERL5LIB:$OPT_SRC/ensembl-vep:$OPT_SRC/ensembl-vep/modules
ENV PERL5LIB $PERL5LIB_TMP:$OPT_SRC/bioperl-live
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV MACHTYPE x86_64
ENV DEPS $OPT_SRC
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
ENV LANG_VAR en_US.UTF-8
ARG BRANCH

# Create vep user
RUN useradd -r -m -U -d "$OPT" -s /bin/bash -c "VEP User" -p '' vep && usermod -a -G sudo vep && mkdir -p $OPT_SRC
USER vep

# Copy downloaded libraries (stage 1) to this image (stage 2)
COPY --chown=vep:vep --from=builder $OPT_SRC $OPT_SRC
#############################################################

# Change user to root for the following complilations/installations
USER root

# Install bioperl-ext, faster alignments for haplo (XS-based BioPerl extensions to C libraries)
WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
RUN perl Makefile.PL && make && make install && rm -f Makefile*

# Install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT_SRC/ensembl-xs
RUN perl Makefile.PL && make && make install && rm -f Makefile* cpanfile

WORKDIR $OPT_SRC
# Install/compile more libraries
RUN export MACHTYPE=$(uname -m) &&\
    ensembl-vep/travisci/build_c.sh && \
    # Remove unused Bio-DB-HTS files
    rm -rf Bio-HTS/cpanfile Bio-HTS/Build.PL Bio-HTS/Build Bio-HTS/_build Bio-HTS/INSTALL.pl && \
    # Install ensembl perl dependencies (cpanm)
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl_cpanfile . && \
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile . && \
    # additional perl dependencies
    cpanm Try::Tiny Module::Build && \
    # Delete bioperl and cpanfiles after the cpanm installs as bioperl will be reinstalled by the INSTALL.pl script
    rm -rf bioperl-live ensembl_cpanfile ensembl-vep/cpanfile && \
    # Configure "locale", see https://github.com/rocker-org/rocker/issues/19
    echo "$LANG_VAR UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=$LANG_VAR && \
    # Copy htslib executables. It also requires the packages 'zlib1g-dev', 'libbz2-dev' and 'liblzma-dev'
    cp $HTSLIB_DIR/bgzip $HTSLIB_DIR/tabix $HTSLIB_DIR/htsfile /usr/local/bin/ && \
    # Remove CPAN cache
    rm -rf /root/.cpanm

ENV LC_ALL $LANG_VAR
ENV LANG $LANG_VAR

# Switch back to vep user
USER vep
ENV PERL5LIB $PERL5LIB_TMP
# Removing FATHMM, FATHMM_MKL & PON_P2 as these plugins require python 2.x which is End-Of-Life
ENV PLUGIN_LIST=AlphaMisssense,AncestralAllele,Blosum62,CADD,Carol,Condel,Conservation,dbNSFP,dbscSNV,DisGeNET,Downstream,Draw,ExAC,ExACpLi,FlagLRG,FunMotifs,G2P,GeneSplicer,gnomADc,GO,HGVSIntronOffset,LD,LocalID,LoFtool,LOVD,Mastermind,MaxEntScan,MPC,MTR,NearestExonJB,NearestGene,neXtProt,Phenotypes,PON_P2,PostGAP,ProteinSeqs,ReferenceQuality,REVEL,SameCodon,satMutMPRA,SingleLetterAA,SpliceAI,SpliceRegion,StructuralVariantOverlap,SubsetVCF,TSSDistance,UTRannotator

# Setup Docker environment for when users run VEP and INSTALL.pl in Docker image:
#   - skip VEP updates in INSTALL.pl
ENV VEP_NO_UPDATE 1
#   - avoid Faidx/HTSLIB installation in INSTALL.pl
ENV VEP_NO_HTSLIB 1
#   - skip plugin installation in INSTALL.pl
ENV VEP_NO_PLUGINS 1
#   - set plugins directory for VEP and INSTALL.pl
ENV VEP_DIR_PLUGINS /plugins
ENV VEP_PLUGINSDIR $VEP_DIR_PLUGINS
WORKDIR $VEP_DIR_PLUGINS

# Update bash profile
WORKDIR $OPT_SRC/ensembl-vep
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile && \
    # Install Ensembl API and plugins
    ./INSTALL.pl --auto ap --plugins all --pluginsdir ${VEP_DIR_PLUGINS} --no_update --no_htslib && \
    # Remove the ensemb-vep tests and travis
    rm -rf t travisci .travis.yml
    # # needed to prevent LoF compilation errors at runtime
    # # only required if LOFTEE is being intalled
    # rm -rf $VEP_DIR_PLUGINS/LoF.pm

# # Copy over LOFTEE build files
# # Only needed if LOFTEE is being installed
# COPY --from=builder /loftee-GRCh37/loftee /plugins/loftee-GRCh37
# COPY --from=builder /loftee-GRCh38/loftee /plugins/loftee-GRCh38

# Install dependencies for VEP plugins:
USER root
ENV PLUGIN_DEPS "https://raw.githubusercontent.com/Ensembl/VEP_plugins/$BRANCH/config"
#   - Ubuntu packages
# Commented out as the packages in question are explicitly installed above.
# Retained to ensure that future updates can more easily reintegrate the commands
# RUN curl -O "$PLUGIN_DEPS/ubuntu-packages.txt" && \
#     apt-get update && apt-get install -y --no-install-recommends \
#     "$(sed s/\#.*//g ubuntu-packages.txt)" && \
#     rm -rf /var/lib/apt/lists/* ubuntu-packages.txt
#   - Perl modules
RUN curl -O "$PLUGIN_DEPS/cpanfile" && \
    cpanm --installdeps --with-recommends . && \
    rm -rf /root/.cpanm cpanfile
#   - Python packages
# Commented out as Python2 has been removed from the container.
# Retained to ensure that future updates can more easily reintegrate the commands
# RUN curl -O "$PLUGIN_DEPS/requirements.txt" && \
#     python -m pip install --no-cache-dir -r requirements.txt && \
#     python -m pip uninstall urllib && \
#     rm requirements.txt

# Set working directory as symlink to $OPT/.vep (containing VEP cache and data)
USER root
RUN ln -s $OPT/.vep /data
USER vep
WORKDIR /data
