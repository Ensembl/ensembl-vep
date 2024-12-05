FROM ubuntu:22.04

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# Install Linux Dependencies
RUN apt-get update && apt-get -y --no-install-recommends install \
    build-essential \
    ca-certificates \
    cpanminus \
    git \
    libdbd-sqlite3-perl \
    liblzma-dev \
    libmysqlclient-dev \
    libncurses5-dev \
    libpng-dev \
    libssl-dev \
    libbz2-dev \
    locales \
    openssl \
    procps \
    unzip \
    wget \
    zlib1g-dev \
    # Installing packages from https://raw.githubusercontent.com/Ensembl/VEP_plugins/$BRANCH/config/ubuntu-packages.txt
    pkg-config \
    libgd-dev \
    uuid-dev \
    sqlite3 \
    python2 \
    python-pip \
    python-setuptools \
    python2-dev && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Perl from source
ENV PERL_VERSION="5.40.0"

RUN wget -q "https://www.cpan.org/src/5.0/perl-${PERL_VERSION}.tar.gz" -O /opt/perl.tar.gz && \
    tar -xvzf /opt/perl.tar.gz -C /opt && \
    cd "/opt/perl-${PERL_VERSION}" && \
    ./Configure -des -Dprefix="/usr/local/perl-${PERL_VERSION}" && \
    make && make install && make clean && \
    rm -rf /opt/perl.tar.gz

ENV OPT="/opt/vep"
ENV OPT_SRC="${OPT}/src"

# Ensemble release version
ENV VEP_VERSION="113"
ENV BRANCH="release/${VEP_VERSION}"
ENV VEP_DIR_PLUGINS="${OPT}/.vep/Plugins"

# Make plugins directory
RUN mkdir -p "${VEP_DIR_PLUGINS}"

# Install Ensembl VEP
RUN git clone -b "${BRANCH}" --depth 1 https://github.com/Ensembl/ensembl-vep.git "${OPT_SRC}/ensembl-vep" && \
    # Clean up of the ensembl-vep libraries to reduce the docker image
    cd "${OPT_SRC}" && \
    rm -rf ensembl-vep/.??* ensembl-vep/docker && \
    chmod u+x "${OPT_SRC}/ensembl-vep"/*.pl && \
    rm ensembl-vep/cpanfile

# Install Bioperl ext
ENV BIOPERL_EXT_VERSION="1-5-1"

RUN wget -q "https://github.com/bioperl/bioperl-ext/archive/refs/tags/bioperl-ext-release-${BIOPERL_EXT_VERSION}.tar.gz" -O /bioperl-ext.tar.gz && \
    tar -xvzf /bioperl-ext.tar.gz -C "${OPT_SRC}" && \
    mv "${OPT_SRC}/bioperl-ext-bioperl-ext-release-${BIOPERL_EXT_VERSION}" "${OPT_SRC}/bioperl-ext" && \
    rm /bioperl-ext.tar.gz && \
    cd "${OPT_SRC}/bioperl-ext/Bio/Ext/Align" && \
    # Adds the fPIC compilation options to the CFLAG environment variable in Makefile.PL
    perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL && \
    perl Makefile.PL && \
    make && make install && make clean && \
    rm Makefile* && \
    # Cleanup of the bioperl libraries to reduce the docker image
    cd "${OPT_SRC}" && \
    rm -rf bioperl-ext/.??* bioperl-ext/Bio/SeqIO bioperl-ext/Bio/Tools bioperl-ext/Makefile.PL bioperl-ext/README* bioperl-ext/t bioperl-ext/examples

# Install Ensembl XS
ENV ENSEMBL_XS_VERSION="2.3.2"

RUN wget -q "https://github.com/Ensembl/ensembl-xs/archive/${ENSEMBL_XS_VERSION}.zip" -O /ensembl-xs.zip && \
    unzip -q /ensembl-xs.zip && \
    mv "ensembl-xs-${ENSEMBL_XS_VERSION}" "${OPT_SRC}/ensembl-xs" && \
    cd "${OPT_SRC}/ensembl-xs" && \
    perl Makefile.PL && \
    make && make install && make clean && \
    rm Makefile* cpanfile && \
    rm /ensembl-xs.zip && \
    cd "${OPT_SRC}" && \
    rm -rf ensembl-xs/.??* ensembl-xs/TODO ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/MANIFEST ensembl-xs/README ensembl-xs/t ensembl-xs/travisci

# Reproduce https://github.com/Ensembl/ensembl-vep/blob/release/112/travisci/get_dependencies.sh

# Install bioperl Bio module (this is deleted later)
ENV BIOPERL_VERSION="1-6-924"

RUN git clone --branch "release-${BIOPERL_VERSION}" --depth 1 https://github.com/bioperl/bioperl-live.git /bioperl-to-delete && \
    mkdir -p "${OPT_SRC}/bioperl-live" && \
    mv /bioperl-to-delete/Bio "${OPT_SRC}/bioperl-live/" && \
    rm -rf /bioperl-to-delete

# Install samtools and htslib
ENV SAMTOOLS_DIR="${OPT_SRC}/samtools"

ENV SAMTOOLS_VERSION="1.9"

ENV HTSLIB_DIR="${OPT_SRC}/samtools/htslib-${SAMTOOLS_VERSION}"

RUN mkdir -p "${SAMTOOLS_DIR}" && \
    wget -q "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" -O /opt/samtools.tar.bz2 && \
    tar -xvjf /opt/samtools.tar.bz2 -C "${OPT_SRC}" && \
    mv "${OPT_SRC}/samtools-${SAMTOOLS_VERSION}"/* "${SAMTOOLS_DIR}" && \
    cd "${SAMTOOLS_DIR}" && \
    ./configure --prefix="${SAMTOOLS_DIR}" && \
    make && make install && make clean && \
    cd "${HTSLIB_DIR}" && \
    make install && \
    rm -rf Makefile ./*.c ./.??* INSTALL NEWS README* test && \
    rm -rf "${OPT_SRC}/samtools-${SAMTOOLS_VERSION}" && \
    rm /opt/samtools.tar.bz2

ENV PATH="${SAMTOOLS_DIR}/bin:${HTSLIB_DIR}:${PATH}"

# Install Ensembl Variation LD C Scripts
# Has to occur after samtools/htslib_dir is configured
RUN git clone -b "${BRANCH}" --depth 1 https://github.com/Ensembl/ensembl-variation.git /ensembl-variation-to-delete && \
    mkdir "${OPT_SRC}/var_c_code" && \
    cp /ensembl-variation-to-delete/C_code/*.c /ensembl-variation-to-delete/C_code/Makefile "${OPT_SRC}/var_c_code/" && \
    chmod u+x "${OPT_SRC}/var_c_code"/* && \
    cd "${OPT_SRC}/var_c_code" && \
    make && \
    rm -f Makefile ./*.c && \
    rm -rf /ensembl-variation-to-delete

# Install BIO-DB-HTS
ENV BIO_DB_HTS_VERSION="2.11"

RUN wget -q "https://github.com/Ensembl/Bio-DB-HTS/archive/refs/tags/${BIO_DB_HTS_VERSION}.tar.gz" -O /opt/bio_db_hts.tar.gz && \
    tar -xvzf /opt/bio_db_hts.tar.gz -C ${OPT_SRC} && \
    mv "${OPT_SRC}/Bio-DB-HTS-${BIO_DB_HTS_VERSION}" "${OPT_SRC}/Bio-HTS" && \
    rm /opt/bio_db_hts.tar.gz && \
    cd "${OPT_SRC}" && \
    rm -rf Bio-HTS/.??* Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/MANIFEST* Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci

# Install Kent base libraries
ENV KENT_VERSION="335"

ENV KENT_SRC="${OPT_SRC}/kent-${KENT_VERSION}_base/src"

RUN wget -q "https://github.com/ucscGenomeBrowser/kent/archive/v${KENT_VERSION}_base.tar.gz" -O /opt/kent.tar.gz && \
    tar -xvzf /opt/kent.tar.gz -C "${OPT_SRC}" && \
    rm -rf "${OPT_SRC}/kent-${KENT_VERSION}_base/java" "${OPT_SRC}/kent-${KENT_VERSION}_base/python" && \
    mv "${OPT_SRC}/kent-${KENT_VERSION}_base" /kent-to-delete && \
    mkdir -p "${KENT_SRC}" && \
    cp -R /kent-to-delete/src/lib /kent-to-delete/src/inc /kent-to-delete/src/jkOwnLib "${KENT_SRC}" && \
    cp /kent-to-delete/src/*.sh "${KENT_SRC}" && \
    rm -rf /kent-to-delete && \
    rm /opt/kent.tar.gz

# Install loftee
ENV LOFTEE_GRCH38_TAG="1.0.4_GRCh38"

RUN wget -q "https://github.com/konradjk/loftee/archive/refs/tags/v${LOFTEE_GRCH38_TAG}.tar.gz" -O /loftee.GRCh38.tar.gz && \
    tar -xvzf /loftee.GRCh38.tar.gz -C /opt && \
    mv "/opt/loftee-${LOFTEE_GRCH38_TAG}" "${VEP_DIR_PLUGINS}" && \
    rm /loftee.GRCh38.tar.gz

ENV PERL5LIB_TMP="${OPT_SRC}/ensembl-vep:${OPT_SRC}/ensembl-vep/modules"
ENV PERL5LIB="${PERL5LIB_TMP}:${OPT_SRC}/bioperl-live"
ENV DEPS="${OPT_SRC}"
ENV PATH="${OPT_SRC}/ensembl-vep:${OPT_SRC}/var_c_code:${PATH}"
ENV LANG_VAR=en_US.UTF-8

# INSTALLING ALL PACKAGES FROM "https://raw.githubusercontent.com/Ensembl/ensembl/$BRANCH/cpanfile" manually
RUN cpanm --notest Bio::Perl@1.6.924 \
    DBI@1.643

## ISSUE: https://forum.bestpractical.com/t/mysql-dependency-error-with-mariadb-and-debian-12/38748/3
## Need to install a version of mysql < 5.X to work with mariadb
## Match version in: https://github.com/Ensembl/ensembl/blob/release/$BRANCH/cpanfile
RUN cpanm --notest DBD::mysql@4.050 \
    HTTP::Tiny@0.088 \
    IO::Compress::Gzip@2.213 \
    URI::Escape@5.29 \
    Config::IniFiles@3.000003 \
    Gzip::Faster@0.21 \
    List::MoreUtils@0.430 \
    Algorithm::Diff@1.201 \
    Tie::IxHash@1.23 \
    DBIx::Class@0.082843 \
    Digest::MD5@2.59 \
    LWP::UserAgent@6.77 \
    List::Util@1.68 \
    Moose@2.2207 \
    SQL::Translator@1.65 \
    Text::CSV@2.04 \
    Text::CSV_XS@1.56 \
    Text::Glob@0.11 \
    URI@5.29 \
    XML::LibXML@2.0210 \
    Text::Unidecode@1.30

#INSTALLING ALL PACKAGED FROM "https://raw.githubusercontent.com/Ensembl/ensembl-vep/$BRANCH/cpanfile"
RUN cpanm --notest Set::IntervalTree@0.12 \
    JSON@4.10 \
    PerlIO::gzip@0.20 \
    IO::Uncompress::Gunzip@2.213 \
    Sereal@5.004 \
    HTML::Lint@2.32 \
    Capture::Tiny@0.48
## INSTRUCTIONS: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#bigfile
## The following code needs to be within 1 RUN block so the exported variables are captured
RUN MACHTYPE=$(uname -m) \
    MYSQLINC="$(mysql_config --include | sed -e 's/^-I//g')" && \
    MYSQLLIBS="$(mysql_config --libs)" && \
    export MYSQLINC MYSQLLIBS && \
    cd "${OPT_SRC}/kent-${KENT_VERSION}_base/src/lib" && \
    echo 'CFLAGS="-fPIC"' > "${OPT_SRC}/kent-${KENT_VERSION}_base/src/inc/localEnvironment.mk" && \
    make && \
    cd "${OPT_SRC}/kent-${KENT_VERSION}_base/src/jkOwnLib" && \
    make && \
    ln -s "${OPT_SRC}/kent-${KENT_VERSION}_base/src/lib/x86_64/"* "${OPT_SRC}/kent-${KENT_VERSION}_base/src/lib" && \
    cpanm --notest Bio::DB::BigFile@1.07 && \
    cd "${OPT_SRC}/kent-${KENT_VERSION}_base/src/lib" && make clean && \
    cd "${OPT_SRC}/kent-${KENT_VERSION}_base/src/jkOwnLib" && make clean
RUN cpanm --notest DBD::SQLite@1.76

# Delete bioperl-live install as no longer needed after cpan installations
RUN cd "${OPT_SRC}" && \
    rm -rf bioperl-live

# Configure "locale", see https://github.com/rocker-org/rocker/issues/19
RUN echo "${LANG_VAR} UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=${LANG_VAR}

ENV LC_ALL="${LANG_VAR}"
ENV LANG="${LANG_VAR}"
# Bioperl-live is deleted so need to update PERL5LIB
ENV PERL5LIB="${PERL5LIB_TMP}"
ENV PLUGIN_URL="https://raw.githubusercontent.com/Ensembl/VEP_plugins"

# INSTALL EXTRA DEPENDENCIES FOR FOLLOWING LAYER
# Needs to go after BigFile installation for BioPerl to be already correctly installed
RUN cpanm --notest Bio::DB::HTS::Tabix@3.01
RUN cpanm --notest Test::Warnings@0.033

ENV SUPPORTED_VEP_PLUGINS_LIST="AVADA,AlphaMissense,AncestralAllele,BayesDel,Blosum62,CADD,CAPICE,CCDSFilter,CSN,Carol,ClinPred,Condel,Conservation,DAS,DeNovo,DisGeNET,DosageSensitivity,Downstream,Draw,EVE,Enformer,FATHMM,FATHMM_MKL,FlagLRG,FunMotifs,G2P,GO,GWAS,GXA,GeneBe,GeneSplicer,Geno2MP,HGVSIntronOffset,HGVSReferenceBase,IntAct,LD,LOEUF,LOVD,LoFtool,LocalID,MPC,MTR,Mastermind,MaveDB,MaxEntScan,NMD,NearestExonJB,NearestGene,NonSynonymousFilter,OpenTargets,PON_P2,Paralogues,PhenotypeOrthologous,Phenotypes,PolyPhen_SIFT,PostGAP,PrimateAI,ProteinSeqs,REVEL,RankFilter,RefSeqHGVS,ReferenceQuality,RiboseqORFs,SameCodon,SingleLetterAA,SpliceAI,SpliceRegion,SpliceVault,StructuralVariantOverlap,SubsetVCF,TSSDistance,TranscriptAnnotator,UTRAnnotator,VARITY,dbNSFP,dbscSNV,gnomADc,mutfunc,neXtProt,pLI,satMutMPRA"

RUN cd "${OPT_SRC}/ensembl-vep" && \
    echo PATH=$PATH:\$PATH >> "${OPT}/.profile" && \
    echo export PATH >> "${OPT}/.profile" && \
    # Run INSTALL.pl and remove the ensemb-vep tests and travis
    # Gwava plugin fails to install currently: https://github.com/Ensembl/ensembl-vep/issues/931
    ./INSTALL.pl \
        -a ap \
        -l -n \
        --pluginsdir "${VEP_DIR_PLUGINS}" \
        --PLUGINURL "${PLUGIN_URL}" \
        -g "${SUPPORTED_VEP_PLUGINS_LIST}" && \
    rm -rf t travisci .travis.yml

# INSTALLING packages from https://raw.githubusercontent.com/Ensembl/VEP_plugins/$BRANCH/config/cpanfile
# DBD::SQLite already installed
RUN cpanm --notest GD@2.77
RUN cpanm --notest Math::CDF@0.1

# Install GeneSplicer binary
RUN wget "ftp://ftp.ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz" && \
    tar -xzf GeneSplicer.tar.gz -C "${VEP_DIR_PLUGINS}" && \
    rm GeneSplicer.tar.gz && \
    cd ${VEP_DIR_PLUGINS}/GeneSplicer/sources && \
    make && \
    mv genesplicer .. && \
    rm -rf GeneSplicer/*/

ENV PATH="${VEP_DIR_PLUGINS}/GeneSplicer:${PATH}"

ENV PERL5LIB="${VEP_DIR_PLUGINS}:${VEP_DIR_PLUGINS}/loftee-${LOFTEE_GRCH38_TAG}:${PERL5LIB}"

# Set working directory as symlink to $OPT/.vep (containing VEP cache and data)
RUN ln -s "${OPT}/.vep" /data
WORKDIR /data

COPY Dockerfile /Dockerfile
