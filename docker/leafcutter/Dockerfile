# Base image for leafcutter,
# uses focal due to the need
# for python2 and R/3.6
FROM ubuntu:20.04

# Depedencies of leafcutter:
# • https://davidaknowles.github.io/leafcutter/articles/Installation.html
# @TLDR:
#    • samtools       apt-get, installs: 1.10-3
#    • regtools       build from source, install: v1.0.0
#    • python<=2.7    apt-get, installs: 2.7.17-2
#    • R<=3.6.0       apt-get, installs: 3.6.3-2
#        • rstan      apt-get, installs: 2.19.2-1
#        • devtools   apt-get, installs: 2.2.2-2
#        • leafcutter (using devtools from github)
LABEL maintainer=kuhnsa@nih.gov

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        locales \
        make \
        python2 \
        r-base \
        r-bioc-dirichletmultinomial \
        r-cran-devtools \
        r-cran-hmisc \
        r-cran-intervals \
        r-cran-optparse \
        r-cran-rstan \
        r-cran-rstantools \
        r-cran-r.utils \
        r-cran-shinycssloaders \
        r-cran-shinyjs \
        samtools \
        unzip \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
RUN ln -sf /usr/bin/python2 /usr/bin/python

############### INSTALL ################
# Install GNU which, increases compatiability
# with which distros included with Debian/Rocky
# linux installations, this version of which
# includes extra options that the ubuntu version 
# does not have. The gnu verison is backwards 
# compatiable with the ubunutu version.
RUN mkdir -p /opt2/ \
        && wget --no-check-certificate -O /opt2/which.tar.gz https://ftp.gnu.org/gnu/which/which-2.21.tar.gz \
        && tar -xvf /opt2/which.tar.gz -C /opt2/ \
        && rm /opt2/which.tar.gz \
        && cd /opt2/which-2.21/ \
        && ./configure --prefix=/opt2/which-2.21 \
        && make \
        && make install

# Install regtools from Github,
# https://github.com/griffithlab/regtools
RUN cd /opt2 \
        && wget https://github.com/griffithlab/regtools/archive/refs/tags/1.0.0.tar.gz -O /opt2/regtools.tar.gz \
        && tar -xvf /opt2/regtools.tar.gz -C /opt2/ \
        && rm /opt2/regtools.tar.gz \
        && cd /opt2/regtools-1.0.0/ \
        && mkdir build \
        && cd build \
        && cmake .. \
        && make 

# Install leafcutter from Github,
# https://davidaknowles.github.io/leafcutter/articles/Installation.html
RUN cd /opt2 \
        && git clone https://github.com/davidaknowles/leafcutter.git \
        && Rscript -e 'install.packages(c("TailRank"), Ncpus = 4, repos="http://cran.r-project.org")' \
        && Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

################ POST #################
# Add Dockerfile and export important 
# environment variables
ADD Dockerfile /opt2/leafcutter_v0-1-0.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2/which-2.21/bin:/opt2/regtools-1.0.0/build:/opt2/leafcutter/clustering:/opt2/leafcutter/scripts:/opt2/leafcutter/leafviz:${PATH}"
WORKDIR /data2