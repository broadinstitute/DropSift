# This is borrowing a bit from existing dockerfile code to use the same OS and R version, but
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install OS stuff that R will need later.
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    software-properties-common \
    dirmngr \
    build-essential \
    gfortran \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libblas-dev \
    liblapack-dev \
    zlib1g-dev \
    libcairo2-dev \
    libpng-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    libfreetype6-dev \
    libjpeg-dev \
    libwebp-dev \
    libuv1-dev \
    cmake \
    pkg-config \
 && rm -rf /var/lib/apt/lists/*

#NOTE: I had a lot of trouble installing a pinned version of R.

RUN . /etc/lsb-release && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu ${DISTRIB_CODENAME}-cran40/" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        r-base \
 && rm -rf /var/lib/apt/lists/*

# Set the repository for all R calls.
RUN echo "options(repos = c(CRAN = 'https://cloud.r-project.org'))" > /etc/R/Rprofile.site
RUN R -e "install.packages('remotes')"

# Install the code.
WORKDIR /app
COPY . /app

RUN R -e "remotes::install_local('/app')"

#Add the testthat package so a user can run validation
RUN R -e "install.packages(c('testthat'))"

#Add presto because it speeds up some parts of the statistics
RUN R -e "remotes::install_github('immunogenomics/presto')"
