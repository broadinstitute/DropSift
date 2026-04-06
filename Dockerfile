# This is borrowing a bit from existing dockerfile code to use the same OS and R version, but
FROM ubuntu:24.04

ARG R_VERSION=4.5.2

ENV DEBIAN_FRONTEND=noninteractive

# Install OS stuff that R will need later.
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
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
    libglpk40 \
    cmake \
    pkg-config \
 && rm -rf /var/lib/apt/lists/*

# Use the pinned version of R
RUN bash -lc '. /etc/lsb-release && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc > /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${DISTRIB_CODENAME}-cran40/" > /etc/apt/sources.list.d/cran-r.list && \
    apt-get update && \
    r_version=${R_VERSION} && \
    r_base_version=${r_version}-1.${DISTRIB_RELEASE/./}.0 && \
    apt-get install -y --no-install-recommends \
        r-base-core=${r_base_version} \
        r-base-dev=${r_base_version} \
    && rm -rf /var/lib/apt/lists/*'


# Use the Posit Package Manager's Ubuntu binary repo
# Sets the HTTP user agent so R can request Linux binaries.
RUN . /etc/lsb-release && \
    cat >/etc/R/Rprofile.site <<EOF
local({
    repo <- getOption("repos")
    repo["CRAN"] <- "https://packagemanager.posit.co/cran/__linux__/${DISTRIB_CODENAME}/latest"
    options(
        repos = repo,
        download.file.method = "curl",
        download.file.extra = paste(
            "-fsSL",
            sprintf(
                '--header "User-Agent: R (%s)"',
                paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])
            )
        )
    )
})
EOF


# Install package dependencies in a separate cached layer.
RUN R -e "install.packages('remotes')"

WORKDIR /app
COPY DESCRIPTION /app/DESCRIPTION
RUN R -e "remotes::install_deps('/app', dependencies = TRUE)"

#Add the testthat package so a user can run validation
RUN R -e "install.packages(c('testthat'))"

#Add presto because it speeds up some parts of the statistics
RUN R -e "remotes::install_github('immunogenomics/presto')"

# Copy the rest of the package source after dependencies are installed.
COPY . /app

# Install the package itself without reinstalling dependencies.
RUN R -e "remotes::install_local('/app', dependencies = FALSE)"