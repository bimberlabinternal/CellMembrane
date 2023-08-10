# Note: this is the last base version supporting ubuntu focal, not jammy
FROM rocker/rstudio:4.2.1

ARG GH_PAT='NOT_SET'

##  Add Bioconductor system dependencies
RUN wget -O install_bioc_sysdeps.sh https://raw.githubusercontent.com/Bioconductor/bioconductor_docker/master/bioc_scripts/install_bioc_sysdeps.sh \
    && bash ./install_bioc_sysdeps.sh 3.16 \
    && rm ./install_bioc_sysdeps.sh

# NOTE: if anything breaks the dockerhub build cache, you will probably need to build locally and push to dockerhub.
# After the cache is in place, builds from github commits should be fast.
# NOTE: locales / locales-all added due to errors with install_deps() and special characters in the DESCRIPTION file for niaid/dsb
# NOTE: the conga rhesus branch should eventually merge, so we'll need to remove -b rhesus. The cd commands downstream are necessary to compile the reimplementation of tcrdist within conga.
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y \
		libhdf5-dev \
		libpython3-dev \
		python3-pip \
        locales \
        locales-all \
    && python3 -m pip install --upgrade pip \
    && pip3 install umap-learn phate scanpy[leiden] \
    && mkdir /conga \
    && cd /conga \
    && git clone -b rhesus https://github.com/phbradley/conga.git \
    && cd conga/tcrdist_cpp \
    && make \
    && cd ../ \
    && pip3 install -e . \
    && cd / \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# NOTE: for some reason 'pip3 install git+https://github.com/broadinstitute/CellBender.git' doesnt work.
# See: https://github.com/broadinstitute/CellBender/issues/93
RUN cd / \
    && git clone https://github.com/broadinstitute/CellBender.git CellBender \
    && chmod -R 777 /CellBender \
    && pip3 install -e CellBender

# NOTE: ggplot2 added to force version 3.4.0, which is needed by ggtree. Otherwise this container is pegged to ./focal/2022-10-28
ENV CRAN='https://packagemanager.rstudio.com/cran/__linux__/focal/latest'

# Let this run for the purpose of installing/caching dependencies
RUN echo "local({r <- getOption('repos') ;r['CRAN'] = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest';options(repos = r);rm(r)})" >> ~/.Rprofile \
    && Rscript -e "install.packages(c('remotes', 'devtools', 'BiocManager'), dependencies=TRUE, ask = FALSE, upgrade = 'always')" \
	&& echo "local({options(repos = BiocManager::repositories('https://packagemanager.rstudio.com/cran/__linux__/focal/latest'))})" >> ~/.Rprofile \
	&& echo "Sys.setenv(R_BIOC_VERSION=as.character(BiocManager::version()));" >> ~/.Rprofile \
	# NOTE: this was added to avoid the build dying if this downloads a binary built on a later R version
	&& echo "Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true');" >> ~/.Rprofile \
    && Rscript -e "print(version)" \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && mkdir /BiocFileCache && chmod 777 /BiocFileCache \
    # See: https://jmarchini.org/software/
    && wget -O /bin/sda_static_linux https://www.dropbox.com/sh/chek4jkr28qnbrj/AADPy1qQlm3jsHPmPdNsjSx2a/bin/sda_static_linux?dl=1 \
    && chmod +x /bin/sda_static_linux

ENV RETICULATE_PYTHON=/usr/bin/python3

# Create location for BioConductor AnnotationHub/ExperimentHub caches:
ENV ANNOTATION_HUB_CACHE=/BiocFileCache
ENV EXPERIMENT_HUB_CACHE=/BiocFileCache

# This should not be cached if the files change
ADD . /CellMembrane

RUN cd /CellMembrane \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT'; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = FALSE);" \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    # NOTE: Related to: https://github.com/satijalab/seurat/issues/7328. Should revert to a release once patched.
    && Rscript -e "remotes::install_github('satijalab/seurat', ref='443ab86684253d9a7290c3d38c2bc1d8db021776');" \
    && R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& unset GITHUB_PAT