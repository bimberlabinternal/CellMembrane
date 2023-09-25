# Note: this is the last base version supporting ubuntu focal, not jammy
FROM rocker/rstudio:4.2.1

ARG GH_PAT='NOT_SET'

## Redo the R installation, since we need a base image using focal, but updated R version:
# This should be removed in favor of choosing a better base image once Exacloud supports jammy
ENV R_VERSION=4.3.1
ENV CRAN=https://packagemanager.posit.co/cran/__linux__/focal/latest
RUN /bin/sh -c /rocker_scripts/install_R_source.sh \
  && /bin/sh -c /rocker_scripts/setup_R.sh

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
        wget \
        git \
    && python3 -m pip install --upgrade pip \
    && pip3 install umap-learn phate scanpy[leiden] \
    && pip3 install git+https://github.com/broadinstitute/CellBender.git \
    && python3 -m pip install --user git+https://github.com/kmayerb/tcrdist3.git@0.2.2 \
    && mkdir /conga \
    && cd /conga \
    && git clone https://github.com/phbradley/conga.git \
    && cd conga/tcrdist_cpp \
    && make \
    && cd ../ \
    && pip3 install -e . \
    ##  Add Bioconductor system dependencies
    && mkdir /bioconductor && cd /bioconductor \
    && wget -O install_bioc_sysdeps.sh https://raw.githubusercontent.com/Bioconductor/bioconductor_docker/master/bioc_scripts/install_bioc_sysdeps.sh \
    && bash ./install_bioc_sysdeps.sh 3.17 \
    && cd / \
    && rm -Rf /bioconductor \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 cache purge \
    # This is to avoid the numba 'cannot cache function' error, such as: https://github.com/numba/numba/issues/5566
    && mkdir /numba_cache && chmod -R 777 /numba_cache \
    && mkdir /mpl_cache && chmod -R 777 /mpl_cache

ENV NUMBA_CACHE_DIR=/numba_cache
ENV MPLCONFIGDIR=/mpl_cache

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
