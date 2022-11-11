FROM bioconductor/bioconductor_docker:latest

# NOTE: if anything breaks the dockerhub build cache, you will probably need to build locally and push to dockerhub.
# After the cache is in place, builds from github commits should be fast.
# NOTE: locales / locales-all added due to errors with install_deps() and special characters in the DESCRIPTION file for niaid/dsb \
# NOTE: libicu-dev added to avoid stringi /  libicui18n.so.66: cannot open shared object file error
# NOTE: libssl-dev libcrypto-dev added due to 'Cannot find libcrypto error'
RUN apt-get update -y \
    && apt-get install -y \
		libhdf5-dev \
		libpython3-dev \
		python3-pip \
        locales \
        locales-all \
    && python3 -m pip install --upgrade pip \
	&& pip3 install umap-learn phate \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# NOTE: for some reason 'pip3 install git+https://github.com/broadinstitute/CellBender.git' doesnt work.
# See: https://github.com/broadinstitute/CellBender/issues/93
RUN cd / \
    && git clone https://github.com/broadinstitute/CellBender.git CellBender \
    && chmod -R 777 /CellBender \
    && pip3 install -e CellBender

# Let this run for the purpose of installing/caching dependencies
RUN Rscript -e "install.packages(c('remotes', 'devtools', 'BiocManager'), dependencies=TRUE, ask = FALSE, upgrade = 'always')" \
	&& echo "local({options(repos = BiocManager::repositories())})" >> ~/.Rprofile \
	&& echo "Sys.setenv(R_BIOC_VERSION=as.character(BiocManager::version()));" >> ~/.Rprofile \
	# NOTE: this was added to avoid the build dying if this downloads a binary built on a later R version
	&& echo "Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true');" >> ~/.Rprofile \
    && Rscript -e "print(version)" \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# This should not be cached if the files change
ADD . /CellMembrane

ENV RETICULATE_PYTHON=/usr/bin/python3

# Create location for BioConductor AnnotationHub/ExperimentHub caches:
ENV ANNOTATION_HUB_CACHE=/BiocFileCache
ENV EXPERIMENT_HUB_CACHE=/BiocFileCache
RUN mkdir /BiocFileCache && chmod 777 /BiocFileCache

RUN cd /CellMembrane \
	&& Rscript -e "BiocManager::install(ask = FALSE, upgrade = 'always');" \
    && Rscript -e "install.packages('rhdf5', force = TRUE, ask = FALSE, upgrade = 'always');" \
RUN	Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \

RUN	R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds