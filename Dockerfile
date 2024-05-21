FROM ghcr.io/bimberlabinternal/discvr-base:latest

ARG GH_PAT='NOT_SET'

# This should not be cached if the files change
ADD . /CellMembrane

RUN cd /CellMembrane \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT'; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = FALSE);" \
    # Pre-install the cpu-only versions of pytorch to try to reduce image size:
    && python3 -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    # TODO: this is to fix the as_cholmod_sparse' not provided by package 'Matrix' errors. This should ultimately be removed
    && Rscript -e "install.packages('Matrix', type = 'source', force = TRUE, repos = 'https://cran.wustl.edu/')" \
    && Rscript -e "install.packages('irlba', type = 'source', force = TRUE, repos = 'https://cran.wustl.edu/')" \
    && R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& unset GITHUB_PAT
