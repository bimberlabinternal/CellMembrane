FROM ghcr.io/bimberlabinternal/discvr-base:latest

ARG GH_PAT='NOT_SET'

# This should not be cached if the files change
ADD . /CellMembrane

RUN cd /CellMembrane \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT'; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = FALSE);" \
    # Pre-install the cpu-only versions of pytorch to try to reduce image size:
    && python3 -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never');" \
    # This should ultimately be removed. It is being used to fix the merge() bug in SeuratObject 5.3.0
    && Rscript -e "remotes::install_version('SeuratObject', version = '5.2.0')" \
    # NOTE: this is added to try to address the "Error in `req_perform1(req, req_prep, path = path, handle = handle, resend_count = n)`: unused arguments (req_prep, resend_count = n)"
    && Rscript -e 'remotes::install_version("httr2", version = "1.2.1")' \
    && R CMD build . \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& unset GITHUB_PAT
