# Headless test image for the wxPython desktop app (StackImages, GetFit, ShowClusters).
FROM continuumio/miniconda3:25.3.1-1
RUN apt-get update \
 && apt-get install -y --no-install-recommends xvfb xauth git zsh \
 && rm -rf /var/lib/apt/lists/* \
 && git config --system --add safe.directory /work
COPY ci/conda-linux-64.lock /tmp/conda-linux-64.lock
RUN conda create -y -q -n BobsAstro --file /tmp/conda-linux-64.lock \
 && conda clean -afy
ENV PATH=/opt/conda/envs/BobsAstro/bin:$PATH \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONDONTWRITEBYTECODE=1
WORKDIR /work
CMD ["xvfb-run", "-a", "python", "-m", "pytest", "-v", "-ra", "-p", "no:cacheprovider", "tests/desktop", "tests/packaging"]
