FROM mambaorg/micromamba:1.5.0

WORKDIR /workspace

COPY . /workspace

RUN micromamba install -y -n base -c conda-forge \
    python=3.11 \
    pip \
    numpy \
    biopython \
    requests \
    requests-toolbelt \
    pymol-open-source \
    mmseqs2 \
    && micromamba clean -afy

RUN pip install --no-cache-dir nextflow

ENV PYTHONUNBUFFERED=1
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

ENTRYPOINT ["/bin/bash"]
