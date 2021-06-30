FROM continuumio/miniconda3:4.9.2-alpine


WORKDIR /opt/app/

COPY . ./

RUN conda env create -f environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "-n", "ADVenv", "/bin/ash", "-c"]

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .





