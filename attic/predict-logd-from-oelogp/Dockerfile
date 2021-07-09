FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY oe_license.txt setup.py print_logD.py environment.yml ./

RUN conda env update -f environment.yml && \
    conda clean --all --yes

RUN pip install .
