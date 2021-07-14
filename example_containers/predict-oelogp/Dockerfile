FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY oe_license.txt setup.py environment.yml ./

RUN conda env update -f environment.yml && \
    conda clean --all --yes

COPY print_logP.py ./

RUN pip install .

ENTRYPOINT ["print-LogP"]
