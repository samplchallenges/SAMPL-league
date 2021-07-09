FROM continuumio/miniconda3:4.9.2-alpine
WORKDIR /opt/app

COPY oe_license.txt setup.py environment-2.yml ./

RUN conda env update -f environment-2.yml && \
    conda clean --all --yes

COPY print_logp.py ./

RUN pip install .

ENTRYPOINT ["print-LogP"]
