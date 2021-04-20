FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY oe_license.txt setup.py print_logP.py environment.yml ./

RUN conda env create -f environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "-n", "SAMPL-hello-world", "/bin/ash", "-c"]

RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "SAMPL-hello-world", "print-LogP"]
