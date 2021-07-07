FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY oe_license.txt setup.py environment.yml ./

RUN conda env create -f environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "-n", "SAMPL-hello-world", "/bin/ash", "-c"]

COPY print_logd.py ./

RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "SAMPL-hello-world", "print-LogD"]
