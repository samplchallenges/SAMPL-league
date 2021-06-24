FROM continuumio/miniconda3:4.9.2-alpine
WORKDIR /opt/app

COPY oe_license.txt setup.py print_logp.py environment.yml ./

RUN conda env create -f environment.yml && \
    conda clean --all --yes
RUN conda install -n SAMPL-hello-world -c conda-forge rdkit


SHELL ["conda", "run", "-n", "SAMPL-hello-world", "/bin/ash", "-c"]

RUN pip install .


ENTRYPOINT ["conda", "run", "-n", "SAMPL-hello-world", "print-LogP"]
