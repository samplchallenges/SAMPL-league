FROM continuumio/miniconda3:4.9.2-alpine

WORKDIR /opt/app

COPY . ./

RUN conda env create -f environment.yml && \
    conda clean --all --yes

SHELL ["conda", "run", "-n", "oepy37", "/bin/ash", "-c"]


RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "oepy37", "oedock"]
