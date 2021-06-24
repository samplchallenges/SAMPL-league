FROM continuumio/miniconda3

WORKDIR /opt/app

COPY setup.py autodockvina.py ./

RUN pip install .

ENTRYPOINT ["adv"]

#RUN conda install -n vina -c bioconda autodock-vina
