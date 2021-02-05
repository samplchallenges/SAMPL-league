=====
Conda
=====

.. note::
   A familiarity of conda is assumed


Example docker file

.. code-block:: docker

   FROM continuumio/miniconda3:4.9.2-alpine
   WORKDIR /opt/app
   COPY environment.yml ./
   RUN conda env create -f environment.yml && \
       conda clean --all --yes
   ENTRYPOINT ["conda", "run", "-n", "test"]

The ``ENTRYPOINT`` directive will ensure that commands passed into the container with ``run`` will use the correct conda environment. 
