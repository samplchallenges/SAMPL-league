#FROM mcs07/rdkit:latest
FROM ghcr.io/megosato/basepyrd:latest

WORKDIR /opt/app

COPY coords.py ./

ENV PATH="/root/.local/bin:$PATH"

#RUN python3 -m pip install --user openbabel
#RUN python3 -m pip install --user .

ENTRYPOINT ["python3", "/opt/app/coords.py"]
