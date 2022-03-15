#FROM mcs07/rdkit:latest
FROM ghcr.io/megosato/basepyrd:latest

WORKDIR /opt/app

COPY coords.py ./

ENV PATH="/root/.local/bin:$PATH"


ENTRYPOINT ["python3", "/opt/app/coords.py"]
