#FROM mcs07/rdkit:latest
FROM ghcr.io/megosato/basepyrd:latest

WORKDIR /opt/app

COPY molwt.py ./

#RUN apt update && apt install -y python3-pip && apt clean

ENV PATH="/root/.local/bin:$PATH"

ENTRYPOINT ["python3", "/opt/app/molwt.py"]
