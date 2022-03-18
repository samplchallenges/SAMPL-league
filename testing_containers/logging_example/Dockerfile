#FROM mcs07/rdkit:latest
FROM ghcr.io/megosato/basepyrd:latest

WORKDIR /opt/app

COPY logging_example.py ./

ENV PATH="/root/.local/bin:$PATH"

ENTRYPOINT ["python3", "/opt/app/logging_example.py"]
