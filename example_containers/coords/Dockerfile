#FROM mcs07/rdkit:latest
FROM robbason/basepyrd:latest

WORKDIR /opt/app

COPY setup.py coords.py ./

ENV PATH="/root/.local/bin:$PATH"

RUN python3 -m pip install --user openbabel
RUN python3 -m pip install --user .

ENTRYPOINT ["calc-coords"]
