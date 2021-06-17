FROM mcs07/rdkit:latest

WORKDIR /opt/app

COPY setup.py molwt.py ./

RUN apt update && apt install -y python3-pip && apt clean

ENV PATH="/root/.local/bin:$PATH"

RUN python3 -m pip install --user .

ENTRYPOINT ["calc-mw"]
