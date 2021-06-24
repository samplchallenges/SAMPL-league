FROM mcs07/rdkit:latest

WORKDIR /opt/app/

COPY . ./

RUN apt update && apt install -y python3-pip && apt clean

ENV PATH="/root/.local/bin:$PATH"

RUN python3 -m pip install --user .

RUN /opt/app/dependencies/mgl/install.sh -c 1

ENV PATH="/opt/app/dependencies/mgl/:$PATH"

RUN /opt/app/dependencies/adv/bin/vina --help

ENTRYPOINT ["dock"]
