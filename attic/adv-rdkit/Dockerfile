FROM osatom/adv-rdkit-base:latest

WORKDIR /opt/app/

COPY run_autodock.py setup.py ./

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .


RUN ls -l /opt/app

ENTRYPOINT ["run-autodock"]

