FROM osatom/adv-rdkit-base:latest

WORKDIR /opt/app/

ENV PATH="/root/.local/bin:$PATH"

COPY run_autodock.py autodock.py setup.py ./

RUN pip install .

ENTRYPOINT ["run-autodock"]

