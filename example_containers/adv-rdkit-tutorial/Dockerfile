FROM osatom/adv-rdkit-base:latest

WORKDIR /opt/app/

ENV PATH="/root/.local/bin:$PATH"

COPY setup.py autodock.py run_autodock.py ./

RUN pip install .

ENTRYPOINT ["run-autodock"]

