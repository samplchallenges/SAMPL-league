FROM oedock-base

WORKDIR /opt/app

COPY oe_io.py docking.py ligand.py receptor.py run_oedock.py setup.py ./

RUN pip install .

ENTRYPOINT ["oedock"]
