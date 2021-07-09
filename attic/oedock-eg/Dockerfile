FROM oedock-base

WORKDIR /opt/app

COPY . ./


RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "oepy37", "oedock"]
