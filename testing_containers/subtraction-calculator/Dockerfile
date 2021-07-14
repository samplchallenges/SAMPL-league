FROM python:3.9-alpine

WORKDIR /opt/app

COPY setup.py subtract.py ./

ENV PATH="/root/.local/bin:$PATH"

RUN python3 -m pip install --user .

ENTRYPOINT ["calc-sub"]
