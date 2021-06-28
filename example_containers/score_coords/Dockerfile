FROM robbason/basepyrd:latest

WORKDIR /opt/app

COPY setup.py score.py ./

ENV PATH="/root/.local/bin:$PATH"

RUN python3 -m pip install --user .

ENTRYPOINT ["score-coords"]
