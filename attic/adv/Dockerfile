FROM adv-base


WORKDIR /opt/app/

COPY . ./

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "ADVenv", "dock"]

