FROM adv-base:0.1


WORKDIR /opt/app/

COPY . ./

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .

RUN /opt/app/dependencies/mgl/install.sh

ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

RUN /opt/app/dependencies/adv/bin/vina --help

ENTRYPOINT ["conda", "run", "-n", "ADVenv", "dock"]

