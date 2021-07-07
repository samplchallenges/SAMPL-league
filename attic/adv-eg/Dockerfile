FROM osatom/adv-rdkit-base

WORKDIR /opt/app/

COPY . ./

ENV PATH="/root/.local/bin:$PATH"

RUN pip install .

ENTRYPOINT ["conda", "run", "-n", "ADVRDenv", "dock"]

