FROM taccsciapps/autodock-vina

ADD vina ~
ADD helloworld.sh ./

RUN sh ./helloworld.sh

RUN /usr/bin/time -p vina

