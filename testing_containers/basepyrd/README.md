# How to build


Run `docker build -t basepyrd:0.1 .` to build the container.
Replace if you want to push the container to the hub, include your username in the tag (for example if your user name is mmh42) you would run `docker build -t mmh42/calc-molwt:0.1 .`

# How to use

Extend this to build other containers more quickly

# More information

See the documentation here: https://docs.docker.com/engine/reference/commandline/run/ for more information on the `docker run` command line options and see the documentation here: https://docs.docker.com/engine/reference/commandline/build/ for more information about the `docker build` command.
