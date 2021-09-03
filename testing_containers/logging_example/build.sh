docker build -t robbason/logging-example:latest .
docker tag robbason/logging-example:latest ghcr.io/robbason/logging-example:latest
docker push ghcr.io/robbason/logging-example:latest
