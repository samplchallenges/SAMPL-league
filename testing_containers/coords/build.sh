docker build -t robbason/calc-coords:latest .
docker tag robbason/calc-coords:latest ghcr.io/robbason/calc-coords:latest
docker push ghcr.io/robbason/calc-coords:latest
