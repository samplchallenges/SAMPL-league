docker build -t robbason/calc-molwt:latest .
docker tag robbason/calc-molwt:latest ghcr.io/robbason/calc-molwt:latest
docker push ghcr.io/robbason/calc-molwt:latest
