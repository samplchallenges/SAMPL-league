import pytest

import ever_given.wrapper


def test_fail(container_engine):
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    file_kwargs = {"molfile": "nonexistent_file"}
    with pytest.raises(ever_given.wrapper.ContainerFailedException):
        results = {
            key: value
            for key, value in ever_given.wrapper.run(
                container_uri, 
                kwargs={},
                container_type="docker",
                engine_name=container_engine,
                file_kwargs=file_kwargs
            )
        }
        print(results)