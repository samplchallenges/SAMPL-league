import os.path
import tempfile

import pytest

import ever_given.wrapper

@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_run_inputfile_only(container_engine):
    test_mdlfile_rel = "data/ChEBI_16716.mdl"
    print(os.path.dirname(__file__))
    test_mdlfile_abs = os.path.join(os.path.dirname(__file__), test_mdlfile_rel)
    container_uri = "ghcr.io/megosato/calc-molwt:latest"
    file_kwargs = {"molfile": test_mdlfile_abs}
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
    
    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(78.046950192)


@pytest.mark.parametrize(["container_engine"], [["docker"], ["singularity"]])
def test_run_outputfile_only(tmp_path, container_engine):
    container_uri = "ghcr.io/megosato/calc-molwt-outfile:latest"
    kwargs = {"smiles": "c1cccnc1"}
    results = {
        key: value
        for key, value in ever_given.wrapper.run(
            container_uri, 
            kwargs=kwargs,
            output_file_keys=["outfile"],
            output_dir=tmp_path,
            container_type="docker",
            engine_name=container_engine,
            file_kwargs={}
        )
    }

    assert set(results.keys()) == {"numAtoms", "numBonds", "molWeight", "outfile"}
    molWeight = float(results["molWeight"].strip())
    assert molWeight == pytest.approx(79.04219916)
    assert os.path.exists(os.path.join(tmp_path, "outfile.txt"))

    with open(os.path.join(tmp_path, "outfile.txt"), 'r') as of:
        assert of.readline().startswith("Keys: molWeight, numAtoms, numBonds")
