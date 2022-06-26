from pathlib import Path
from unittest.mock import patch

import pytest

from .. import filecache, models

TEST_DATA_PATH = Path(__file__).parent / "data"


@pytest.mark.parametrize(["remote_storage"], [[False], [True]])
def test_delete_local_copy(challenge, benzene_from_mol, settings, remote_storage):

    filename = "Conformer3D_CID_241.mdl"
    output_path = TEST_DATA_PATH / filename
    file_value = models.FileValue.from_string(
        output_path, challenge=challenge, input_element=benzene_from_mol
    )
    file_value.save()
    field_file = file_value.value
    with patch("os.remove") as mock_remove:
        if remote_storage:
            # Note that we're not fully configuring S3 storage for the test
            # So we have to save the FileValue before overriding settings
            settings.DEFAULT_FILE_STORAGE = "storages.backends.s3boto3.S3Boto3Storage"
        filecache.delete_local_cache(field_file)
        assert mock_remove.called == remote_storage
