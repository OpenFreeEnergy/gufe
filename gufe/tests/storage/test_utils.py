import pytest
from gufe.storage.utils import (
    REMAPPED_CLASSES, import_qualname, get_defaults_from_init
)

def test_remapped_classes():
    # While remapped classes is empty, as it originally is, this test is a
    # do-nothing. However, once it contains mappings, it will ensure that we
    # and old object name stays importable, even if the object is moved more
    # than once.
    for old_mod, old_cls in REMAPPED_CLASSES:
        try:
            import_qualname(old_mod, old_cls, REMAPPED_CLASSES)
        except Exception as e:
            raise AssertionError(
                "An error occurred while attempting to import the object "
                f"formerly know as '{old_mod}.{old_cls}'. Got: "
                f"{e.__class__}: {str(e)}"
            )

def test_import_qualname():
    pytest.skip()


def test_get_defaults_from_init():
    pytest.skip()
