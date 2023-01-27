import string
import pytest
from pathlib import Path

from gufe.protocols.protocolunit import ProtocolUnit, Context, ProtocolUnitResult, ProtocolUnitFailure
from gufe.tests.test_tokenization import GufeTokenizableTestsMixin

class DummyUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):

        if an_input !=2:
            raise ValueError("`an_input` should always be 2(!!!)")

        return {"foo": "bar"}

@pytest.fixture
def dummy_unit():
    return DummyUnit(name="qux")

class TestProtocolUnit(GufeTokenizableTestsMixin):
    cls = DummyUnit
    key = "predetermined"

    @pytest.fixture
    def instance(self, dummy_unit):
        dummy_unit._set_key("predetermined")
        return dummy_unit

    def test_key_differs(self):
        u1 = DummyUnit()
        u2 = DummyUnit()
        assert u1.key != u2.key


    def test_execute(self, tmpdir):
        with tmpdir.as_cwd():
            u: ProtocolUnitFailure = DummyUnit().execute(shared=Path('.'), an_input=3)
            assert u.exception[0] == "ValueError"

            # now try actually letting the error raise on execute
            with pytest.raises(ValueError, match="should always be 2"):
                DummyUnit().execute(shared=Path('.'), raise_error=True, an_input=3)

    def test_normalize(self, dummy_unit):
        thingy = dummy_unit.key

        assert thingy.startswith('DummyUnit-')
        assert all(t in string.hexdigits for t in thingy.partition('-')[-1])
