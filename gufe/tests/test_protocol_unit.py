import pytest
from gufe.protocols.protocolunit import ProtocolUnit, Context
from gufe.tests.test_tokenization import GufeTokenizableTestsMixin

class DummyUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, **inputs):
        return {"foo": "bar"}

@pytest.fixture
def dummy_unit():
    return DummyUnit(name="qux")

class TestProtocolUnit(GufeTokenizableTestsMixin):
    cls = DummyUnit
    key = "predetermined"

    @pytest.fixture
    def instance(self, dummy_unit):
        # TODO: switch this to _set_key
        dummy_unit._set_key("predetermined")
        return dummy_unit

    def test_key_differs(self):
        u1 = DummyUnit()
        u2 = DummyUnit()
        assert u1.key != u2.key

