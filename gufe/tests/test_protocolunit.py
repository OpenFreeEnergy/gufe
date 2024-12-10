import string
from pathlib import Path

import pytest

from gufe.protocols.protocolunit import Context, ProtocolUnit, ProtocolUnitFailure, ProtocolUnitResult
from gufe.tests.test_tokenization import GufeTokenizableTestsMixin


class DummyUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):

        if an_input != 2:
            raise ValueError("`an_input` should always be 2(!!!)")

        return {"foo": "bar"}


class DummyKeyboardInterruptUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):

        if an_input != 2:
            raise KeyboardInterrupt

        return {"foo": "bar"}


@pytest.fixture
def dummy_unit():
    return DummyUnit(name="qux")


class TestProtocolUnit(GufeTokenizableTestsMixin):
    cls = DummyUnit
    key = "predetermined"
    repr = "DummyUnit(qux)"

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

            unit = DummyUnit()

            shared = Path("shared") / str(unit.key)
            shared.mkdir(parents=True)

            scratch = Path("scratch") / str(unit.key)
            scratch.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch)

            u: ProtocolUnitFailure = unit.execute(context=ctx, an_input=3)
            assert u.exception[0] == "ValueError"

            unit = DummyUnit()

            shared = Path("shared") / str(unit.key)
            shared.mkdir(parents=True)

            scratch = Path("scratch") / str(unit.key)
            scratch.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch)

            # now try actually letting the error raise on execute
            with pytest.raises(ValueError, match="should always be 2"):
                unit.execute(context=ctx, raise_error=True, an_input=3)

    def test_execute_KeyboardInterrupt(self, tmpdir):
        with tmpdir.as_cwd():

            unit = DummyKeyboardInterruptUnit()

            shared = Path("shared") / str(unit.key)
            shared.mkdir(parents=True)

            scratch = Path("scratch") / str(unit.key)
            scratch.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch)

            with pytest.raises(KeyboardInterrupt):
                unit.execute(context=ctx, an_input=3)

            u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

            assert u.outputs == {"foo": "bar"}

    def test_normalize(self, dummy_unit):
        thingy = dummy_unit.key

        assert thingy.startswith("DummyUnit-")
        assert all(t in string.hexdigits for t in thingy.partition("-")[-1])
