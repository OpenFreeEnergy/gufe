import string
from pathlib import Path

import pytest

from gufe.protocols.errors import ExecutionInterrupt
from gufe.protocols.protocolunit import Context, ProtocolUnit, ProtocolUnitFailure, ProtocolUnitResult
from gufe.tests.test_tokenization import GufeTokenizableTestsMixin


class DummyUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):

        # Write mock subprocess stdout and stderr for multiple
        # "processes".  Do this before raising any exceptions to check
        # that the ProtocolUnitFailure can contain stderr or stdout.
        for (output_type, output_dir) in [("stderr", ctx.stderr), ("stdout", ctx.stdout)]:
            for process_number in range(1, 3):
                filename = output_dir / f"dummy_execute_{output_type}_process_{process_number}"
                output = f"Sample {output_type} from process {process_number}".encode()
                with open(filename, "wb") as f:
                    f.write(output)

        if an_input != 2:
            raise ValueError("`an_input` should always be 2(!!!)")

        return {"foo": "bar"}


class DummyKeyboardInterruptUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):
        if an_input != 2:
            raise KeyboardInterrupt

        return {"foo": "bar"}


class DummyExecutionInterruptUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx: Context, an_input=2, **inputs):
        if an_input != 2:
            raise ExecutionInterrupt

        return {"foo": "bar"}


class TestProtocolUnit(GufeTokenizableTestsMixin):
    cls = DummyUnit
    key = "predetermined"
    repr = "DummyUnit(qux)"

    @pytest.fixture
    def instance(self):
        return DummyUnit(name="qux")

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

            stderr = Path("stderr") / str(unit.key)
            stderr.mkdir(parents=True)

            stdout = Path("stdout") / str(unit.key)
            stdout.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)

            u: ProtocolUnitFailure = unit.execute(context=ctx, an_input=3)
            assert u.exception[0] == "ValueError"

            for output_type in ("stderr", "stdout"):
                data = u.__getattribute__(output_type)
                for process_number in range(1, 3):
                    entry = f"dummy_execute_{output_type}_process_{process_number}"
                    output = f"Sample {output_type} from process {process_number}".encode()
                    assert data[entry] == output

            # now try actually letting the error raise on execute
            with pytest.raises(ValueError, match="should always be 2"):
                unit.execute(context=ctx, raise_error=True, an_input=3)

    def test_execute_ExecutionInterrupt(self, tmpdir):
        with tmpdir.as_cwd():
            unit = DummyExecutionInterruptUnit()

            shared = Path("shared") / str(unit.key)
            shared.mkdir(parents=True)

            scratch = Path("scratch") / str(unit.key)
            scratch.mkdir(parents=True)

            stderr = Path("stderr") / str(unit.key)
            stderr.mkdir(parents=True)

            stdout = Path("stdout") / str(unit.key)
            stdout.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)

            with pytest.raises(ExecutionInterrupt):
                unit.execute(context=ctx, an_input=3)

            u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

            assert u.outputs == {"foo": "bar"}

    def test_execute_KeyboardInterrupt(self, tmpdir):
        with tmpdir.as_cwd():
            unit = DummyKeyboardInterruptUnit()

            shared = Path("shared") / str(unit.key)
            shared.mkdir(parents=True)

            scratch = Path("scratch") / str(unit.key)
            scratch.mkdir(parents=True)

            stderr = Path("stderr") / str(unit.key)
            stderr.mkdir(parents=True)

            stdout = Path("stdout") / str(unit.key)
            stdout.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)

            with pytest.raises(KeyboardInterrupt):
                unit.execute(context=ctx, an_input=3)

            u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

            assert u.outputs == {"foo": "bar"}

    def test_normalize(self, instance):
        thingy = instance.key

        assert thingy.startswith("DummyUnit-")
        assert all(t in string.hexdigits for t in thingy.partition("-")[-1])
