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
        for output_type, output_dir in [("stderr", ctx.stderr), ("stdout", ctx.stdout)]:
            if output_dir is None:
                continue
            for process_number in range(1, 3):
                filename = Path(output_dir) / f"dummy_execute_{output_type}_process_{process_number}"
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

    @pytest.mark.parametrize("capture_stderr_stdout", [False, True])
    def test_execute(self, tmp_path, capture_stderr_stdout):
        unit = DummyUnit()

        shared = Path(tmp_path / "shared") / str(unit.key)
        shared.mkdir(parents=True)

        scratch = Path(tmp_path / "scratch") / str(unit.key)
        scratch.mkdir(parents=True)

        if capture_stderr_stdout:
            stderr = Path(tmp_path / "stderr") / str(unit.key)
            stderr.mkdir(parents=True)

            stdout = Path(tmp_path / "stdout") / str(unit.key)
            stdout.mkdir(parents=True)

            ctx = Context(shared=shared, scratch=scratch, stderr=stderr, stdout=stdout)
        else:
            ctx = Context(shared=shared, scratch=scratch)

        u: ProtocolUnitFailure = unit.execute(context=ctx, an_input=3)
        assert u.exception[0] == "ValueError"

        for output_type in ("stderr", "stdout"):
            data = getattr(u, output_type)
            if not capture_stderr_stdout:
                assert data == {}
                continue
            for process_number in range(1, 3):
                entry = f"dummy_execute_{output_type}_process_{process_number}"
                output = f"Sample {output_type} from process {process_number}".encode()
                assert data[entry] == output

        # now try actually letting the error raise on execute
        with pytest.raises(ValueError, match="should always be 2"):
            unit.execute(context=ctx, raise_error=True, an_input=3)

    def test_execute_ExecutionInterrupt(self, tmp_path):
        unit = DummyExecutionInterruptUnit()

        shared = Path(tmp_path / "shared") / str(unit.key)
        shared.mkdir(parents=True)

        scratch = Path(tmp_path / "scratch") / str(unit.key)
        scratch.mkdir(parents=True)

        ctx = Context(shared=shared, scratch=scratch, stderr=None, stdout=None)

        with pytest.raises(ExecutionInterrupt):
            unit.execute(context=ctx, an_input=3)

        u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

        assert u.outputs == {"foo": "bar"}

    def test_execute_KeyboardInterrupt(self, tmp_path):
        unit = DummyKeyboardInterruptUnit()

        shared = Path(tmp_path / "shared") / str(unit.key)
        shared.mkdir(parents=True)

        scratch = Path(tmp_path / "scratch") / str(unit.key)
        scratch.mkdir(parents=True)

        ctx = Context(shared=shared, scratch=scratch, stderr=None, stdout=None)

        with pytest.raises(KeyboardInterrupt):
            unit.execute(context=ctx, an_input=3)

        u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

        assert u.outputs == {"foo": "bar"}

    def test_normalize(self, instance):
        thingy = instance.key

        assert thingy.startswith("DummyUnit-")
        assert all(t in string.hexdigits for t in thingy.partition("-")[-1])
