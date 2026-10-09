import string
from pathlib import Path

import pytest

from gufe.protocols.errors import ExecutionInterrupt
from gufe.protocols.protocolunit import Context, ProtocolUnit, ProtocolUnitFailure, ProtocolUnitResult
from gufe.storage.externalresource import MemoryStorage
from gufe.tests.test_tokenization import GufeTokenizableTestsMixin


@pytest.fixture
def scratch_storage(tmpdir):
    """Fixture to provide a scratch directory for ProtocolUnit tests."""
    with tmpdir.as_cwd():
        scratch = Path("scratch")
        scratch.mkdir(parents=True)
        yield scratch


@pytest.fixture
def shared_storage():
    yield MemoryStorage()


@pytest.fixture
def permanent_storage():
    yield MemoryStorage()


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
    def test_execute(self, tmpdir, scratch_storage, shared_storage, permanent_storage, capture_stderr_stdout):
        unit = DummyUnit()

        if capture_stderr_stdout:
            stderr = Path("stderr") / str(unit.key)
            stderr.mkdir(parents=True)

            stdout = Path("stdout") / str(unit.key)
            stdout.mkdir(parents=True)

            ctx = Context(
                scratch=scratch_storage,
                dag_label="test",
                unit_label=unit.key,
                stderr=stderr,
                stdout=stdout,
                shared_storage=shared_storage,
                permanent_storage=permanent_storage,
            )
        else:
            ctx = Context(
                scratch=scratch_storage,
                dag_label="test",
                unit_label=unit.key,
                shared_storage=shared_storage,
                permanent_storage=permanent_storage,
            )

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

    def test_execute_ExecutionInterrupt(self, scratch_storage, shared_storage, permanent_storage):
        unit = DummyExecutionInterruptUnit()

        shared = Path("shared") / str(unit.key)
        shared.mkdir(parents=True)

        scratch = Path("scratch") / str(unit.key)
        scratch.mkdir(parents=True)

        ctx = Context(
            shared_storage=shared_storage,
            permanent_storage=permanent_storage,
            dag_label="test",
            unit_label=unit.key,
            scratch=scratch_storage,
        )

        with pytest.raises(ExecutionInterrupt):
            unit.execute(context=ctx, an_input=3)

        u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

        assert u.outputs == {"foo": "bar"}

    def test_execute_KeyboardInterrupt(self, scratch_storage, permanent_storage, shared_storage):
        unit = DummyKeyboardInterruptUnit()

        shared = Path("shared") / str(unit.key)
        shared.mkdir(parents=True)

        scratch = Path("scratch") / str(unit.key)
        scratch.mkdir(parents=True)

        ctx = Context(
            shared_storage=shared_storage,
            permanent_storage=permanent_storage,
            dag_label="test",
            unit_label=unit.key,
            scratch=scratch_storage,
        )

        with pytest.raises(KeyboardInterrupt):
            unit.execute(context=ctx, an_input=3)

        u: ProtocolUnitResult = unit.execute(context=ctx, an_input=2)

        assert u.outputs == {"foo": "bar"}

    def test_normalize(self, instance):
        thingy = instance.key

        assert thingy.startswith("DummyUnit-")
        assert all(t in string.hexdigits for t in thingy.partition("-")[-1])


class TestContext:
    """Test the Context class context manager functionality."""

    def test_context_manager_enter_exit(
        self, scratch_storage, shared_storage: MemoryStorage, permanent_storage: MemoryStorage
    ):
        """Test that Context can be used as a context manager."""
        ctx = Context(
            dag_label="test",
            unit_label="test_unit",
            scratch=scratch_storage,
            shared_storage=shared_storage,
            permanent_storage=permanent_storage,
        )
        file_text = b"Hello World!"

        # Test __enter__
        with ctx as context:
            assert context is ctx
            assert ctx.shared.scratch_dir == scratch_storage
            assert ctx.permanent.scratch_dir == scratch_storage
            filename = "test.txt"
            test_file = context.scratch / filename
            context.shared.register(filename)
            context.permanent.register(filename)
            with test_file.open("b+w") as f:
                f.write(file_text)

        # Test __exit__ - should transfer the file to a namespaced location in shared_storage
        assert shared_storage.exists("test/test_unit/test.txt")
        with shared_storage.load_stream("test/test_unit/test.txt") as item:
            out = item.read()
            assert out == file_text

        assert permanent_storage.exists("test/test_unit/test.txt")
        with permanent_storage.load_stream("test/test_unit/test.txt") as item:
            out = item.read()
            assert out == file_text

    def test_context_manager_cleanup_stdout_stderr(
        self, scratch_storage, shared_storage: MemoryStorage, permanent_storage: MemoryStorage, tmp_path
    ):
        stdout: Path = tmp_path / "stdout"
        stdout.mkdir()
        # We write something into stdout
        stderr: Path = tmp_path / "stderr"
        stderr.mkdir()

        ctx = Context(
            dag_label="test",
            unit_label="test_unit",
            scratch=scratch_storage,
            shared_storage=shared_storage,
            permanent_storage=permanent_storage,
            stdout=stdout,
            stderr=stderr,
        )
        # Validate the directory is empty
        assert any(Path(stdout).iterdir()) == False
        assert any(Path(stderr).iterdir()) == False
        file_text = b"Hello world"
        with ctx as context:
            filename = "test.txt"
            stdout_file = context.stdout / filename
            stderr_file = context.stderr / filename
            with stdout_file.open("b+w") as f:
                f.write(file_text)
            with stderr_file.open("b+w") as f:
                f.write(file_text)
            assert any(Path(stdout).iterdir()) == True
            assert any(Path(stderr).iterdir()) == True
        # Validate we cleanup
        assert not Path(stderr).exists()
        assert not Path(stdout).exists()
