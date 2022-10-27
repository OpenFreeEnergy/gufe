# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Optional, Iterable
import io
import contextlib
import warnings

from openff.toolkit.utils.serialization import Serializable
from ..tokenization import GufeTokenizable

from ..chemicalsystem import ChemicalSystem
from ..protocols import Protocol, ProtocolDAG, ProtocolResult, ProtocolDAGResult
from ..mapping import ComponentMapping


class ensure_filelike:
    """Context manager to convert pathlike or filelike to filelike.

    This makes it so that methods can allow a range of user inputs.

    Parameters
    ----------
    fn : PathLike or FileLike
        The user input to normalize.
    mode : str or None
        The mode, if ``fn`` is pathlike. If ``fn`` is filelike, a warning
        will be emitted and the mode will be ignored.
    force_close : bool, default False
        Whether to forcibly close the stream on exit. For pathlike inputs,
        the stream will always be closed. Filelike inputs will close
        if this parameter is True.
    """
    def __init__(self, fn, mode=None, force_close=False):
        filelikes = (io.TextIOBase, io.RawIOBase, io.BufferedIOBase)
        if isinstance(fn, filelikes):
            if mode is not None:
                warnings.warn(
                    f"mode='{mode}' specified with {fn.__class__.__name__}."
                    " User-specified mode will be ignored."
                )
            self.to_open = None
            self.do_close = force_close
            self.context = fn
        else:
            if mode is None:
                mode = "r"
            self.to_open = fn
            self.do_close = True
            self.context = None

        self.fn = fn

        self.mode = mode

    def __enter__(self):
        if self.to_open is not None:
            self.context = open(self.to_open, mode=self.mode)

        return self.context

    def __exit__(self, type, value, traceback):
        if self.do_close:
            self.context.close()

# def ensure_filelike(fn, mode=None):
#     # TODO: move this somewhere reusable
#     if isinstance(fn, (io.TextIOBase, io.RawIOBase, io.BufferedIOBase)):
#         if mode is not None:
#             warnings.warn(
#                 f"mode='{mode}' specified with {fn.__class__.__name__}. "
#                 "User-specified mode will be ignored."
#             )
#         return fn
#     else:
#         if mode is None:
#             mode = 'r'  # effectively given the same implicit behavior
#         return open(fn, mode=mode)


class Transformation(GufeTokenizable):
    """An edge of an alchemical network.

    Connects two `ChemicalSystem`s, with directionality.

    Attributes
    ----------
    stateA : ChemicalSystem
        The starting `ChemicalSystem` for the transformation.
    stateB : ChemicalSystem
        The ending `ChemicalSystem` for the transformation.
    protocol : Protocol
        The protocol used to perform the transformation.
        Includes all details needed to perform required
        simulations/calculations and encodes the alchemical pathway used.
        May also correspond to an experimental result.
    mapping : Optional[Dict[str, ComponentMapping]]
        Mapping of e.g. atoms between the `stateA` and `stateB`
        `ChemicalSystem`s.
    name : Optional[str]
        Optional identifier for the transformation; set this to a unique value
        if adding multiple, otherwise identical transformations to the same
        `AlchemicalNetwork` to avoid deduplication

    """

    def __init__(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        protocol: Protocol,
        mapping: Optional[dict[str, ComponentMapping]] = None,
        name: Optional[str] = None,
    ):

        self._stateA = stateA
        self._stateB = stateB
        self._mapping = mapping
        self._name = name

        self._protocol = protocol

    def _defaults(self):
        return super()._defaults()

    def __repr__(self):
        return f"{self.__class__.__name__}(stateA={self.stateA}, "\
               f"stateB={self.stateB}, protocol={self.protocol})"

    @property
    def stateA(self):
        """The starting `ChemicalSystem` for the transformation."""
        return self._stateA

    @property
    def stateB(self):
        """The ending `ChemicalSystem` for the transformation."""
        return self._stateB

    @property
    def protocol(self):
        """The protocol for sampling the transformation to derive free energy
        differences between `stateA` and `stateB` `ChemicalSystem`s.

        """
        return self._protocol

    @property
    def mapping(self):
        """The mapping between atoms in `stateA` to `stateB`"""
        return self._mapping

    @property
    def name(self):
        """User-specified for the transformation; used as part of its hash."""
        return self._name

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._name != other.name:
            return False
        if self._protocol != other.protocol:
            return False
        if self._mapping != other.mapping:
            return False
        if self._stateB != other.stateB:
            return False
        if self._stateA != other.stateA:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                self._stateA,
                self._stateB,
                self._mapping,
                self._name,
                self._protocol,
            )
        )

    def _to_dict(self) -> dict:
        return {
            "stateA": self.stateA,
            "stateB": self.stateB,
            "protocol": self.protocol,
            "mapping": self.mapping,
            "name": self.name,
        }

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(**d)

    def create(self) -> ProtocolDAG:
        """Returns a `ProtocolDAG` executing this `Transformation.protocol`."""
        return self.protocol.create(
            stateA=self.stateA,
            stateB=self.stateB,
            mapping=self.mapping,
            name=str(self.key),
        )

    def gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> ProtocolResult:
        """Gather multiple `ProtocolDAGResult`s into a single `ProtocolResult`.

        Parameters
        ----------
        protocol_dag_results : Iterable[ProtocolDAGResult]
            The `ProtocolDAGResult`s to assemble aggregate quantities from.

        Returns
        -------
        ProtocolResult
            Aggregated results from many `ProtocolDAGResult`s from a given `Protocol`.

        """
        return self.protocol.gather(protocol_dag_results=protocol_dag_results)

    def dump(self, file):
        """Dump this Transformation to a JSON file.
        """
        with ensure_filelike(file, mode='w') as f:
            json.dump(self.to_dict(), f, cls=JSON_HANDLER.encoder,
                      sort_keys=True)

    @classmethod
    def load(cls, file):
        """Create a Transformation from a JSON file.
        """
        with ensure_filelike(file, mode='r') as f:
            dct = json.load(f, cls=JSON_HANDLER.decoder)

        return cls.from_dict(dct)



# we subclass `Transformation` here for typing simplicity
class NonTransformation(Transformation):
    """A non-alchemical edge of an alchemical network.

    A "transformation" that performs no transformation at all.
    Technically a self-loop, or an edge with the same `ChemicalSystem` at
    either end.

    Functionally used for applying a dynamics protocol to a `ChemicalSystem`
    that performs no alchemical transformation at all. This allows e.g.
    equilibrium MD to be performed on a `ChemicalSystem` as desired alongside
    alchemical protocols between it and and other `ChemicalSystem`s.

    Attributes
    ----------
    system : ChemicalSystem
    protocol : Protocol
        The protocol used to perform the dynamics. Includes all details needed
        to perform required simulations/calculations.
    name : Optional[str]
        Optional identifier for the nontransformation; set this to a unique
        value if adding multiple, otherwise identical transformations to the
        same `AlchemicalNetwork` to avoid deduplication

    """

    def __init__(
        self,
        system: ChemicalSystem,
        protocol: Protocol,
        name: Optional[str] = None,
    ):

        self._system = system
        self._name = name
        self._protocol = protocol

    @property
    def stateA(self):
        return self._system

    @property
    def stateB(self):
        return self._system

    @property
    def system(self):
        return self._system

    @property
    def protocol(self):
        """The protocol for sampling dynamics of the `ChemicalSystem`."""
        return self._protocol

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._name != other.name:
            return False
        if self._protocol != other.protocol:
            return False
        if self._system != other.system:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                self._name,
                self._protocol,
                self._system,
            )
        )

    def _to_dict(self) -> dict:
        return {
            "system": self.system,
            "protocol": self.protocol,
            "name": self.name,
        }

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(**d)

    def create(self) -> ProtocolDAG:
        """Returns a `ProtocolDAG` executing this `Transformation.protocol`."""
        return self.protocol.create(
            stateA=self.system, stateB=self.system, name=str(self.__hash__())
        )
