# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from importlib.metadata import version

from . import tokenization

from . import visualization

from .components import (
    Component,
    SmallMoleculeComponent,
    ProteinComponent,
    SolventComponent
)

from .chemicalsystem import ChemicalSystem

from .mapping import (
    ComponentMapping,  # how individual Components relate
    AtomMapping, AtomMapper,  # more specific to atom based components
    LigandAtomMapping,
)

from .settings import Settings

from .protocols import (
    Context,
    Protocol,  # description of a method
    ProtocolUnit,  # the individual step within a method
    ProtocolDAG,  # many Units forming a workflow
    ProtocolUnitResult,  # the result of a single Unit
    ProtocolDAGResult,  # the collected result of a DAG
    ProtocolResult,  # potentially many DAGs together, giving an estimate
)

from .transformations import Transformation, NonTransformation

from .network import AlchemicalNetwork
from .ligandnetwork import LigandNetwork

__version__ = version("gufe")
