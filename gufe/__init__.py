# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from importlib.metadata import version

from . import tokenization, visualization
from .chemicalsystem import ChemicalSystem
from .components import Component, ProteinComponent, SmallMoleculeComponent, SolventComponent
from .ligandnetwork import LigandNetwork
from .mapping import (
    AtomMapper,  # more specific to atom based components
    AtomMapping,
    ComponentMapping,  # how individual Components relate
    LigandAtomMapping,
)
from .network import AlchemicalNetwork
from .protocols import (  # potentially many DAGs together, giving an estimate
    Context,
    Protocol,  # description of a method
    ProtocolDAG,  # many Units forming a workflow
    ProtocolDAGResult,  # the collected result of a DAG
    ProtocolResult,
    ProtocolUnit,  # the individual step within a method
    ProtocolUnitResult,  # the result of a single Unit
)
from .settings import Settings
from .transformations import NonTransformation, Transformation

__version__ = version("gufe")
