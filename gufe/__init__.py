# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from importlib.metadata import version

from . import tokenization, visualization
from .chemicalsystem import ChemicalSystem
from .components import Component, ProteinComponent, SmallMoleculeComponent, SolventComponent
from .ligandnetwork import LigandNetwork
from .mapping import AtomMapper  # more specific to atom based components
from .mapping import ComponentMapping  # how individual Components relate
from .mapping import AtomMapping, LigandAtomMapping
from .network import AlchemicalNetwork
from .protocols import Protocol  # description of a method
from .protocols import ProtocolDAG  # many Units forming a workflow
from .protocols import ProtocolDAGResult  # the collected result of a DAG
from .protocols import ProtocolUnit  # the individual step within a method
from .protocols import ProtocolUnitResult  # the result of a single Unit
from .protocols import Context, ProtocolResult  # potentially many DAGs together, giving an estimate
from .settings import Settings
from .transformations import NonTransformation, Transformation

__version__ = version("gufe")
