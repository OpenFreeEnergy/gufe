# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from . import _version
__version__ = _version.get_versions()['version']

from .components import (
    Component, 
    SmallMoleculeComponent, 
    ProteinComponent, 
    SolventComponent
)

from .chemicalsystem import ChemicalSystem

from .protocols import (
    Protocol,  # description of a method
    ProtocolUnit,  # the individual step within a method
    ProtocolDAG,  # many Units forming a workflow
    ProtocolUnitResult,  # the result of a single Unit
    ProtocolDAGResult,  # the collected result of a DAG
    ProtocolResult,  # potentially many DAGs together, giving an estimate
)
from .mapping import (
    AtomMapping, AtomMapper
)
from .transformations import Transformation, NonTransformation

from .network import AlchemicalNetwork
