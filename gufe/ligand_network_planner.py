import abc
from typing import Iterable

from . import SmallMoleculeComponent, LigandNetwork
from . import AtomMapper
from .mapping.atom_mapping_scorer import AtomMappingScorer

class NetworkPlanner(abc.ABC):
    """A generic class for calculating :class:`.LigandNetworks`.

    Implementations of this class can require an arbitrary and non-standardised
    number of input arguments to create.

    Implementations of this class provide the :meth:`.get_score` method

    """

    def __init__(self, mapper: AtomMapper, scorer:AtomMappingScorer, *args, *kwargs):
        """ Generate a Ligand Network Planner. This class in general needs a mapper and a scorer.

        Parameters
        ----------
        mapper: AtomMapper
        scorer: AtomMappingScorer
        args
        kwargs
        """
        self.mapper = mapper
        self.scorer =  scorer


    def __call__(self, ligands: Iterable[SmallMoleculeComponent])-> LigandNetwork:
        return self.generate_ligand_network(*args, **kwargs)

    @abc.abstractmethod
    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent])->LigandNetwork:
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        ligands : Iterable[SmallMoleculeComponent]
        the ligands to include in the Network

        Returns
        -------
        LigandNetwork
            A Network, that connects all ligands with each other.
        """
