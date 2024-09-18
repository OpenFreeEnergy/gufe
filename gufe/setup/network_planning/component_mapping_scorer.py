# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
from gufe.tokenization import GufeTokenizable

from  .component_mapping import ComponentMapping

class ComponentMappingScorer(GufeTokenizable):
    """A generic class for scoring Atom mappings.
    this class can be used for example to build graph algorithm based networks.
    Implementations of this class can require an arbitrary and non-standardised
    number of input arguments to create.
    Implementations of this class provide the :meth:`.get_score` method
    """

    def __call__(self, mapping: ComponentMapping) -> float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self, mapping: ComponentMapping) -> float:
        """ calculate the score for an  :class:`.AtomMapping`
            the scoring function returns a value between 0 and 1.
            a value close to 1.0 indicates a small change, a score close to zero indicates a large cost/change.
        Parameters
        ----------
        mapping: AtomMapping
            the mapping to be scored
        args
        kwargs
        Returns
        -------
        float
            a value between [0,1] where zero is a very bad score and one a very good one.
        """
        raise NotImplementedError("This function was not implemented.")
