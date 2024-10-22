# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import abc

from ..tokenization import GufeTokenizable
from .. import AlchemicalNetwork


class AlchemicalNetworkPlanner(GufeTokenizable, abc.ABC):
    """
    this abstract class defines the interface for the alchemical Network Planners.
    """

    @abc.abstractmethod
    def __call__(self, *args, **kwargs) -> AlchemicalNetwork:
        raise NotImplementedError()
