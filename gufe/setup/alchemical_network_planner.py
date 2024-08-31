# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import abc
from typing import Iterable
from .. import AlchemicalNetwork


class AlchemicalNetworkPlanner(abc.ABC):
    """
    this abstract class defines the interface for the alchemical Network Planners.
    """

    @abc.abstractmethod
    def __call__(self, *args, **kwargs) -> AlchemicalNetwork:
        raise NotImplementedError()
