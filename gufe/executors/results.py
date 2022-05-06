import abc
from typing import Iterable, Tuple

from ..transformations import Transformation


class ResultStoreClient(abc.ABC):

    def estimate(self, transformations: Iterable[Transformation]) -> Tuple[Tuple[float],Tuple[float],Tuple[float]]:
        """Get free energy estimates, uncertainties, and rates of convergence
        for the given transformations.

        Parameters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.

        Returns
        -------
        dG : Iterable[float]
            Free energy estimates for the given transformations, in order.
        ddG : Iterable[float]
            Uncertainties in dG for the given transformations, in order.
        rate_of_convergence : Iterable[float]
            Rate of convergence for dG for the given transformations, in order.
        """
        ...

    def estimate_dG(self, transformations: Iterable[Transformation]) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Parameters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.

        Returns
        -------
        dG : Iterable[float]
            Free energy estimates for the given transformations, in order.
        """
        ...

    def estimate_uncertainty(self, transformations: Iterable[Transformation]) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Parameters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.

        Returns
        -------
        ddG : Iterable[float]
            Uncertainties in dG for the given transformations, in order.
        """
        ...

    def estimate_rate_of_convergence(self, transformations: Iterable[Transformation]) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Requires `client` to be defined on this `AlchemicalNetwork`.

        Parameters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.

        Returns
        -------
        rate_of_convergence : Iterable[float]
            Rate of convergence for dG for the given transformations, in order.
        """
        ...


class LocalResultStoreClient(ResultStoreClient):
    ...
