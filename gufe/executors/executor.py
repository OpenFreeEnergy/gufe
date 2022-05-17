import datetime
from typing import Iterable, Union

from .results import ResultStoreClient


class Executor:
    """Destination for `AlchemicalNetwork` submission and result requests via `Client`.
    
    The `Executor` exposes a RESTful API for incoming `Client` requests.
    It also coordinates all activity between `ResultStore` and connected `Scheduler`s.

    Submitted `AlchemicalNetwork`s and `Strategy` are applied to produce
    `ProtocolDAG`s, which are made available for retrieval by connected `Scheduler`s.
    As `Scheduler`s complete `ProtocolDAG`s, they return `ProtocolResult`s to the `Executor`.
    The `Executor` deposits these in the `ResultStore` for long-term persistence.

    Attributes
    ----------
    results : ResultStoreClient
        Client interface to storage resource for all depositing and retrieving
        all `ProtocolResult`s for submitted `AlchemicalNetwork`s.
    state : str
        URI for state store; this is a log-oriented DB that can be used to
        recapitulate current `Executor` state, or any past state.
    config : str
        URI for Executor configuration.

    """
    def __init__(self,
            results: ResultStoreClient,
            state: str,
            ):
        ...


    def export_state(self, timestamp: datetime.datetime):
        """Export complete state as it was at given `timestamp`.

        This method is useful for checkpointing the `Executor` state for off-site backup.

        """
        ...


class LocalExecutor(Executor):
    """In-process, synchronous implementation of an `Executor`.

    Useful for testing and single-machine execution.

    May also work for multi-machine execution on a network filesystem, but may
    fail depending on the filesystem due to insufficent locking.

    """
    ...
