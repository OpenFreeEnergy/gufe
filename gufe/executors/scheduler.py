from typing import Iterable

from .client import Client


class Scheduler:
    """A process that polls one or more `Executor`s for `ProtocolDAG`s to
    execute on its resources.

    May include its own implementations of certain `ProtocolUnit`s, specific to
    its compute platform requirements.

    """

    def __init__(self, clients: Iterable[Client]):
        ...
