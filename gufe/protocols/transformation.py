
from typing import List

class Transformation:
    """Spawns units of work, gives status of units.

    """

    def __init__(
            self,
            initial: ChemicalSystem,
            final: ChemicalSystem,
            protocol: Protocol,
            mapping: AtomMapping = None,
            ):
        ...

    def prepare(self) -> List[WorkUnit]:
        """

        """
        self._work_units = work_units

    def execute(self):
        ...

    def estimate(self) -> Result:
        ...

    def status(self):
        ...
