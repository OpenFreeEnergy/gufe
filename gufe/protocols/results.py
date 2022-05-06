
class ProtocolUnitResult:
    """Result for a single `ProtocolUnit` execution.

    """
    ...


class ProtocolDAGResult:
    """Result for a single `ProtocolDAG` execution.

    There may be many of these in a given `ResultStore` for a given `Transformation.protocol`.

    """
    ...

class ProtocolResult:
    """Container for all `ProtocolDAGResult`s for a given `Transformation.protocol`.

    There will be exactly one of these in a given `ResultStore` for a given `Transformation.protocol`.

    """
    ...
