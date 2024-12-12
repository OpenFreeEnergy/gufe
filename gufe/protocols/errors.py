# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

class BaseGufeError(Exception):
    """The base gufe error that other errors should subclass."""

# Protocol Errors
class ProtocolValidationError(BaseGufeError):
    """Error when the protocol setup or settings can not be validated."""


class ProtocolSetupError(BaseGufeError):
    """Error when executing the setup unit of the protocol."""


class ProtocolExecutionError(BaseGufeError):
    """Error when executing the production unit of the protocol."""


class ProtocolAnalysisError(BaseGufeError):
    """Error when trying to perform some analyses after the protocol has been executed."""


# Protocol Results Errors
class MissingUnitResultError(BaseGufeError):
    """Error when a ProtocolDAGResult is missing a unit result."""


class ProtocolUnitFailureError(BaseGufeError):
    """Error when a ProtocolDAGResult contains a failed protocol unit."""
