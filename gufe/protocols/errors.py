# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe


# Protocol Errors
class ProtocolValidationError(Exception):
    """Error when the protocol setup or settings can not be validated."""

class ProtocolSetupError(Exception):
    """Error when executing the setup unit of the protocol."""

class ProtocolExecutionError(Exception):
    """Error when executing the production unit of the protocol."""

class ProtocolAnalysisError(Exception):
    """Error when trying to perform some analyses after the protocol has been executed."""


# Protocol Results Errors
class MissingProtocolUnitError(Exception):
    """Error when a ProtocolDAGResult is massing a protocol unit."""

class ProtocolUnitFailureError(Exception):
    """Error when a ProtocolDAGResult contains an failed protocol unit."""
