# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe


class GufeProtocolError(Exception):
    """The base gufe error that other errors should subclass."""

# Protocol Errors
class ProtocolValidationError(GufeProtocolError):
    """Error when the protocol setup or settings can not be validated."""

class ProtocolUnitExecutionError(GufeProtocolError):
    """Error when executing a protocol unit."""

# Protocol Results Errors
class ProtocolDAGResultError(GufeProtocolError):
    """Base error when dealing with DAG results."""

class MissingUnitResultError(ProtocolDAGResultError):
    """Error when a ProtocolDAGResult has no ProtocolUnitResult(s) for a given ProtocolUnit."""

class ProtocolUnitFailureError(ProtocolDAGResultError):
    """Error when a ProtocolDAGResult contains a failed protocol unit."""
