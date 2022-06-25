# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

class ExternalResourceError(Exception):
    """Base class for errors due to problems with external resources"""
    # TODO: is it necessary to have a base class here? Would you ever have
    # one catch that handles both subclass errors?


class MissingExternalResourceError(ExternalResourceError):
    """Error when the external resource could not be loaded"""


class ChangedExternalResourceError(ExternalResourceError):
    """Error when there's a metadata mismatch with an external resource"""
