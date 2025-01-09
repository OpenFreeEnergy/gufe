from .explicitmoleculecomponent import ExplicitMoleculeComponent



class ExplicitRepeatedComponent(ExplicitMoleculeComponent):


    def __init__(self):
        """An explicit molecule that is repeated many times, with coordinates
        specified for each instance.

        """
        ...


class ExplicitSolventCompoent(ExplicitRepeatedComponent):
    ...


class MembraneComponent(ExplicitRepeatedComponent):
    ...


class IonsComponent(ExplicitRepeatedComponent):
    ...
