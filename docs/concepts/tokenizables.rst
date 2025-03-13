
GufeTokenizables
================

Nearly all objects in gufe are subclasses of :class:`.GufeTokenizable`.
This base class enforces common behavior and sets requirements necessary
to guarantee performance and reproducibility between all downstream packages that use gufe.

By definition, a ``GufeTokenizable`` must be:

1. immutable
2. hashable
3. serializable




.. toctree::

    immutability
    gufe_keys
    serialization