# Vendored from https://github.com/openforcefield/openff-interchange/blob/main/openff/interchange/components/_packmol.py
import numpy
from openff.toolkit import Quantity


def _box_vectors_are_in_reduced_form(box_vectors: Quantity) -> bool:
    """
    Return ``True`` if the box is in OpenMM reduced form; ``False`` otherwise.

    These conditions are shared by OpenMM and GROMACS and greatly simplify
    working with triclinic boxes. Any periodic system can be represented in this
    form by rotating the system and lattice reduction.
    See http://docs.openmm.org/latest/userguide/theory/05_other_features.html#periodic-boundary-conditions
    """
    assert box_vectors.shape == (3, 3)
    a, b, c = box_vectors.m
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    return (
        [ay, az] == [0, 0]
        and bz == 0
        and ax > 0
        and by > 0
        and cz > 0
        and ax >= 2 * numpy.abs(bx)
        and ax >= 2 * numpy.abs(cx)
        and by >= 2 * numpy.abs(cy)
    )
