**Added:**

* ``gufe.utils.get_bonds``, which returns a molecule's bonds in index order in linear time (`PR #834 <https://github.com/OpenFreeEnergy/gufe/pull/834>`_).

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Walking a whole molecule's bonds is now linear rather than quadratic in the number of bonds, since ``Mol.GetBonds()`` is no longer used to do so. This affects serialization, deserialization, OpenMM topology conversion, and mapping visualization; for a solvated membrane system of ~154k atoms and ~126k bonds, ``ProteinMembraneComponent.from_dict`` drops from ~290 s to ~2.8 s and ``to_dict`` from ~79 s to ~0.9 s. This is a saving in time rather than memory: peak RSS falls ~8% for ``from_dict`` and rises ~8% for ``to_dict``. Serialized output and ``GufeKey``\s are unchanged (`PR #834 <https://github.com/OpenFreeEnergy/gufe/pull/834>`_).
* ``SolvatedPDBComponent._from_dict`` no longer builds an intermediate ``ProteinComponent``, halving the number of whole-object serializations a deserialization costs (`PR #834 <https://github.com/OpenFreeEnergy/gufe/pull/834>`_).

**Security:**

* <news item>
