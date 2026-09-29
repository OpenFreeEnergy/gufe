**Added:**

* <news item>

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Performance fix: serializing and deserializing ``ProteinComponent``\s and subclasses such as ``SolvatedPDBComponent`` and ``ProteinMembraneComponent`` is now much faster, as is converting them to OpenMM topologies; ``SmallMoleculeComponent`` serialization and mapping visualization benefit likewise. The time taken to walk a whole molecule's bonds is now linear rather than quadratic in the number of bonds, since ``Mol.GetBonds()`` is no longer used to do so. (`PR #834 <https://github.com/OpenFreeEnergy/gufe/pull/834>`_).
* Performance fix: deserializing a ``SolvatedPDBComponent`` or ``ProteinMembraneComponent`` no longer builds a throwaway intermediate ``ProteinComponent``, halving the number of whole-object serializations a deserialization costs (`PR #834 <https://github.com/OpenFreeEnergy/gufe/pull/834>`_).

**Security:**

* <news item>
