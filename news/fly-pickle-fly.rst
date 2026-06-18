**Added:**

* Unpickling gufe tokenizables now supports the gufe token registry (PR #797 <https://github.com/OpenFreeEnergy/gufe/pull/797>_).

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Fixed pickling of ``SmallMoleculeComponent``/``ExplicitMoleculeComponent`` objects so RDKit molecule properties, such as ``_Name``, are preserved across serialization round-trips (PR #797 <https://github.com/OpenFreeEnergy/gufe/pull/797>_).
  This will fix issues from using multiprocessing with these objects.


**Security:**

* <news item>
