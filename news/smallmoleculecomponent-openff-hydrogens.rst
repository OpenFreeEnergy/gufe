**Added:**

* <news item>

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Under some rare circumstances calling ``SmallMoleculeComponents.to_openff()`` may have lead to hydrogens being re-assigned when converted to OpenFF Molecules (e.g. during Protocol execution). ``SmallMoleculeComponents`` now explicitly pass the ``hydrogens_are_explicit=True`` flag on OpenFF Molecule creation to avoid this issue.

**Security:**

* <news item>
