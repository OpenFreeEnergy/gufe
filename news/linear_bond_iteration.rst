**Added:**

* ``gufe.utils.iter_bonds``, which returns a molecule's bonds in index order in linear time. ``Mol.GetBonds()`` reaches each bond through ``Mol.GetBondWithIdx``, which RDKit implements as an ``O(n_bonds)`` search because bonds are graph edges rather than an indexable container; iterating it is therefore ``O(n_bonds ** 2)``.
* ``ProteinComponent._rdkit_from_dict``, which rebuilds the ``rdkit`` molecule and name from a dict representation without constructing a ``ProteinComponent``. Subclasses whose constructors take additional arguments can use it to avoid instantiating a throwaway intermediate.

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* ``ProteinComponent``, ``SolvatedPDBComponent``, ``ProteinMembraneComponent``, and ``SmallMoleculeComponent`` no longer iterate ``Mol.GetBonds()`` when serializing, deserializing, or converting to an OpenMM topology, so those operations are now linear rather than quadratic in the number of bonds. For a solvated membrane system of ~154k atoms and ~126k bonds, ``from_dict`` drops from 321 s to 3.1 s, ``to_dict`` from 115 s to 1.0 s, and ``to_openmm_topology`` from 115 s to 1.8 s. Serialized output and ``GufeKey``\s are unchanged.
* ``SolvatedPDBComponent._from_dict`` no longer builds an intermediate ``ProteinComponent``. Constructing one computes its key, which serializes the object in full, so deserializing previously paid for two whole-object serializations rather than one.

**Security:**

* <news item>
