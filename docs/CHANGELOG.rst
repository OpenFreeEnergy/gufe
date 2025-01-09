===============
gufe Change Log
===============

.. current developments

v1.2.0
====================

**Added:**

* Added `Protocol` errors hierarchy
* Added `AtomMappingError`
* Added LigandNetwork.trim_graph
* Added warning when pickling an ``ExplicitMoleculeComponent`` that RDKit mol properties not preserved by default.
* JSON encoder now uses `zstandard compression <https://github.com/OpenFreeEnergy/gufe/pull/438>`_ .



v1.1.0
====================

**Added:**

* Use rever to manage changelog
* KeyedChain GufeTokenizable representation was added, allowing
  GUFE objects to be deduplicated when serializing GufeTokenizables
  (PR #286).
* Added `to_json` and `from_json` convenience methods to GufeTokenizables
  to more easily convert to a JSON keyed chain representation (PR #368).

**Changed:**

* Minimum Python version has been raised to v3.10 (PR #340)

**Fixed:**

* Fixed an issue where partial charges were not being read from rdkit
  Molecules where atom level properties were not set. This occured
  mainly when reading from an SDF file with partial charge tags (PR #312).
* Fixed an issue where ProtocolDAG DAG order & keys were unstable /
  non-deterministic between processes under some circumstances (PR #315).
* Fixed a bug where edge annotations were lost when converting a ``LigandNetwork`` to graphml, all JSON codec types are now supported.
