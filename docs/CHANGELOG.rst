===============
CHANGELOG
===============

.. current developments

v1.3.0
====================

**Added:**

* Added ExecutionInterrupt, a special exception that does not get handled as a ``ProtocolUnitFailure``.

**Changed:**

* Nodes and edges are now sorted by inchikey before being added to a networkx graph in the ``LigandNetwork.graph`` property. This replaces gufekey sorting which is not stable between versions and should result in reproducible network generation.
* The message stating that partial charges are already present in an ``ExplicitMoleculeComponent`` is now included in ``logger.info``, rather than as a warning message. This should make output significantly less noisy for some users.
* ``Protocol`` subclasses now require that the ``_settings_cls``
  attribute is set to the intended ``SettingsBaseModel``
  subclass. This attribute is validated during ``Protocol``
  instantiation.
* ``GufeTokenizable.from_json`` now falls back to loading ``dict`` representation if ``from_keyed_chain`` fails
* ``ExplicitMoleculeComponent`` now uses ``GufeTokenizable`` ``to_json`` and ``from_json`` methods via inheritance

**Deprecated:**

* ``Transformation.dump``, ``Transformation.load`` are now deprecated, use ``Transformation.to_json`` and ``Transformation.from_json`` instead

**Fixed:**

* Fixed bug where an error was only being raised if the difference between the sum of partial charges and the small molecule's net charge was a positive value. Behavior has been fixed such that negative discrepancies now raise an error as well.
* Under some rare circumstances calling ``SmallMoleculeComponents.to_openff()`` may have led to hydrogens being re-assigned when converted to OpenFF Molecules (e.g. during Protocol execution). ``SmallMoleculeComponents`` now explicitly pass the ``hydrogens_are_explicit=True`` flag on OpenFF Molecule creation to avoid this issue.



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
