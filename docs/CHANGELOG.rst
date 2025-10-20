===============
CHANGELOG
===============

.. current developments

v1.7.0
====================

**Added:**

* Added method ``LigandAtomMapping.get_alchemical_charge_difference()``. This replaces the functionality of the now-deprecated method ``openfe.utils.ligand_utils.get_alchemical_charge_difference()`` in ``openfe`` (`PR #602 <https://github.com/OpenFreeEnergy/gufe/pull/602>`_).
* Added method ``ChemicalSystem.contains()`` to check if a ``Component`` or ``Component`` type is present in the ``ChemicalSystem`` (`PR #608 <https://github.com/OpenFreeEnergy/gufe/pull/608>`_).
* Added method ``ChemicalSystem.get_components_of_type()`` to return a list of ``Component``\s that match the ``Component`` type in the ``ChemicalSystem`` (`PR #608 <https://github.com/OpenFreeEnergy/gufe/pull/608>`_).
* Added ``LigandAtomMapping.view_3d()`` method (previously implemented as ``openfe.utils.visualization_3D.view_mapping()`` (`PR #646 <https://github.com/OpenFreeEnergy/gufe/pull/646>`_).
* Added option for protocol developers to specify paths for storing stderr and stdout  (`PR #600 <https://github.com/OpenFreeEnergy/gufe/pull/600>`_ and `PR #638 <https://github.com/OpenFreeEnergy/gufe/pull/638>`_).

**Changed:**

* The default short range cutoff ``nonbonded_cutoff`` in OpenMMSystemGeneratorFFSettings has been reduced to 0.9 nm in line with best practices. (`PR #648 <https://github.com/OpenFreeEnergy/gufe/pull/648>`_).
* The default small molecule force field version has been updated from openff-2.1.1 to openff-2.2.1 (`PR #601 <https://github.com/OpenFreeEnergy/gufe/pull/601>`_).
* ``FloatQuantity`` is no longer supported. Instead, use ``GufeQuantity`` and ``specify_quantity_units()`` to make a ``TypeAlias``. See `the how-to guide <https://gufe.openfree.energy/en/v1.7.0/how-tos/custom_quantities.html>`_ for a small example of how to define a custom ``Quantity``.(`PR #584 <https://github.com/OpenFreeEnergy/gufe/pull/584>`_).
* System generator setting ``nonbonded_cutoff`` no longer attempts to coerce ambiguous inputs to ``unit.nanometer``. Instead, a length unit is required, e.g. ``2.2 * unit.nanometer`` or ``"2.2 nm"`` (`PR #584 <https://github.com/OpenFreeEnergy/gufe/pull/584>`_).
* ``ThermoSettings`` parameters ``pressure`` and ``temperature`` no longer attempt to coerce ambiguous inputs to unts. Instead, the units must be passed explicitly, e.g. ``1.0 * units.bar`` or ``"1 bar"`` for pressure, and ``300 * unit.kelvin`` or ``"300 kelvin"`` for temperature (`PR #584 <https://github.com/OpenFreeEnergy/gufe/pull/584>`_ and `PR #637 <https://github.com/OpenFreeEnergy/gufe/pull/637>`_).


**Fixed:**

* We now ensure that all ``ProtocolUnit``\s in ``protocol_units`` for a ``ProtocolDAGResult`` are included in ``self._unit_result_mapping``, ensuring that calls to ``self.ok()`` do not raise a ``KeyError`` in cases where a ``ProtocolUnitResult`` was never executed for a given ``ProtocolUnit`` (`PR #622 <https://github.com/OpenFreeEnergy/gufe/pull/622>`_).
* Fixed bug where pdb files containing phosphorous would cause an error when creating a ``ProteinComponent`` (`PR #639 <https://github.com/OpenFreeEnergy/gufe/pull/639>`_)



v1.6.1
====================

**Changed:**

* ``ambertools`` is no longer a dependency (`PR #620 <https://github.com/OpenFreeEnergy/gufe/pull/620>_`)



v1.6.0
====================

**Added:**

* Added support for python 3.13, including support for serialization/deserialization between python 3.12 and 3.13 (`PR #577 <https://github.com/OpenFreeEnergy/gufe/pull/577>`_).

**Changed:**

* Consolidated the contents of ``gufe.custom_codecs`` and ``gufe.custom_json`` into the ``gufe.serialization.json`` module (`PR #532 <https://github.com/OpenFreeEnergy/gufe/pull/532>`_).



v1.5.0
====================

**Added:**

* Added support OpenMM 8.2 (PR `#539 <https://github.com/OpenFreeEnergy/gufe/pull/539>`_)



v1.4.1
====================

**Fixed:**

* Fixed a typo in the gufe 2D visualization code which affected bond highlighting (`PR #545 <https://github.com/OpenFreeEnergy/gufe/pull/545>`_)



v1.4.0
====================

**Added:**

* ``GufeTokenizable`` objects now have the ``to_msgpack`` and ``from_msgpack`` methods for MessagePack serialization and deserialization. `(PR #372) <https://github.com/OpenFreeEnergy/gufe/issues/372>`_
* ``Protocol`` objects now have the ``_validate`` abstract method and ``validate`` method. If ``_validate`` is implemented by a ``Protocol`` author, calls to ``validate`` will validate the inputs provided without needing to create a ``ProtocolDAG``. `(PR # 412) <https://github.com/OpenFreeEnergy/gufe/issues/412>`_

**Changed:**

* ``ProteinComponent.to_openmm_topology`` now also adds some OpenMM bond type information as objects, if detected. (#501)
* An expanded set of ions has been added to ``components/proteincomponent.py`` based on amber (https://github.com/Amber-MD/AmberClassic/blob/42e88bf9a2214ba008140280713a430f3ecd4a90/dat/leap/lib/atomic_ions.lib#L1C1-L68C6).
  ``ions_dict`` provides names and charges. Now, net charges for ions are taken from it,
  instead of ``GetDefaultValence()`` from ``rdkit.Chem.GetPeriodicTable()`` function.
* ``openff.models`` is no longer supported, so we now vendor it.
  This will not require any changes in user code.
  See `issue #397 <https://github.com/OpenFreeEnergy/openfe/issues/397>`_ and `PR #535 <https://github.com/OpenFreeEnergy/openfe/pull/535>`_ for more details.

**Fixed:**

* Fixed behavior of ``ProteinComponent.to_openmm_topology`` to correctly assign ``int`` or ``None`` types to bond order for OpenMM, instead of objects. (#501)



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
  gufe objects to be deduplicated when serializing GufeTokenizables
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
