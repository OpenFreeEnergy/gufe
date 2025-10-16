**Added:**

* <news item>

**Changed:**

* ``FloatQuantity`` is no longer supported. Instead, use `GufeQuantity` and `specify_quantity_units()` to make a `TypeAlias`.
* System generator setting ``nonbonded_cutoff`` no longer attempts to coerce ambiguous inputs to ``unit.nanometer``. Instead, a length unit is required, e.g. ``2.2 * unit.nanometer`` or ``"2.2 nm"``.
* ``ThermoSettings`` parameters ``pressure`` and ``temperature`` no longer attempt to coerce ambiguous inputs to unts. Instead, the units must be passed explicitly, e.g. ``1.0 * units.bar`` or ``"1 bar"`` for pressure, and ``300 * unit.kelvin`` or ``"300 kelvin"`` for temperature.

**Deprecated:**

*

.. TODO: add a link to docs
**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
