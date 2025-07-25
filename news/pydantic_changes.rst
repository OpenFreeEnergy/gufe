**Added:**

* <news item>

**Changed:**

* system generator setting ``nonbonded_cutoff`` no longer attempts to coerce ambiguous inputs to ``unit.nanometer``. Instead, a length unit is required, e.g. ``2.2 * unit.nanometer`` or ``"2.2 nm"``.
* ``ThermoSettings`` parameter ``temperature`` no longer attempts to coerce ambiguous inputs to ``unit.kelvin``. Instead, the units must be passed explicitly, e.g. ``300 * unit.kelvin`` or ``"300 kelvin"``.
* system generator setting ``nonbonded_method`` now is case sensitive and must be one of ``"CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "LJPME", "NoCutoff", "PME"``.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
