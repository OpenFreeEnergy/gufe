**Added:**

* <news item>

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* We now ensure that all ``ProtocolUnit``\s in ``protocol_units`` for a ``ProtocolDAGResult`` are included in ``self._unit_result_mapping``, ensuring that calls to ``self.ok()`` do not raise a ``KeyError`` in cases where a ``ProtocolUnitResult`` was never executed for a given ``ProtocolUnit``.

**Security:**

* <news item>
