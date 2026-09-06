**Added:**

* ``StorageManager`` for managing files during unit execution. Allows for auto transfer to storage mediums following execution of a unit.

**Changed:**

* ``Context`` now uses ``StorageManager`` to handle operations in protocol units. Units must now register files in scratch to be transferred shared or permanent storage.
* ``ProtocolUnit``s now have permanent storage which can be used to outlast a running ``ProtocolDAG``.


**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
