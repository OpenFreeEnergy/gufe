**Added:**

* KeyedChain GufeTokenizable representation was added, allowing
  GUFE objects to be deduplicated when serializing GufeTokenizables
  (PR #286).
* Added `to_json` and `from_json` convenience methods to GufeTokenizables
  to more easily convert to a JSON keyed chain representation (PR #368).

**Changed:**

* Minimum Python version has been raised to v3.10 (PR #340)

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Fixed an issue where partial charges were not being read from rdkit
  Molecules where atom level properties were not set. This occured
  mainly when reading from an SDF file with partial charge tags (PR #312).
* Fixed an issue where ProtocolDAG DAG order & keys were unstable /
  non-deterministic between processes under some circumstances (PR #315).

**Security:**

* <news item>
