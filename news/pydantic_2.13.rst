**Added:**

* Added support for **Pydantic 2.13** and **Python 3.14**.

**Changed:**

*

**Deprecated:**

* Dropped support for **Python 3.11**

**Removed:**

* <news item>

**Fixed:**

* When creating a pydantic JSON serialization schema, for a ``settings`` (e.g.``OpenMMSystemGeneratorFFSettings.model_json_schema(mode="serialization")``), default settings that contain ``openff.units`` are now included.

**Security:**

* <news item>
