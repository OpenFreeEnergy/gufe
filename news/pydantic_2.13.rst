**Added:**

* Added support for **Pydantic 2.13** and **Python 3.14**.

**Changed:**

* For all ``gufe.settings`` models, ``model_dump()`` uses ``polymorphic_serialization``, meaning that sub-models will be included in the output.
  See `Pydantic's docs on the new polymorphic serialization behavior <https://pydantic.dev/articles/pydantic-v2-13-release#polymorphic-serialization>`_ for more info.

**Deprecated:**

* Dropped support for **Python 3.11**

**Removed:**

* <news item>

**Fixed:**

* When creating a pydantic JSON serialization schema, for a ``settings`` (e.g.``OpenMMSystemGeneratorFFSettings.model_json_schema(mode="serialization")``), default settings that contain ``openff.units`` are now included.

**Security:**

* <news item>
