So I just yanked what we needed from from https://github.com/openforcefield/openff-models/tree/077ed7b

Some changes:

* Instead of using:
```python
try:
    from pydantic.v1 import BaseModel
except ImportError:
    from pydantic import BaseModel  # type: ignore[assignment]
```
We are going to just use from `pydantic.v1 import BaseModel` directly since we depend on the pydantic 1.x version where that was added.

Then for `types.py`, `models.py`, and `exceptions.py` I ran our formatting hooks + pyupgrade --py310-plus.
