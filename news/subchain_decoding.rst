**Added:**

* ``KeyedChain.decode_subchains(func)`` allows decoding constituent ``GufeTokenizable`` objects whenever the passed function evaluates to a truthy value when operating on the keyed dicts of those gufe tokenizables within the ``KeyedChain``.
* ``KeyedChain.to_gufe`` now optionally accepts a dict-like object used for ``GufeTokenizable`` caching when decoding the ``KeyedChain``. This is useful when decoding multiple constituent ``GufeTokenizable`` objects from the same ``KeyedChain``. It is recommended to use ``KeyedChain.decode_subchains``, which automatically takes advantage of this feature, instead of this mechanism directly.

**Changed:**

* <news item>

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
