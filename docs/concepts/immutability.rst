Immutability Requirements of GufeTokenizables
---------------------------------------------

One important restriction on :class:`.GufeTokenizable` subclasses is that
they must be functionally immutable. That is, they must have no mutable
attributes that change their functionality. In most cases, this means that
the object must be strictly immutable.

When an object is immutable, that means that none of its attributes change
after initialization. So all attributes should be set when you create an
object, and never changed after that. If your object is immutable, then it
is suitable to be a :class:`.GufeTokenizable`.

There is a special case of mutability that is also allowed, which is if the
object is functionally immutable.  As an example, consider a flag to turn on
or off usage of a cache of input-output pairs for some deterministic method.
If the cache is turned on, you first try to return the value from it, and
only perform the calculation if the inputs don't have a cached output
associated. In this case, the flag is mutable, but this has no effect on the
results. Indeed, the cache itself may be implemented as a mutable attribute
of the object, but again, this would not change the results that are
returned. It would also be recommended that an attribute like a cache, which
is only used internally, should be marked private with a leading underscore.

On the other hand, a flag that changes code path in a way that might
change the results of any operation would mean that the object cannot be a
:class:`.GufeTokenizable`.
