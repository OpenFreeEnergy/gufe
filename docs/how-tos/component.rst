.. _howto-component:

How to define a new ``Component``
=================================

The **gufe**  :ref:`Component <component>` is a :class:`.GufeTokenizable` intended to be used as an *extensible point* of the library, such that you can describe a simulated system in terms of GufeTokenizables that are compatible with **gufe** and the rest of the OpenFE ecosystem.

Step 1: Choose a Component to Extend
-------------------------------------

In many cases, you will likely want your custom ``Component`` to inherit from one of the following extensible points that themselves inherit from the ``Component`` base class - meaning that all of the following are ``Component``\s, but come with additional functionality:

* :class:`.ExplicitMoleculeComponent`
* :class:`.ProteinComponent`
* :class:`.SmallMoleculeComponent`
* :class:`.BaseSolventComponent`
* :class:`.SolventComponent`

In the rare case where you want as much custom implementation as possible, you could inherit directly from :class:`.Component` itself.
However, many OpenFE ``Protocols`` require specific Component types, in which case sub-classing from a compatible class is recommended.

Step 2: Write Tests for Expected Behavior
-----------------------------------------

As when defining any new ``GufeTokenizable``, you are encouraged to use the ``GufeTokenizableTestsMixin`` pytest fixture to ensure that your new ``Component`` works as intended.

.. note::

    A key benefit of inheriting from a fully implemented class, such as ``SmallMoleculeComponent`` in this example, is that you only need to implement *new* functionality.

For example:

.. TODO: if we expose the other testing mixins, reference ExplicitMoleculeComponentMixin here.

.. code-block:: python
    :caption: custom_component.py

    from gufe import SmallMoleculeComponent

    # this is all the code you need to get the GufeTokenizableTestsMixin to pass,
    # since SmallMoleculeComponent fully implemented
    class CustomComponent(SmallMoleculeComponent):
        pass

.. code-block:: python
    :caption: test_custom_component.py

    from gufe.tests import GufeTokenizableTestsMixin
    from .custom_component import CustomComponent
    import pytest
    from rdkit import Chem


    class TestCustomComponent(GufeTokenizableTestsMixin):
        cls = CustomComponent
        repr = "CustomComponent(name=ethane)"

        @pytest.fixture()
        def instance(self):
            mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
            Chem.AllChem.Compute2DCoords(mol)
            return CustomComponent(rdkit=mol, name="ethane")


Step 3: Define Required Methods
-----------------------------------

When inheriting from abstract base classes, such as ``Component``, you will need to define anything that is an ``abstractmethod``.
This includes both in ``Component`` itself, as well as any ``abstractmethod``\s it inherits from ``GufeTokenizable`` (since component is a subclass of GufeTokenizable).

In other cases, such as when inheriting from ``ExplicitMoleculeComponent``, you will only need to define methods specifically not implemented - in this case ``_to_dict()`` and ``_from_dict()``.

.. TODO: point to gufe tokenizable how-to for _to_dict, _from_dict examples, and/or source code examples.


Step 4: Define Additional Functionality
---------------------------------------

While the code in Step 2 is technically correct, it doesn't actually add anything new; it merely creates a new class identical to ``SmallMoleculeComponent`` with a new name.
To add functionality *in addition to* ``SmallMoleculeComponent``'s existing functionality, you can add new attributes, such as ``custom_attribute``, and new methods, such as ``print_custom_attribute``:

.. code-block:: python
    :caption: custom_component.py

    from gufe import SmallMoleculeComponent
    from rdkit import Chem


    class CustomComponent(SmallMoleculeComponent):
        def __init__(self, rdkit: Chem.rdchem.Mol, name: str = "", custom_attribute: int = 4):
            self.custom_attribute = custom_attribute
            super().__init__(rdkit=rdkit, name=name)

        def custom_functionality(self) -> str:
            return f"my custom attribute is {self.custom_attribute}"

        # Since we added a new attribute, must include that attribute in serialization
        # by defining _to_dict and _from_dict
        def _to_dict(self):
            # first, use the parent classes' implementation
            d = super()._to_dict()

            # now, add our custom attribute
            d["custom_attribute"] = self.custom_attribute
            return d

        @classmethod
        def _from_dict(cls, d:dict):
            # first, use the parent classes' implementation
            obj = super()._from_dict(d)
            # now, pass through to construct our custom attribute
            return cls(rdkit=obj._rdkit, name=obj.name, custom_attribute=d.get("custom_attribute"))

Just make sure you test all of your new features!

.. code-block:: python
    :caption: test_custom_component.py

    from gufe.tests import GufeTokenizableTestsMixin
    from .custom_component import CustomComponent
    import pytest
    from rdkit import Chem


    class TestCustomComponent(GufeTokenizableTestsMixin):
        cls = CustomComponent
        repr = "CustomComponent(name=ethane)"

        @pytest.fixture()
        def instance(self):
            mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
            Chem.AllChem.Compute2DCoords(mol)

            # it's important we test defining the custom attribute to make sure it round-trips correctly
            test_instance = CustomComponent(rdkit=mol, name="ethane", custom_attribute=8)
            return test_instance

        def test_print_custom_attribute_default(self, instance):
            # test using the instance fixture created above
            assert instance.print_custom_attribute() == "my custom attribute is 4"

        def test_print_custom_attribute_user_defined(self):
            # test by creating a new instance
            mol = Chem.AddHs(Chem.MolFromSmiles("CCC"))
            Chem.AllChem.Compute2DCoords(mol)
            custom_component = CustomComponent(
                rdkit=mol, name="propane", custom_attribute=7
            )
            assert custom_component.print_custom_attribute() == "my custom attribute is 7"
