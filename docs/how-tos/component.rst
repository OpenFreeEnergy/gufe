.. _howto-component:

How to define a new ``Component``
=================================

The **gufe**  :ref:`Component<component>` is a GufeTokenizable intended to be used as an *extensible point* of the library, so that you can describe a simulated system in terms of GufeTokenizables that are compatible with **gufe** and the rest of the OpenFE ecosystem.


Step 1: Choose a Component to Extend
-------------------------------------

In many cases, you will likely want your custom ``Component`` to inherit from one of the following extensible points that inherit from the ``Component`` base class - meaning that all of the following are ``Component``\s, but come with additional functionality:

* :class:`.ExplicitMoleculeComponent`
* :class:`.ProteinComponent`
* :class:`.SmallMoleculeComponent`
* :class:`.SolventComponent`

In the rare case where you want to custom-define as much as possible, you can inherit directly from :class:`.Component` itself.

Step 2: Write Tests for Expected Behavior
-----------------------------------------

As when defining any new ``GufeTokenizable``, you are encouraged to use the ``GufeTokenizableTestsMixin`` pytest fixture to ensure that your new ``Component`` works as intended.

For example:

.. code-block:: python
    :caption: custom_component.py

    from gufe import SmallMoleculeComponent

    # this is all the code you need to get the GufeTokenizableTestsMixin to pass,
    # since SmallMoleculeComponent fully implemented
    class CustomComponent(SmallMoleculeComponent):
        pass

.. code-block:: python
    :caption: test_custom_component.py

    from .custom_component import CustomComponent


    class TestCustomComponent(GufeTokenizableTestsMixin):
        cls = CustomComponent
        repr = "CustomComponent(name=ethane)"

        @pytest.fixture()
        def instance(self):
            mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
            Chem.AllChem.Compute2DCoords(mol)
            return CustomComponent(rdkit=mol, name="ethane")


Step 3: Define any Required Methods
-----------------------------------

When inheriting from abstract base classes, such as ``Component``, you will need to define anything that is an ``abstractmethod``, both in ``Component`` itself, as well as the ``abstractmethod``\s that it inherits from ``GufeTokenizable`` (since component is a subclass of GufeTokenizable).

In other cases, such as inheriting from ``ExplicitMoleculeComponent``, you will only need to define methods specifically not implemented - in this case ``to_dict()`` and ``from_dict()``.

.. TODO: point to gufe tokenizable how-to for to_dict, from_dict examples, and/or source code examples.

A key benefit of instead inheriting from a class that is already a completely implemented class, as shown next, since you will only need to implement the extra functionality you want.



Step 4: Define Additional Functionality
---------------------------------------


While the code in Step 2 is technically correct, it doesn't actually add anything new.
To add functionality in addition to ``SmallMoleculeComponent``'s existing functionality, you can add new attributes, such as ``custom_attribute``, and new methods, such as ``print_custom_attribute``:

.. code-block:: python
    :caption: custom_component.py

    from gufe import SmallMoleculeComponent
    from rdkit import Chem


    class CustomComponent(SmallMoleculeComponent):
        def __init__(self, rdkit: Chem.rdchem.Mol, name: str = "", custom_attribute: int = 4):
            self.custom_attribute = custom_attribute
            super().__init__(rdkit=rdkit, name=name)

        def print_custom_attribute(self) -> str:
            return f"my custom attribute is {self.custom_attribute}"


Just make sure you test any added features!

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
            test_instance = CustomComponent(rdkit=mol, name="ethane")
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
