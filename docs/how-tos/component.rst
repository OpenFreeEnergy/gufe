.. _howto-component:

How to define a new ``Component``
=================================

The **gufe**  :ref:`Component <component>` is a :ref:`GufeTokenizable <understanding_gufetokenizables>` intended to be used as an *extensible point* of the library, so that you can describe a simulated system in terms of ``GufeTokenizable``\s that are compatible with **gufe** and the rest of the OpenFE ecosystem.

Our recommendation in this how-to is to start simple by extending an existing component class rather than building from scratch.


Overview: What makes a ``Component``
-------------------------------------

A ``Component`` represents an individual chemical entity within a :ref:`ChemicalSystem <chemicalsystem>`. All components share these characteristics:

* **Immutable and hashable**: Components can be used as dictionary keys and in sets
* **Serializable**: Components can be saved to and loaded from JSON or MessagePack format via :ref:`GufeTokenizable <understanding_gufetokenizables>`
* **Named**: Each component has a ``name`` property for identification
* **Charge-aware**: Components expose a ``total_charge`` property (if applicable)

The base :class:`.Component` class is abstract and requires implementing serialization methods (``_to_dict()``, ``_from_dict()``, ``_gufe_tokenize()``). However, **gufe** provides several fully-implemented subclasses:

* :class:`.SmallMoleculeComponent` - Small organic molecules with explicit coordinates
* :class:`.ProteinComponent` - Biological macromolecules
* :class:`.SolventComponent` - Solvent systems
* :class:`.ExplicitMoleculeComponent` - Base class for molecules with explicit coordinates

When creating custom components, you'll typically extend one of these concrete classes rather than the abstract :class:`.Component` base.


When to Create Custom Components
---------------------------------

You may want to create a custom ``Component`` subclass when:

* **Adding domain-specific metadata**: You need to attach additional properties to molecular components, such as experimental binding affinities, synthesis routes, or pharmacokinetic data that should be preserved through serialization.

* **Implementing specialized file format support**: You're working with a molecular representation format not currently supported by the existing components (e.g., specialized force field formats, quantum chemistry outputs, or domain-specific file types).

* **Extending molecules with computed properties**: You want to cache expensive calculations (like conformer energies, partial charges, or molecular descriptors) as part of the component itself, ensuring they persist across network serialization.

* **Creating non-standard chemical entities**: You need to represent chemical systems beyond the scope of current components, such as coarse-grained models, quantum mechanical systems, or multi-scale representations.

* **Building protocol-specific requirements**: Your custom :ref:`Protocol <protocol>` requires components with specific attributes or methods that aren't provided by the base classes.

In most cases, you'll extend an existing fully-implemented component like :class:`.SmallMoleculeComponent` rather than starting from the abstract :class:`.Component` base class. This approach lets you leverage existing serialization, hashing, and conversion functionality while adding only the specific features you need.


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

.. TODO: if we expose the other testing mixins, reference ExplicitMoleculeComponentMixin here.

.. code-block:: python
    :caption: custom_component.py

    from gufe import SmallMoleculeComponent

    # this is all the code you need to get the GufeTokenizableTestsMixin to pass,
    # since SmallMoleculeComponent is fully implemented
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


Step 3: Understand Required Methods
------------------------------------

The methods you need to implement depend on which base class you're extending:

**Extending fully-implemented classes** (like :class:`.SmallMoleculeComponent`):
  No required methods - you can add functionality without implementing anything. However, if you add custom attributes, you must implement ``_to_dict()``, ``_from_dict()``, and ``_gufe_tokenize()`` to handle serialization.

**Extending partially-implemented classes** (like :class:`.ExplicitMoleculeComponent`):
  Must implement ``_to_dict()`` and ``_from_dict()`` for serialization. You typically don't need to implement ``_gufe_tokenize()`` as it's handled by the parent.

**Extending abstract** :class:`.Component` **directly**:
  Must implement all abstract methods from both :class:`.Component` and :class:`.GufeTokenizable`:

  * ``name`` property
  * ``total_charge`` property
  * ``_to_dict()`` method
  * ``_from_dict()`` class method
  * ``_gufe_tokenize()`` method

For most use cases, extending :class:`.SmallMoleculeComponent` is the simplest path, as shown in the following step.

.. TODO: point to gufe tokenizable how-to for to_dict, from_dict examples, and/or source code examples.



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

        def print_custom_attribute(self) -> str:
            return f"my custom attribute is {self.custom_attribute}"


.. warning::
   If you add custom attributes like ``custom_attribute`` above, you must also implement
   ``_to_dict()``, ``_from_dict()``, and ``_gufe_tokenize()`` methods to properly handle
   serialization and hashing. See the :ref:`serialization how-to <howto-serialization>`
   for details. The example above is simplified to demonstrate adding methods; in practice,
   you may want to add only methods (not attributes) to avoid serialization complexity.

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


Using your Component
--------------------

Once you've created your custom ``Component``, you can use it anywhere a standard component is used:

.. code-block:: python

    from gufe import ChemicalSystem
    from rdkit import Chem

    # Create your custom component
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    Chem.AllChem.Compute2DCoords(mol)
    my_component = CustomComponent(rdkit=mol, name="ethanol", custom_attribute=42)

    # Use it in a ChemicalSystem
    system = ChemicalSystem({"ligand": my_component})

    # Access the custom functionality
    print(my_component.print_custom_attribute())  # "my custom attribute is 42"

    # Serialize and deserialize
    serialized = my_component.to_dict()
    restored = CustomComponent.from_dict(serialized)
    assert restored.custom_attribute == 42

    # Use in transformations and networks
    from gufe import Transformation
    system_a = ChemicalSystem({"ligand": my_component})
    system_b = ChemicalSystem({"ligand": another_component})
    transformation = Transformation(
        stateA=system_a,
        stateB=system_b,
        protocol=my_protocol,
        mapping=None
    )


Best practices and tips
-----------------------

1. **Extend, don't start from scratch**: Inherit from :class:`.SmallMoleculeComponent`, :class:`.ProteinComponent`, or other concrete classes rather than the abstract :class:`.Component` base. This gives you serialization, hashing, and conversion methods for free.

2. **Be conservative with custom attributes**: Adding instance attributes requires implementing ``_to_dict()``, ``_from_dict()``, and ``_gufe_tokenize()``. Consider whether you can achieve your goals by adding methods instead.

3. **Maintain immutability**: Components should be immutable after creation. Don't add setters or methods that modify state. This ensures they remain hashable and safe to use as dictionary keys.

4. **Test serialization thoroughly**: Use :class:`~gufe.tests.GufeTokenizableTestsMixin` to ensure your component serializes correctly. Test that ``component.to_dict()`` and ``ComponentClass.from_dict()`` round-trip properly.

5. **Document your additions**: Add clear docstrings explaining what your custom attributes and methods do, especially if they'll be used by others.

6. **Consider backwards compatibility**: If you modify serialization, consider implementing ``serialization_migration()`` to handle loading old serialized versions.

7. **Handle RDKit carefully**: If working with RDKit molecules, remember:

   * Set ``Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)`` to preserve molecule properties
   * Ensure molecules have at least one conformer before creating a component
   * Be aware that only the first conformer is used if multiple are present

8. **Name your components meaningfully**: The ``name`` property is used in ``repr()`` and helps with debugging. Choose descriptive names.

9. **Validate inputs in ``__init__``**: Check that inputs are valid before calling ``super().__init__()``. Raise clear exceptions for invalid inputs.


Common pitfalls
---------------

* **Forgetting to implement serialization**: If you add custom attributes, you must implement ``_to_dict()``, ``_from_dict()``, and ``_gufe_tokenize()``, or serialization will fail or lose data.

* **Breaking immutability**: Adding methods that modify state breaks the component's hashability and can cause subtle bugs in networks and sets.

* **Not preserving RDKit properties**: RDKit doesn't preserve molecule properties during pickling by default. Use ``Chem.PropertyPickleOptions.AllProps`` if you need to preserve partial charges or other properties.

* **Inconsistent hashing**: If ``_gufe_tokenize()`` returns different values for equivalent objects, you'll get duplicate components in memory.

* **Missing tests**: Skipping the :class:`~gufe.tests.GufeTokenizableTestsMixin` tests means serialization bugs may only appear at runtime, often in production.
