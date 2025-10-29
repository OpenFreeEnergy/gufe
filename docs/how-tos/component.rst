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
    :caption: catmoleculecomponent.py

    from gufe import SmallMoleculeComponent

    # this is all the code you need to get the GufeTokenizableTestsMixin to pass,
    # since SmallMoleculeComponent fully implemented
    class CatMoleculeComponent(SmallMoleculeComponent):
        pass

.. code-block:: python
    :caption: test_catmoleculecomponent.py

    from .catcomponent import CatMoleculeComponent


    class TestCatMoleculeComponent(GufeTokenizableTestsMixin):
        cls = CatMoleculeComponent
        repr = "CatMoleculeComponent(name=ethane)"

        @pytest.fixture()
        def instance(self):
            mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
            Chem.AllChem.Compute2DCoords(mol)
            return CatMoleculeComponent(rdkit=mol, name="ethane")

.. TODO: how to include a new attribute in dict while keeping parent class behavior?

Step 3: Define the Required Methods
-----------------------------------

- you will need to define anything that is an ``abstractmethod``, both in ``Component`` itself, as well as the ``abstractmethod``\s that it inherits from ``GufeTokenizable`` (since component is a subclass of GufeTokenizable)
- this is a key benefit of instead inheriting from a class that is not an abstract base class, as many more methods will come pre-defined. For example, if you inherit from `ExplicitMoleculeComponent` you only have to define to/from  dict. For example:




Step 4: Define Additional Functionality
---------------------------------------

.. code-block:: python
    :caption: catmoleculecomponent.py

    class CatMoleculeComponent(SmallMoleculeComponent):
        def __init__(self, rdkit: Chem.rdchem.Mol, name:str="", is_hydrophobic:bool=True):
            self.is_hydrophobic = is_hydrophobic
            super().__init__(rdkit=rdkit, name=name)

.. code-block:: python
    :caption: test_catmoleculecomponent.py

    from .catcomponent import CatMoleculeComponent

    @pytest.fixture
    def cat_ethane():
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        Chem.AllChem.Compute2DCoords(mol)
        return CatMoleculeComponent(rdkit=mol, name="ethane")

    class TestCatMoleculeComponent(GufeTokenizableTestsMixin):
        cls = CatMoleculeComponent
        repr = "CatMoleculeComponent(name=ethane)"

        @pytest.fixture()
        def instance(self, cat_ethane):
            return cat_ethane

    def test_is_hydrophobic_default(cat_ethane):
        assert cat_ethane.is_hydrophobic is True
