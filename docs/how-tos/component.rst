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

As when defining any new ``GufeTokenizable``, you are encouraged to use the ``GufeTokenizableTestsMixin`` to ensure that your new ``Component`` works as intended.

Step 3: Define the Template Methods
-----------------------------------


Step 4: Define Additional Functionality
---------------------------------------
