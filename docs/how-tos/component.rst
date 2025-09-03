.. _howto-component:

How to define a new ``Component``
=================================

The **gufe**  :ref:`Component<component>` is a GufeTokenizable intended to be used as an *extensible point* of the library, so that you can describe a simulated system in terms of GufeTokenizables that are compatible with **gufe** and the rest of the OpenFE ecosystem.

In many cases, you will likely want your custom Component to inherit from one of the following extensible points that inherit from the Component base class - meaning that all of the following are Components, but come with additional functionality:


Step 1: Choose a Component to Extend
-------------------------------------

.. TODO: link to gufe components toctree here

Step 2: Write Tests for Expected Behavior
-----------------------------------------

Step 3: Define the Template Methods
-----------------------------------


Step 4: Define Additional Functionality
---------------------------------------
