.. _howto-protocol:

How to define a new ``Protocol``
================================

The **gufe** :ref:`Protocols <protocol>` system is designed as an *extensible point* of the library,
allowing developers to write their own methods for performing free energy calculations in a form that can take full advantage of the OpenFE ecosystem.

Our recommendation in this how-to is to start simple; we will call out areas where you can iterate further on making your ``Protocol`` more sophisticated.

.. image:: ../_static/gufe_protocol_diagram.svg
    :alt: The ``gufe`` protocol system.


Designing and defining your ``ProtocolDAG``
-------------------------------------------

First, subclass :class:`.Protocol`, and create the :meth:`.Protocol._create` method:

.. code-block:: python

    from gufe import ChemicalSystem, Protocol, ComponentMapping, ProtocolDAGResult, ProtocolUnit
    
    class MyProtocol(Protocol):

        def _create(
            self,
            stateA: ChemicalSystem,
            stateB: ChemicalSystem,
            mapping: Optional[Union[ComponentMapping, list[ComponentMapping]]] = None,
            extends: Optional[ProtocolDAGResult] = None,
        ) -> list[ProtocolUnit]:
        ...


This method should define how the :ref:`ProtocolDAG <protocoldag>` is constructed.
For example:

.. code-block:: python

    class MyProtocol(Protocol):
        
        def _create self, stateA, stateB, mapping, extends):

        # convert protocol inputs into starting points for independent simulations
        alpha = InitializeUnit(
            name="the beginning",
            settings=self.settings,
            stateA=stateA,
            stateB=stateB,
            mapping=mapping,
            start=starting_point,
            some_dict={"a": 2, "b": 12},
        )

        # create several units that would each run an independent simulation
        simulations: list[ProtocolUnit] = [
            SimulationUnit(settings=self.settings, name=f"sim {i}", window=i, initialization=alpha)
            for i in range(self.settings.n_repeats)  # type: ignore
        ]

        # gather results from simulations, finalize outputs
        omega = FinishUnit(settings=self.settings, name="the end", simulations=simulations)

        # return all `ProtocolUnit`s we created
        return [alpha, *simulations, omega]






Defining the ``ProtocolUnit``\s of your ``ProtocolDAG``
-------------------------------------------------------


Turning ``ProtocolDAGResult``\s into free energy estimates
----------------------------------------------------------


Finishing up
------------
* :meth:`.Protocol._default_settings`


Putting it all together
-----------------------

How to write tests for this...
