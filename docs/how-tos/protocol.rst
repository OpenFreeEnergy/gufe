.. _howto-protocol:

How to define a new ``Protocol``
================================

The **gufe** :ref:`Protocol <protocol>` system is designed as an *extensible point* of the library,
allowing developers to write their own methods for performing calculations in a form that can take full advantage of the OpenFE ecosystem.

Our recommendation in this how-to is to start simple;
we will call out areas where you can iterate further on making your ``Protocol`` more sophisticated.

.. image:: ../_static/gufe_protocol_diagram.svg
    :alt: The ``gufe`` protocol system.


Overview: What makes a ``Protocol``
-----------------------------------

A complete ``Protocol`` implementation requires the following interconnected components:

1. **Protocol class**: The main class that creates computational workflows
2. **Settings class**: A `Pydantic <https://docs.pydantic.dev/latest/>`_ model defining configuration parameters
3. **ProtocolUnit classes**: Individual computational steps in the workflow
4. **ProtocolResult class**: Container for aggregated results from multiple runs

The ``Protocol`` creates a :ref:`ProtocolDAG <protocoldag>` (directed acyclic graph) of :ref:`ProtocolUnit <protocolunit>` objects that define the computational workflow.
Multiple ``ProtocolDAG`` executions can be aggregated into a single :ref:`ProtocolResult <protocolresult>` to provide statistical estimates.

If this is your first time writing a ``Protocol``, we recommend defining these components in the following step order.
You can of course iterate on these components as you write each one, improving the implementation as you go.


Step 1: Define your Settings
-----------------------------

First, create a settings class that inherits from :class:`.Settings`.
This defines all the configuration parameters your protocol needs:

.. code-block:: python

    from pydantic import InstanceOf

    from gufe.settings import Settings, SettingsBaseModel
    from gufe.settings.typing import NanosecondQuantity


    class SimulationSettings(SettingsBaseModel):
        """Settings for simulation control, including lengths, etc...

        """
        minimization_steps: int = 5000
        """Number of minimization steps to perform. Default 5000."""

        equilibration_length: NanosecondQuantity
        """Length of the equilibration phase in units of time."""

        production_length: NanosecondQuantity
        """Length of the production phase in units of time."""


    class AlchemicalSettings(SettingsBaseModel):
        """Settings specific to the particulars of this alchemical method.

        """
        ...


    class MyProtocolSettings(Settings):
        """Settings for an example Protocol.

        """
        # the number of independent simulations to perform
        n_repeats: int = 5

        simulation_settings: SimulationSettings
        alchemical_settings: AlchemicalSettings

        # by subclassing from ``Settings``, we also get:
        # - forcefield_settings : parameters for force field choices
        # - thermo_settings : thermodynamics parameters, e.g. pressure, temperature


Some notes on the above:

1. **gufe** includes several :class:`~gufe.settings.typing.GufeQuantity` types,
   including the :class:`~gufe.settings.typing.NanosecondQuantity` used above.
   We recommend using these for settings fields that carry units,
   and you can easily make your own if necessary by following the pattern in the :mod:`gufe.settings.typing` module.

2. We definied a couple :class:`~gufe.settings.models.SettingsBaseModel` subclasses to group together related settings,
   such as the number of steps to use for various portions of the simulation in ``SimulationSettings``.
   It is common practice to break a ``Protocol`` ``Settings`` object up in this way to make them more modular and easier to work with.

3. Our :class:`~gufe.settings.models.Settings` subclass ``MyProtocolSettings`` will then feature a hierarchy of settings:
    - ``simulation_settings``: 
    - ``alchemical_settings``
    - ``forcefield_settings``
    - ``thermo_settings``


Step 2: Define your ProtocolResult
----------------------------------

Create a :ref:`ProtocolResult <protocolresult>` subclass that defines how to compute estimates and uncertainties from your ``Protocol``'s outputs:

.. code-block:: python

    from gufe import ProtocolResult
    from openff.units import unit
    import numpy as np

    class MyProtocolResult(ProtocolResult):
        
        # required method
        # return ``None`` if Protocol doesn't produce an estimate
        def get_estimate(self) -> unit.Quantity:
            """Calculate the free energy estimate from all runs."""
            # extract the key results from all completed runs
            free_energies = self.data["free_energies"]

            # get unit of the first value
            u = free_energies[0].u
            
            # return the mean as our best estimate, converting to same units
            return np.mean(np.asarray([dG.to(u).m for dG in free_energies])) * u
        
        # required method
        # return ``None`` if Protocol doesn't produce an estimate
        def get_uncertainty(self) -> unit.Quantity:
            """Calculate the uncertainty from all runs."""
            free_energies = self.data["free_energies"]

            # get unit of the first value
            u = free_energies[0].u
            
            # return the standard error as our uncertainty, converting to the same units
            std_dev = np.std(np.asarray([dG.to(u).m for dG in free_energies])) * u
            std_err = std_dev / np.sqrt(len(free_energies))
            return std_err


It's okay for the implementations of these methods to be a guess for now.
We will define how the contents of ``MyProtocolResult.data`` are assembled in :ref:`howto-protocol-protocol-class`.

Some additional notes:

1. The example above returns a single estimate by taking the sample mean of the individual :ref:`ProtocolDAGResult <protocoldagresult>` estimates,
   and reports the uncertainty in that single estimate as the standard error.
   This is not the only choice available.
   Some ``Protocol``\s choose to report their uncertainty as the standard deviation.
   Other ``Protocol``\s don't report the estimate as a sample mean at all, instead aggregating trajectories of reduced potentials
   or nonequilibrium works and deriving a single estimate from these using estimators such as BAR or MBAR.
   The choice is yours as to what is most appropriate for your ``Protocol``.

2. Although less common, you can also write ``Protocol``\s for :ref:`NonTransformations <nontransformation>`.
   These operate on a single :ref:`ChemicalSystem <chemicalsystem>`,
   and typically perform some form of sampling, such as equilibrium molecular dynamics.
   In this case, the :meth:`~gufe.protocols.protocol.ProtocolResult.get_estimate` and :meth:`~gufe.protocols.protocol.ProtocolResult.get_uncertainty` methods typically lack meaning, and can be writtent to return ``None``.
   It may make sense to create other methods returning quantities of interest from this sampling, however.


Step 3: Define your ProtocolUnits
----------------------------------

Create the computational units that will perform the actual work.
Each unit should inherit from :ref:`ProtocolUnit <protocolunit>` and implement an ``_execute`` method.

.. important :: 

   Use ``ctx.shared`` for large objects that need to be passed between units.
   This avoids serialization issues and improves performance by keeping file paths in the return objects instead of the large objects themselves.

.. code-block:: python

    from gufe import ProtocolUnit

    class SetupUnit(ProtocolUnit):
        """Prepare the system for simulation."""
        
        @staticmethod
        def _execute(ctx, *, stateA, stateB, mapping, settings, **inputs):
            """Set up the alchemical system."""
            import pickle
            from pathlib import Path
            
            # ctx provides scratch and shared directories
            # Use ctx.shared to write files that other units will need
            shared_dir = Path(ctx.shared)
            
            # Your setup logic here...
            # - Create alchemical system from stateA/stateB
            # - Apply the atom mapping 
            # - Set up force field parameters
            prepared_system = ...  # Your setup code here
            topology = ...         # Your topology creation
            coordinates = ...      # Your coordinate preparation
            
            # Write large objects to shared directory instead of returning them
            system_file = shared_dir / "system.pkl"
            topology_file = shared_dir / "topology.pkl" 
            coords_file = shared_dir / "initial_coords.pkl"
            
            with open(system_file, 'wb') as f:
                pickle.dump(prepared_system, f)
            with open(topology_file, 'wb') as f:
                pickle.dump(topology, f)
            with open(coords_file, 'wb') as f:
                pickle.dump(coordinates, f)
            
            return {
                "system_file": str(system_file),
                "topology_file": str(topology_file),
                "initial_coordinates_file": str(coords_file),
                "log": "System setup completed"
            }

    class SimulationUnit(ProtocolUnit):
        """Run an individual simulation."""
        
        @staticmethod 
        def _execute(ctx, *, setup_result, window, settings, **inputs):
            """Execute a single alchemical window simulation."""
            import pickle
            from pathlib import Path
            
            # Load large objects from files written by setup unit
            with open(setup_result.outputs["system_file"], 'rb') as f:
                system = pickle.load(f)
            with open(setup_result.outputs["topology_file"], 'rb') as f:
                topology = pickle.load(f)
            with open(setup_result.outputs["initial_coordinates_file"], 'rb') as f:
                coordinates = pickle.load(f)
            
            # Your simulation logic here...
            # - Run equilibration for `settings.equilibration_length`
            # - Run production simulation for `settings.simulation_length`
            # - Calculate free energy contribution
            dg_contribution = ...  # Your free energy calculation
            final_coords = ...     # Your final coordinates
            
            # Write output files to shared directory
            shared_dir = Path(ctx.shared)
            final_coords_file = shared_dir / f"final_coords_window_{window}.pkl"
            
            with open(final_coords_file, 'wb') as f:
                pickle.dump(final_coords, f)
            
            return {
                "free_energy_contribution": dg_contribution,
                "final_coordinates_file": str(final_coords_file),
                "window": window,
                "log": f"Simulation completed for window {window}"
            }

    class AnalysisUnit(ProtocolUnit):
        """Analyze results from all simulations."""
        
        @staticmethod
        def _execute(ctx, *, simulation_results, settings, **inputs):
            """Combine results from all simulation windows."""
            import pickle
            from pathlib import Path
            
            # simulation_results will be a list of ProtocolUnitResult objects
            total_free_energy = 0.0
            logs = []
            final_coords = {}
            
            for sim_result in simulation_results:
                # Extract numerical results directly
                total_free_energy += sim_result.outputs["free_energy_contribution"]
                logs.append(sim_result.outputs["log"])
                
                # Load coordinate files if needed for analysis
                window = sim_result.outputs["window"]
                coords_file = sim_result.outputs["final_coordinates_file"]
                with open(coords_file, 'rb') as f:
                    coords = pickle.load(f)
                final_coords[window] = coords
            
            # Write final analysis results to shared directory
            shared_dir = Path(ctx.shared)
            analysis_file = shared_dir / "final_analysis.pkl"
            
            analysis_data = {
                "total_free_energy": total_free_energy,
                "final_coordinates": final_coords,
                "simulation_logs": logs
            }
            
            with open(analysis_file, 'wb') as f:
                pickle.dump(analysis_data, f)
            
            return {
                "total_free_energy": total_free_energy,
                "analysis_file": str(analysis_file),
                "log": "Analysis completed"
            }


.. _howto-protocol-protocol-class:

Step 4: Implement your Protocol class
-------------------------------------

Now create your main custom ``Protocol`` class that inherits from :ref:`Protocol <protocol>` and ties everything together:

.. code-block:: python

    from gufe import Protocol, ChemicalSystem, ComponentMapping, ProtocolDAGResult, ProtocolUnit
    from typing import Optional, Union, List, Iterable, Any

    class MyProtocol(Protocol):
        # Required class attributes
        result_cls = MyProtocolResult
        _settings_cls = MyProtocolSettings
        
        @classmethod
        def _default_settings(cls) -> MyProtocolSettings:
            """Provide sensible default settings."""
            return MyProtocolSettings(
                # Include any thermodynamic settings your protocol needs
                # thermo_settings=ThermoSettings(temperature=298.15 * unit.kelvin),
                # Protocol-specific defaults
                n_repeats=5,
                simulation_length=10.0 * unit.nanosecond,
                equilibration_length=1.0 * unit.nanosecond
            )
        
        def _create(
            self,
            stateA: ChemicalSystem,
            stateB: ChemicalSystem,
            mapping: Optional[Union[ComponentMapping, List[ComponentMapping]]] = None,
            extends: Optional[ProtocolDAGResult] = None,
        ) -> List[ProtocolUnit]:
            """Create the computational workflow."""
            
            # Handle extension from previous results if needed
            if extends is not None:
                # Extract useful information from the previous run
                # This might be final coordinates, equilibrated structures, etc.
                starting_point = extends.protocol_unit_results[-1].outputs
            else:
                starting_point = None
            
            # Create the setup unit (runs once)
            setup = SetupUnit(
                name="system_setup",
                stateA=stateA,
                stateB=stateB,
                mapping=mapping,
                settings=self.settings,
                starting_point=starting_point
            )
            
            # Create multiple independent simulation units
            simulations = []
            for i in range(self.settings.n_repeats):
                sim_unit = SimulationUnit(
                    name=f"simulation_{i}",
                    setup_result=setup,  # This creates the dependency
                    window=i,
                    settings=self.settings
                )
                simulations.append(sim_unit)
            
            # Create analysis unit that depends on all simulations
            analysis = AnalysisUnit(
                name="final_analysis",
                simulation_results=simulations,  # Depends on all simulations
                settings=self.settings
            )
            
            # Return all units - dependencies are implicit from constructor args
            return [setup, *simulations, analysis]
        
        def _gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> dict[str, Any]:
            """Aggregate results from multiple ProtocolDAG executions."""
            # This method combines results from multiple independent protocol runs
            # into data that the ProtocolResult can use to compute estimates
            
            free_energies = []
            all_logs = []
            
            for dag_result in protocol_dag_results:
                # Find the terminal (final) unit results
                for unit_result in dag_result.terminal_protocol_unit_results:
                    if unit_result.name == "final_analysis":
                        free_energies.append(
                            unit_result.outputs["total_free_energy"]
                        )
                        all_logs.extend(unit_result.outputs["simulation_logs"])
            
            return {
                "free_energies": free_energies,
                "logs": all_logs
            }


Step 5: Add validation (optional)
----------------------------------

You can add custom validation to check that inputs are compatible with your protocol:

.. code-block:: python

    class MyProtocol(Protocol):
        # ... other methods ...
        
        def _validate(
            self, 
            *, 
            stateA: ChemicalSystem, 
            stateB: ChemicalSystem, 
            mapping: Optional[Union[ComponentMapping, List[ComponentMapping]]] = None,
            extends: Optional[ProtocolDAGResult] = None
        ):
            """Validate inputs for this protocol."""
            from gufe.protocols.errors import ProtocolValidationError
            
            # Check that both states have the required components
            if not stateA.components or not stateB.components:
                raise ProtocolValidationError("Empty chemical systems not supported")
            
            # Check that we have a mapping if needed
            if mapping is None:
                raise ProtocolValidationError("This protocol requires atom mappings")
            
            # Check ability to extend from given results
            if extends and not extends.ok():
                raise ProtocolValidationError("Cannot extend from failed ProtocolDAGResult")


Understanding ProtocolUnit dependencies
---------------------------------------

Dependencies between ``ProtocolUnit`` objects are established implicitly by passing one unit as a constructor argument to another:

.. code-block:: python

    # setup runs first (no dependencies)
    setup = SetupUnit(name="setup", ...)
    
    # simulation depends on setup (setup passed as argument)
    simulation = SimulationUnit(name="sim", setup_result=setup, ...)
    
    # analysis depends on simulation (simulation passed as argument)
    analysis = AnalysisUnit(name="analysis", simulation_results=[simulation], ...)

``ProtocolUnit`` objects can also be nested in dictionaries and lists, and dependencies will still be detected:

.. code-block:: python

    # Dependencies work when units are in lists
    simulations = [sim1, sim2, sim3]
    analysis = AnalysisUnit(name="analysis", simulations=simulations, ...)
    
    # Dependencies work when units are in dictionaries  
    unit_dict = {"equilibration": eq_unit, "production": prod_unit}
    final_unit = FinalUnit(name="final", inputs=unit_dict, ...)

The ``ProtocolDAG`` automatically determines the execution order from these dependencies.
Units with no dependencies run first, followed by units whose dependencies have completed.


Putting it all together: A complete example
--------------------------------------------

Here's a simplified but complete protocol implementation:

.. code-block:: python

    from gufe import Protocol, ProtocolUnit, ProtocolResult
    from gufe.settings import Settings
    from openff.units import unit
    from typing import Iterable, Any, List
    import numpy as np

    # Settings
    class SimpleProtocolSettings(Settings):
        n_repeats: int = 3

    # Result
    class SimpleProtocolResult(ProtocolResult):
        def get_estimate(self):
            return np.mean(self.data["values"]) * unit.kilocalorie_per_mole
        
        def get_uncertainty(self):
            values = self.data["values"]
            if len(values) < 2:
                return 0.0 * unit.kilocalorie_per_mole
            return np.std(values) / np.sqrt(len(values)) * unit.kilocalorie_per_mole

    # Units
    class SimpleUnit(ProtocolUnit):
        @staticmethod
        def _execute(ctx, **inputs):
            # Simulate a calculation that returns a random result
            result = np.random.normal(5.0, 1.0)  # Mean=5, std=1
            return {"result": result}

    # Protocol
    class SimpleProtocol(Protocol):
        result_cls = SimpleProtocolResult
        _settings_cls = SimpleProtocolSettings
        
        @classmethod
        def _default_settings(cls):
            return SimpleProtocolSettings(n_repeats=3)
        
        def _create(self, stateA, stateB, mapping=None, extends=None) -> List[ProtocolUnit]:
            # Create n_repeats independent units
            units = [
                SimpleUnit(name=f"calc_{i}", replica=i, settings=self.settings)
                for i in range(self.settings.n_repeats)
            ]
            return units
        
        def _gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> dict[str, Any]:
            values = []
            for dag_result in protocol_dag_results:
                for unit_result in dag_result.protocol_unit_results:
                    values.append(unit_result.outputs["result"])
            return {"values": values}


Using your Protocol
-------------------

Once implemented, your protocol can be used like any other **gufe** protocol:

.. code-block:: python

    # Create protocol with custom settings
    settings = MyProtocolSettings(n_repeats=10, simulation_length=20*unit.nanosecond)
    protocol = MyProtocol(settings)
    
    # Create a ProtocolDAG for specific chemical systems
    dag = protocol.create(
        stateA=chem_system_a,
        stateB=chem_system_b,
        mapping=atom_mapping
    )
    
    # Execute on a scheduler (not shown)
    # dag_result = scheduler.execute(dag)
    
    # Gather multiple results into final estimate
    # final_result = protocol.gather([dag_result1, dag_result2, ...])


Best practices and tips
-----------------------

1. **Start simple**: Begin with a minimal working implementation and add complexity gradually.

2. **Handle errors gracefully**: Use ``try``/``except`` in ``_execute`` methods and return meaningful error information.

3. **Use the context effectively**: The ``ctx`` parameter provides ``scratch`` (temporary, persists over execution of a single ``ProtocolUnit``) and ``shared`` (persists over execution of the ``ProtocolDAG``) directories.
   Use ``ctx.shared`` for large objects that need to pass between units; store file paths in return objects, not the objects themselves.

4. **Test thoroughly**: Write unit tests for your ``ProtocolUnit`` classes early in development.

5. **Document your settings**: Use Pydantic's `Field() function <https://docs.pydantic.dev/latest/concepts/fields/>`_ with descriptions to document what each setting does.

6. **Consider serialization**: All your classes should be serializable - avoid complex objects that can't be serialized with ``GufeTokenizable.to_json``.

7. **Resource management**: Clean up temporary files in your ``_execute`` methods when possible.

8. **Validate early**: Implement ``_validate`` to catch configuration problems before expensive computations begin.


Testing your Protocol
----------------------

Create unit tests for each component:

.. code-block:: python

    def test_protocol_creation():
        """Test that the protocol can be created with default settings."""
        protocol = MyProtocol(MyProtocol.default_settings())
        assert isinstance(protocol.settings, MyProtocolSettings)

    def test_dag_creation(sample_chemical_systems):
        """Test ProtocolDAG creation."""
        protocol = MyProtocol(MyProtocol.default_settings())
        dag = protocol.create(
            stateA=sample_chemical_systems[0],
            stateB=sample_chemical_systems[1],
            mapping=sample_mapping
        )
        
        assert len(dag.protocol_units) > 0
        # Test that dependencies are set up correctly
        
    def test_unit_execution():
        """Test individual ProtocolUnit execution."""
        from gufe.protocols.protocolunit import Context
        
        unit = SimpleUnit(name="test", replica=0, settings=SimpleProtocolSettings())
        
        # Mock context and inputs
        ctx = Context(scratch="/tmp", shared="/tmp")
        result = unit._execute(ctx, replica=0)
        
        assert "result" in result
        assert isinstance(result["result"], float)
