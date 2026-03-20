.. _understanding_logging:

Logging in ``GufeTokenizable``\s
=================================

All :class:`.GufeTokenizable` objects in **gufe** have built-in logging support through the ``logger`` property.
This logging mechanism is designed to make debugging and runtime introspection easier by automatically tagging log messages with the object's :class:`.GufeKey`.

Basic Usage
-----------

Every ``GufeTokenizable`` has a ``logger`` property that returns a :class:`logging.LoggerAdapter` configured specifically for that instance:

.. code-block:: python

    from gufe import SolventComponent

    sc = SolventComponent()
    sc.logger.info("Creating solvent component")
    sc.logger.debug("Detailed information about the component")
    sc.logger.warning("Something to be aware of")
    sc.logger.error("An error occurred")

How It Works
------------

The logging system uses Python's standard :mod:`logging` module with a custom :class:`logging.LoggerAdapter`.
Each ``GufeTokenizable`` instance gets its own logger adapter that automatically adds the object's ``GufeKey`` to the logging context.

Logger Naming
^^^^^^^^^^^^^

Loggers are named hierarchically based on the object's module and class:

.. code-block:: python

    gufekey.<module>.<qualname>

For example, a :class:`.SolventComponent` has the logger name:

.. code-block:: python

    gufekey.gufe.components.solventcomponent.SolventComponent

This hierarchical naming allows you to control logging at different levels of granularity (e.g., all ``gufe`` objects, all components, or just ``SolventComponent`` instances).

GufeKey in Log Messages
^^^^^^^^^^^^^^^^^^^^^^^^

The logger adapter adds a ``gufekey`` field to each log record, containing the full :class:`.GufeKey` of the object (e.g., ``SolventComponent-26b4034ad9dbd9f908dfc298ea8d449f``).

This allows you to use ``%(gufekey)s`` in your `logging formatters <https://docs.python.org/3/library/logging.html#logging.Formatter>`_ to include the object's key in log messages.

Configuring Logging
-------------------

To see log messages from **gufe** objects, you need to configure the Python logging system.
Here's a complete example:

.. code-block:: python

    import logging
    import time
    from gufe import SolventComponent

    # Get the root logger for all gufe objects
    logger = logging.getLogger('gufekey')

    # Create a formatter that includes the GufeKey
    formatter = logging.Formatter(
        "[%(asctime)s] [%(gufekey)s] [%(levelname)s] %(message)s"
    )
    formatter.converter = time.gmtime  # Use UTC time for timestamps

    # Create a handler and attach the formatter
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    # Add the handler to the logger and set the level
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    # Now log messages will be visible
    sc = SolventComponent()
    sc.logger.info("Solvent component created")

This produces output like:

.. code-block:: text

    [2025-10-29 23:36:27,220] [SolventComponent-26b4034ad9dbd9f908dfc298ea8d449f] [INFO] Solvent component created

Logging Levels
--------------

The standard Python logging levels apply:

* ``DEBUG``: Detailed diagnostic information
* ``INFO``: General informational messages
* ``WARNING``: Warning messages for potentially problematic situations
* ``ERROR``: Error messages for serious problems
* ``CRITICAL``: Critical messages for very serious errors

You control which messages appear by setting the logger level:

.. code-block:: python

    # Show only warnings and errors
    logger.setLevel(logging.WARNING)

    # Show everything, including debug messages
    logger.setLevel(logging.DEBUG)

Best Practices
--------------

1. **Use appropriate log levels**: Reserve ``DEBUG`` for detailed diagnostic information, ``INFO`` for general progress updates, and ``WARNING``/``ERROR`` for problems.

2. **Don't check verbosity flags**: Unlike some older patterns, you should not check verbosity settings before logging.
   Instead, log freely at the appropriate level and let the logging configuration control what gets shown.

   .. code-block:: python

       # Good - just log at the appropriate level
       self.logger.debug("Detailed information")
       self.logger.info("Progress update")

       # Bad - don't do this
       if verbose:
           self.logger.info("Progress update")

3. **Configure at execution time**: The user or execution engine should configure logging levels when running code, not the library code itself.

4. **Use the GufeKey in formatters**: Include ``%(gufekey)s`` in your log formatters to track which object generated each message.
   This is especially useful when many objects are logging simultaneously.

5. **Log to the object's logger**: Always use ``self.logger`` within your ``GufeTokenizable`` methods, not a module-level logger or ``print()`` statements.

Advanced: Hierarchical Logger Control
--------------------------------------

Because logger names are hierarchical, you can control logging at different levels of specificity:

.. code-block:: python

    # Control all gufe logging
    logging.getLogger('gufekey').setLevel(logging.INFO)

    # Control only Protocol logging
    logging.getLogger('gufekey.gufe.protocols').setLevel(logging.DEBUG)

    # Control only a specific Protocol subclass
    logging.getLogger('gufekey.gufe.protocols.protocolunit').setLevel(logging.WARNING)

This allows fine-grained control over which objects produce log output without modifying the code itself.

Note: Logger Initialization Timing
-----------------------------------

When a ``GufeTokenizable`` is first created, its ``GufeKey`` may not yet be computed.
If you log a message before the key is available, the ``gufekey`` field will show ``"UNKNOWN"`` instead of the actual key.

.. code-block:: python

    class MyComponent(GufeTokenizable):
        def __init__(self, value):
            # This logs before the key is computed
            self.logger.info("Starting initialization")  # Shows "UNKNOWN"
            self.value = value
            # This logs after the key is computed
            self.logger.info("Initialization complete")  # Shows actual GufeKey

In practice, this is rarely an issue since most logging happens after initialization is complete.
