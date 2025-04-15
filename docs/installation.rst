Installation
============

.. toctree::

Installation with mamba
-----------------------

**gufe** is available through conda-forge and can be installed with mamba (or conda):

.. code-block:: bash

    $ mamba install -c conda-forge gufe


Developer Installation
----------------------

If you're a developer, you will want to build a local editable installation.

1. clone the repository:

.. code-block:: bash

    $ git clone https://github.com/OpenFreeEnergy/gufe.git
    $ cd gufe


2. create and activate a new environment:

.. code-block:: bash

    $ mamba create -f environment.yml
    $ mamba activate gufe

3. build an editable installation:

.. code-block:: bash

    $ python -m pip install --no-deps -e .
