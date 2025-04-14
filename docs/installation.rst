Installation
============

.. toctree::

**gufe** is available through conda-forge and can be installed with mamba (or conda):

.. code:: bash

    $ mamba install -c conda-forge gufe


If you're a developer, you will want to build an editable install:

Clone the repository:

.. code:: bash
    $ git clone https://github.com/OpenFreeEnergy/gufe.git
    $ cd gufe


Create and activate a new environment:

.. TODO: use mamba update -f environment.yaml if they need to build in an existing environment?

.. code:: bash
    $ mamba create -f environment.yml
    $ mamba activate gufe

Build an editable installation:

.. code:: bash
    $ python -m pip install --no-deps -e .