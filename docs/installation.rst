Installation
============

.. toctree::

Installation with mamba
-----------------------

**gufe** is available through *conda-forge* and can be installed with mamba (or conda):

.. code-block:: bash

    $ mamba install -c conda-forge gufe


Optional Visualization Dependencies
-----------------------------------

**gufe** has optional software dependences which can be installed to visualize atom mappings.
This dependences are optional to keep **gufe** lightweight when used purely for compute.

**py3Dmol** is used to create three dimensional views of atom mappings and can be installed with:

.. code-block:: bash

    $ mamba install -c conda-forge py3dmol


**ipywidgets** is used to create a widget in a jupyter notebook to view atom mappings and can be installed with:

.. code-block:: bash

    $ mamba install -c conda-forge ipywidgets

For an optimal installation experience we recommend installing the optional packages at the same time you install **gufe**:

.. code-block:: bash

    $ mamba install -c conda-forge gufe py3dmol ipywidgets

While they can be installed after **gufe** is installed, we find that **mamba** has an easier time solving the environment when everything is installed at the same time.

Developer Installation
----------------------

If you're a developer, you will likely want to create a local editable installation.

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
