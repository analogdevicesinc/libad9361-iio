libad9361-iio Python Bindings
==================================

Python bindings for the `libad9361-iio`. A device-specific library for AD936X transceivers.

Installation
==================

The libad9361-iio python bindings can be installed from pip

.. code-block:: bash

  (sudo) pip install pylibad9361

or by grabbing the source directly

.. code-block:: bash

  git clone https://github.com/analogdevicesinc/libad9361-iio.git
  cd bindings/python
  (sudo) python3 setup.py install

.. note::

  On Linux the libad9361-iio python bindings are sometimes installed in locations not on path. On Ubuntu this is a common fix

  .. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:/usr/lib/python{python-version}/site-packages

.. toctree::
   :maxdepth: 1
   :caption: Contents:

Components
==================
.. toctree::
   :maxdepth: 1

   functions

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
