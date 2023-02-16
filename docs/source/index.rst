.. CasCy documentation master file, created by
   sphinx-quickstart on Thu Feb  9 01:24:05 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CasCy's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Plane-cylinder geometry
-----------------------

.. autoclass:: CasCy.plane_cylinder.plane_cylinder_system

Cylinder-cylinder geometry
--------------------------

.. autoclass:: CasCy.cylinder_cylinder.cylinder_cylinder_system


Cylinder reflection
-------------------

.. automodule:: CasCy.cylinder_reflection

.. autofunction:: CasCy.cylinder_reflection.pwrc_TMTM
.. autofunction:: CasCy.cylinder_reflection.pwrme
.. autofunction:: CasCy.cylinder_reflection.reflection_matrix

.. automodule:: CasCy.scattering_amplitude.scattering_amplitude

Bessel functions
----------------
.. automodule:: CasCy.bessel

.. autofunction:: CasCy.bessel.bessel_array

Quadrature rules
----------------
.. automodule:: CasCy.quadratures
.. autofunction:: CasCy.quadratures.fcqs_semiinfinite
.. autofunction:: CasCy.quadratures.fcqs_combined

Matrix operations
-----------------
.. automodule:: CasCy.matrix_operations
.. autofunction:: CasCy.matrix_operations.logdet1m



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
