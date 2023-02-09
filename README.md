# CasCy
A package to calculate the **Cas**imir energy between two parallel, dielectric **cy**linders of arbitrary radii (R1 and R2) immersed in an electrolyte solution.
The Casimir energy is calculated within the scattering theory using plane waves as a basis for the electromagnetic modes.
In this code, the so-called zero-frequency contribution to the Casimir energy is calculated, which becomes exact at large temperatures or
separations larger than the thermal wave-length.
The separation d between the cylinders is assumed to be much smaller than their lengths L, so that edge effects can be neglected.


<p align="center">
  <img src="images/geometry.svg" height="60%" width="60%" >
</p>

If one cylinder radius becomes infinitely large, it becomes a plane. This package also offers a code to calculate the Casimir energu in the limit of the cylinder-plane geometry.

## Installation

The package is written in python and cython.
It can be easily installed with the package installer [pip](https://pypi.org/project/pip/). First navigate to package folder
```
cd path_to/cascy_ded
```
and install the package with the command
```
pip install -e .
```

## Usage

The Casimir energy in units of k_B T (T is the system's temperature) between two cylinders of radius 3nm and length 15um at a separation of 6nm can be calculated with the following python code:

```
from cascy_ded.cylinder_cylinder import cylinder_cylinder_system
s = cylinder_cylinder_system(d=6.e-9, R1=3.e-9, R2=3.e-9, L=15.e-6)
print(s.calculate_casimir_energy())

>>> -5.12100852917314
```

The Casimir energy in units of k_B T (T is the system's temperature) between a cylinder of radius 3nm and length 15um at a separation of 6nm can be calculated with the following python code:

```
from cascy_ded.plane_cylinder import plane_cylinder_system
s = plane_cylinder_system(d=26.e-9, R=12.e-9, L=50.e-6)
print(s.calculate_casimir_energy())

>>> -16.777232984640037
```
Further examples can be found in the `examples/` folder. Notice that the package `matplotlib` is required to plot the results there.

## Documentation

(wip)
* make sphinx documentation, finish some of the docstrings
