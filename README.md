# CasCy
A package to calculate the **Cas**imir energy between two parallel, dielectric **cy**linders of arbitrary radii (R1 and R2) immersed in an electrolyte solution.
The Casimir energy is calculated within the scattering theory using plane waves as a basis for the electromagnetic modes.
In this code, the so-called zero-frequency contribution to the Casimir energy is calculated, which becomes exact at large temperatures or
separations larger than the thermal wave-length.
The separation d between the cylinders is assumed to be much smaller than their lengths L, so that edge effects can be neglected.


<p align="center">
  <img src="images/geometry.svg" height="60%" width="60%" >
</p>

If one cylinder radius becomes infinitely large, it becomes a plane. This package also offers a code to calculate the Casimir energy in the limit of the cylinder-plane geometry.

## Installation

CasCy is written in Python and uses Cython for performance improvements. Note that Cython requires a C compiler to be present on the system.
The package can be easily installed with the package installer [pip](https://pypi.org/project/pip/). First navigate to package folder
```
cd path/to/CasCy
```
If you do not wish to install the package into your base python library, you may want to create a [virtual environment](https://docs.python.org/3/tutorial/venv.html) before installing the package. A virtual environment can be created with
```
python -m venv env
```
To activate the environment on Unix/MacOS, run:
```
source env/bin/activate
```
or on Windows, run:
```
env/Scripts/activate
```

Finally, install the package with the command
```
pip install .
```

## Usage

The Casimir energy in units of k_B T (with the system's temperature T) between two cylinders of radius 3nm and length 15um at a separation of 6nm can be calculated with the following python code:

```
from CasCy.cylinder_cylinder import cylinder_cylinder_system
s = cylinder_cylinder_system(d=6.e-9, R1=3.e-9, R2=3.e-9, L=15.e-6)
print(s.calculate_casimir_energy())

>>> -5.12100852917314
```
By default, the resulting value is accurate to about 1%. Thera are optional arguments for the function
`calculate_casimir_energy()` which can be tuned to improve the numerical accuracy. For more information see the
documentation.


The Casimir energy in units of k_B T (with the system's temperature T) between a plane and a cylinder of radius 3nm and length 15um at a separation of 6nm can be calculated with the following python code:

```
from CasCy.plane_cylinder import plane_cylinder_system
s = plane_cylinder_system(d=26.e-9, R=12.e-9, L=50.e-6)
print(s.calculate_casimir_energy())

>>> -16.777232984640037
```
Further examples can be found in the `examples` folder. Notice that the package `matplotlib` is required to plot the results there.

## Documentation

The documentation contains more information about the classes and functions defined in this package.
It can be built using [sphinx](https://www.sphinx-doc.org/en/master/).
First navigate to the documentation folder
```
cd docs/
```
and execute the command
```
make html
```
to build the documentation. The html files can be found in the folder `docs/build/html` and they can be viewed with any standard web browser.

## Citation

If you use this software in your research, please cite the associated paper:

B. Spreng, H. Berthoumieux, A. Lambrecht, A.-F. Bitbol, P. A. Maia Neto, S. Reynaud, "Universal {{Casimir}} Attraction between Filaments at the Cell Scale," New J. Phys. 26, 013009 (2024). DOI: [10.1088/1367-2630/ad1846](https://doi.org/10.1088/1367-2630/ad1846)

### BibTeX

For those who wish to cite the paper in LaTeX:

```bibtex
@article{SprengUniversalCasimirAttraction2024,
  title={Universal {{Casimir}} Attraction between Filaments at the Cell Scale},
  author={Spreng, B. and Berthoumieux, H. and Lambrecht, A. and Bitbol, A.-F. and Maia Neto, P. A. and Reynaud, S.},
  journal={New J. Phys.},
  volume={26},
  number={1},
  pages={013009},
  year={2024},
  doi={10.1088/1367-2630/ad1846}
}

