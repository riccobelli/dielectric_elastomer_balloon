# Stability of a dielectric elastomer balloon covered by a passive elastic layer
Source code for the stability analysis of a spherical shell composed of a dielectric elastomer and an elastic layer.

## Dependencies

The code is written in Python 3 and tested with version 3.10. The following additional libraries are required, in the parentheses we indicate the version used in the simulations reported in the paper:
* FEniCS (https://fenicsproject.org/)
* Numpy (https://numpy.org/)
* Scipy (https://scipy.org/)
* gmsh (https://gmsh.info/)
* BiFEniCS (https://github.com/riccobelli/bifenics)

## Repository structure

The repository is organized as follows:
* `spherical_de.py` contains the main class which contains all the information about the numerical problem.
* `fixed_potential.py` performs first a ramp to increase the voltage and then applies an arclength continuation to investigate the buckling induced by an applied pression.
* `mesh.geo` is a file containing the geometrical parameters of the mesh. Use gmsh and the tool `dolfin-convert` to generate the file `mesh.xml` used by FEniCS to import the mesh.
```
gmsh -2 -format msh2 mesh.geo
dolfin-convert mesh.msh mesh.xml
```

## Citing

If you find this code useful for your work, please cite [[1]](#1)

## Licence

The source code contained in this repository is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

## References
<a id="1">[1]</a>
Y. Su, D. Riccobelli, Y. Chen, W. Chen, P. Ciarletta (2023), *under review*.

## Author
Davide Riccobelli, MOX - Politecnico di Milano (<davide.riccobelli@polimi.it>)
