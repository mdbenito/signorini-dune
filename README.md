# Signorini-DUNE

This is a finite element implementation of the problem of contact of two
elastic bodies in the context of classical linear elasticity [6]. We use a
mortar method [2] for the discretization and an active/inactive dual set
strategy [1] for the resolution of the non-linearity stemming from the 
non-penetration condition.

The former was first introduced in the 1990s as a domain decomposition
technique for non-matching grids enabling the gluing of different
discretizations at the common interface, and has since proved 
to have wide applicability in contact problems. The latter is an iterative
scheme to reduce the non-linearity at the contact zone to a choice between
Neumann and Dirichlet boundary conditions at each node.


## Dependencies

The project was developed using Apple's XCode, but it should be easy
to write a CMake project for it. Any modern C++11 compiler should work.
You will need / want:

* The open source C++ framework [DUNE](https://www.dune-project.org/),
  in particular its modules Dune-Grid and Dune-Istl [3].
* SuperLU as linear solver [4].
* [Paraview](www.paraview.org) for visualization and postprocessing.
* [Gmsh](http://gmsh.info/) for generating meshes [5].

## Comments on the examples chosen

For simplicity we considered simple geometries, two two-dimensional
rectangular plates, one resting over the other, and solved the equations
for isotropic materials with Lamé coefficients determined from Young's
modulus and Poisson's ratio. We chose Q1 elements in uniform meshes,
although given the facilities provided by Dune, it poses no problem to use
much more complicated meshes once any issues due to the realization of boundary
conditions are taken into account (for instance, the constraints may impose a
different choice for A\_0, as explained in [1, §6.2]).

## References

[1] S. Hüeber and B. I. Wohlmuth. “A primal-dual active set strategy for
    non-linear multibody contact problems.” Computer Methods in Applied
    Mechanics and Engineering, 194(27-29):3147-3166, jul 2005.

[2] B. I. Wohlmuth. “A Mortar Finite Element Method Using Dual Spaces for
    the Lagrange Multiplier.” SIAM Journal on Numerical Analysis,
    38(3):989-1012, jan 2001.

[3] M. Blatt and P. Bastian, “The iterative solver template library,”
    in Proceedings of the 8th international conference on applied parallel
    computing: state of the art in scientific computing, Berlin, Heidelberg,
    2007, pp. 666–675.

[4] X. S. Li, “An overview of SuperLU: Algorithms, implementation, and user
    interface,” ACM Transactions on Mathematical Software, vol. 31, no. 3,
    pp. 302–325, 2005.

[5] C. Geuzaine and J.-F. Remacle. “Gmsh: a three-dimensional finite element
    mesh generator with built-in pre- and post-processing facilities.”
    International Journal for Numerical Methods in Engineering 79(11),
    pp. 1309-1331, 2009.

[6] M. de Benito Delgado, “On the contact between two linearly elastic bodies,”
    Master’s thesis, Technische Universität München, Munich, 2013.

## License

This software falls under the GNU general public license version 3 or later.
It comes without **any warranty whatsoever**.
For details see http://www.gnu.org/licenses/gpl-3.0.html.
