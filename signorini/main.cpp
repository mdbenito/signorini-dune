/******************************************************************************
 * main.cpp                                                                   *
 ******************************************************************************/

/*
 Here is a Python script for visualization of displacements in ParaView:
 
from paraview import vtk
 
pdi  = self.GetInput()
pdo  = self.GetOutput()
solX = pdi.GetPointData().GetArray('u[0]')
solY = pdi.GetPointData().GetArray('u[1]')
#sum  = pdi.GetPointData().GetArray('Sum') # in case we use the calculator
 
newPoints = vtk.vtkPoints()
for i in range(0, pdi.GetNumberOfPoints()):
  coord = pdi.GetPoint(i)
  x, y, z = coord[:3]
  x = x + solX.GetValue(i)
  y = y + solY.GetValue(i)
  newPoints.InsertPoint(i, x, y, z)
 
pdo.SetPoints(newPoints)

 */

  //// Dune includes

#include "config.h"                    // file constructed by ./configure script
#include <dune/grid/sgrid.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/alugrid.hh>

  //// Project includes

#include "benchmark.hpp"
#include "functors.hpp"
#include "shapefunctions.hpp"
#include "penaltymethod.hpp"
#include "activeset.hpp"
#include "twobodies.hpp"

  //// stdlib includes

//#include <string>
//using std::string;

int main (int argc, char** argv)
{
  const int          dim = 2;

  const double         E = 8.0e10;
  const double        E2 = 8.0e10;
  const double        nu = 0.3;
  const double       nu2 = 0.3;
  
  //const double         E = 5.0e9;       // Young's modulus (in Pa) [FV05, p.35]
  //const double       eps = 1.0e-5 / E;  // See [KO88, p.140]
  //const double       eps = 1.0e-14 / E;  // See [KO88, p.140]
  //const double tolerance = 1.0e-5;      // For the iterative penalty method
  //const int     maxsteps = 10;

    //// Grid setup

  typedef SGrid<dim, dim>          grid_t;
//  typedef ALUSimplexGrid<dim, dim> grid_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;

  FieldVector<int, dim> N(1);
  coord_t       origin (0.0);
  coord_t     topright (1.0);
  
  grid_t gridSlave (N, origin, topright);
  gridSlave.globalRefine (3);
  
    // Careful changing this! Check the supports of the boundary functors!!
  origin[1] = -1;
  topright[1] = 0;
  grid_t gridMaster (N, origin, topright);
  gridMaster.globalRefine (3);

  /*
  std::string path ("/Users/miguel/Universidad/Two bodies/");
  grid_t gridMaster, gridSlave;
  GmshReader<grid_t> gmshReader;
  
  gmshReader.read (gridMaster, path * "cylinder master.msh");
  gmshReader.read (gridSlave, path * "cylinder slave.msh");
  */
    //// Problem setup

  typedef VolumeLoad<ctype, dim>          VolumeF;
  typedef Tractions<ctype, dim>         BoundaryF;
  typedef NormalGap<ctype, dim>               Gap;
  typedef DirichletFunctor<ctype, dim>  Dirichlet;
  typedef HookeTensor<ctype, dim>          HookeT;

  typedef MLinearShapeFunction<ctype, dim>             ShapeF;
  typedef Q1ShapeFunctionSet<ctype, dim, ShapeF>     ShapeSet;
  typedef LagrangeSpaceShapeFunction<ctype, dim>     LSShapeF;
  typedef Q1ShapeFunctionSet<ctype, dim, LSShapeF> LSShapeSet;
  
  typedef SignoriniFEPenalty<GV, HookeT, VolumeF, BoundaryF, Gap, Dirichlet, ShapeSet>
          PMSolver;
  typedef SignoriniIASet<GV, HookeT, VolumeF, BoundaryF, Gap, ShapeSet, LSShapeSet>
          IASolver;
  typedef TwoBodiesIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, Gap, ShapeSet, LSShapeSet>
          TwoSolver;

    //Should be:
    //SignoriniIASet<GV, ProblemData, ShapeSet> IASolver;
  
//  testShapes<ctype, dim, ShapeSet>("/tmp/basis");
//  testShapes<ctype, dim, LSShapeSet>("/tmp/multbasis");
//  test_basisOfPlaneNormalTo<ctype, dim>();
//  exit(1);

  HookeT    a (E, nu);
  HookeT    a2 (E2, nu2);
//  VolumeF   f (coord3 (0.0, 0.0, 0.0));
//  VolumeF   f2 (coord3 (0.0, 0.0, 0.0));
//  Dirichlet d (coord3 (0.0, -0.07, 0.0));
//  Dirichlet d2 (coord3 (0.0, 0.0, 0.0));
//  BoundaryF p (coord3 (-3e6, -7e6, 0.0));
//  BoundaryF p2 (coord3 (3e6, -4e6, 0.0));
////  Gap       g (0.0, 0.0);
////  Gap       g2 (0.0, 0.0);
//  Gap g (0.0, 0.05);
  VolumeF   f (coord2 (0.0, 0.0));
  VolumeF   f2 (coord2 (0.0, 0.0));
  Dirichlet d (coord2 (0.0, -0.07));
  Dirichlet d2 (coord2 (0.0, 0.0));
  BoundaryF p (coord2 (-3e6, -7e6));
  BoundaryF p2 (coord2 (3e6, -4e6));
    //  Gap       g (0.0, 0.0);
    //  Gap       g2 (0.0, 0.0);
  Gap g (0.0, 0.05);
  
  
    //PMSolver  fem (gridSlave.leafView(), a, f, p, g, d, eps);
  IASolver fem2 (gridSlave.leafView(), a, f, p, g);
//  TwoSolver twoFem (gridMaster.leafView(), gridSlave.leafView(),
//                    a, a2, f, f2, d, d2, p, p2, g, g2, 4);

    //// Solution
  
  try {   // Pokemon Exception Handling
    cout << "Gremlin population: " << gridSlave.size(dim) << "\n";
//      // Penalty method, one body
//      fem.initialize();
//      fem.solve (maxsteps, tolerance);
//    
      // Active / inactive, one body
      fem2.solve ();
//      // Active / inactive, two bodies
//      twoFem.solve();
    return 0;
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
