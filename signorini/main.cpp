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
//#include <dune/grid/alugrid.hh>

  //// Project includes

#include "benchmark.hpp"
#include "functors.hpp"
#include "shapefunctions.hpp"
#include "penaltymethod.hpp"
#include "activeset.hpp"

  //// stdlib includes

#include <string>
using std::string;

int main (int argc, char** argv)
{
  const int          dim = 2;
    //const double         E = 5.0e9;       // Young's modulus (in Pa) [FV05, p.35]
    //const double        nu = 0.3;         // Poisson's ratio         [FV05, p.35]
  const double         E = 200;         // See [HW04, p.3154]
  const double        nu = 0.3;         // See [HW04, p.3154]
  const double       eps = 1.0e-14 / E;  // See [KO88, p.140]
  const double tolerance = 1.0e-5;      // For the iterative penalty method
  const int     maxsteps = 10;

    //// Grid setup (SGrid, test)

  typedef SGrid<dim, dim>          grid_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;

  FieldVector<int, dim> N(1);
  coord_t       origin (0.0);
  coord_t     topright (1.0);
  
  grid_t grid (N, origin, topright);
  grid.globalRefine (3);

  const GV& gv = grid.leafView();

  /* This won't work... (AluGrid seems not to be properly configured)

  typedef AluGrid<3, 3> grid_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;
  
  grid_t grid;
  std::string gridName = "/Users/miguel/cube-fine.msh";
  GmshReader<grid_t> gmshReader;
  gmshReader.read (gridName);
  
  gridName.erase (0 , gridName . rfind ("/")+1);
  gridName.erase(gridName.find(".",0), gridName.length());
  */

    //// Problem setup

  typedef VolumeLoad<ctype, dim>          VolumeF;
  typedef Tractions<ctype, dim>         BoundaryF;
  typedef NormalGap<ctype, dim>               Gap;
  typedef HookeTensor<ctype, dim>          HookeT;
  typedef MLinearShapeFunction<ctype, dim>             ShapeF;
  typedef Q1ShapeFunctionSet<ctype, dim, ShapeF>     ShapeSet;
  typedef LagrangeSpaceShapeFunction<ctype, dim>     LSShapeF;
  typedef Q1ShapeFunctionSet<ctype, dim, LSShapeF> LSShapeSet;
  
  typedef SignoriniFEPenalty<GV, HookeT, VolumeF, BoundaryF, Gap, ShapeSet>
          PMSolver;
  typedef SignoriniIASet<GV, HookeT, VolumeF, BoundaryF, Gap, ShapeSet, LSShapeSet>
          IASolver;
    //Should be:
    //SignoriniIASet<GV, ProblemData, ShapeSet> IASolver;
  
    //testShapes<ctype, dim, LSShapeSet>();
    //testShapes<ctype, dim, ShapeSet>();
  
  HookeT  a (E, nu);
  VolumeF         f;
  BoundaryF       p;
  Gap             g;
  
    //PMSolver  fem (gv, a, f, p, g, eps);
  IASolver fem2 (gv, a, f, p, g);

    //// Misc.

    //PostProcessor<GV, ProblemData, ShapeSet> post (gv, data);

  
    //// Solution
  
  try {   // Pokemon Exception Handling!!
    cout << "Gremlin population: " << grid.size(dim) << "\n";
      //fem.initialize();
      //fem.solve (maxsteps, tolerance);
    fem2.initialize ();
    fem2.solve ();
    return 0;
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
