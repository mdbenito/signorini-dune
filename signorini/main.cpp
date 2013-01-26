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

  // Dune includes
#include "config.h"                    // file constructed by ./configure script
#include <dune/grid/sgrid.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/gmshreader.hh>
  //#include <dune/grid/alugrid.hh>

#include "traverse.hpp"
#include "functors.hpp"
#include "integration.hpp"
#include "finiteelements.cpp"
#include "benchmark.hpp"
#include "postprocessor.cpp"

#include <sstream>

int main (int argc, char** argv)
{
  static const int dim = 2;

    //// Grid creation (SGrid, test)

  typedef SGrid<dim, dim>          grid_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;

  FieldVector<int, dim> N(1);
  coord_t       origin (0.0);
  coord_t     topright (1.0);
  
  grid_t grid (N, origin, topright);
  grid.globalRefine (7);
  
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

  typedef VolumeLoad<ctype, dim>     VolumeForces;
  typedef Tractions<ctype, dim>    BoundaryForces;
  typedef NormalGap<ctype, dim>               Gap;
  typedef HookeTensor<ctype, dim>           Hooke;
  typedef Q1ShapeFunctionSet<ctype, dim> ShapeSet;
  
  Hooke          a (5e9, 0.3);  // E = 5x10^9 Pa, nu = 0,3 (cf. [FV05, p.35])
  VolumeForces   f;
  BoundaryForces p;
  Gap            g;
  
  const GV& gv = grid.leafView();
  
  
  SignoriniFEPenalty<GV, Hooke, VolumeForces, BoundaryForces, Gap> p1 (gv, a, f, p, g, 4);
  
    //// 
  
  Benchmark bench;
  PostProcessor<GV, Hooke, ShapeSet> post (gv, a);
  
    //// Solution
  
  try {   // Pokemon Exception Handling!!
    
      //testShapes<ctype, dim, Q1ShapeFunctionSet>();

    cout << "Gremlin population: " << grid.size(dim) << "\n";

    bench.start ("Initialization");
    p1.initialize();
    p1.assembleMain();
    bench.stop ("Initialization");
    
      //printmatrix(std::cout, p1.A, "Stiffness matrix","");
  
    int         step = 0;
    double       eps = 1.0e-3 / 5.0e9;  // See [KO88, p.140]
    double tolerance = 1.0e-5;
    double     error = 1.0;
    
    bench.start ("Resolution");
    while (step++ < 10 && (error > tolerance || (error <= tolerance && step < 3))) {
        //while (step++ < 10 && error > tolerance) {
      bench.start ("Iteration");
      
      cout << "Iteration " << step << ":\n";

      if (1) {
        bench.start ("Penalty matrix assembly");
        p1.assemblePenalties (eps);
        //printmatrix(std::cout, p1.P, "Penalty matrix","");
        bench.stop ("Penalty matrix assembly");
      }
      
      bench.start ("System resolution");
      p1.solve();
      bench.stop ("System resolution");
      
      bench.start ("Postprocessing");
      post.check (p1.solution());
      error = post.computeError (p1.solution ());
      post.computeVonMises();
      cout << "\t\tNew solution diverged by: " << error << "\n";
      cout << "\t\tOutputting... ";
      post.writeVTKFile ("/tmp/SignoriniFEM", step);
      cout << "Done.\n";
      bench.stop ("Postprocessing");
      
      bench.stop ("Iteration");
    }
    bench.stop ("Resolution");
    return 0;
  } catch (MatrixBlockError& e) {
    cout << "DEAD! " << e.what() << "\n";// << p1.A[e.r][e.c] << "\n";
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
