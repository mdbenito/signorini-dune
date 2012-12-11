/******************************************************************************
 * main.cpp                                                                   *
 
 Python script for visualization of displacements in ParaView:
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
 x = x + solX.GetValue(i)/5
 y = y + solY.GetValue(i)/5

 newPoints.InsertPoint(i, x, y, z)
 
 pdo.SetPoints(newPoints)
 
 ******************************************************************************/

  // Dune includes
#include "config.h"                    // file constructed by ./configure script
#include <dune/grid/sgrid.hh>
#include <dune/common/mpihelper.hh>    // include mpi helper class
#include <dune/common/exceptions.hh>   // We use exceptions
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
  //#include <dune/grid/alugrid.hh>

#include "traverse.hpp"
#include "functors.hpp"
#include "integration.hpp"
#include "finiteelements.cpp"

#include <ctime>
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

  typedef VolumeLoad<ctype, dim>   VolumeForces;
  typedef Tractions<ctype, dim>  BoundaryForces;
  typedef NormalGap<ctype, dim>             Gap;
  typedef HookeTensor<ctype, dim>         Hooke;
  
  Hooke a (5e9, 0.3);  // E = 5x10^9 Pa, nu = 0,3 (rock material)
  VolumeForces f;
  BoundaryForces p;
  Gap g;
  
  const GV& gv = grid.leafView();
  
  SignoriniFEPenalty<GV, Hooke, VolumeForces, BoundaryForces, Gap> p1 (gv, a, f, p, g);
  
    //// Solution
  
  try {   // Pokemon Exception Handling!!
    
      //testShapes<ctype, dim, Q1ShapeFunctionSet>();
    
    cout << "-----------------------------------\n";
    cout << "Number of unknowns: " << grid.size(dim) << "\n";

    p1.initialize();

    p1.assembleMain();
      //printmatrix(std::cout, p1.A, "Stiffness matrix","");
  
    int     step = 0;
    double   eps = 1.0e-3;
    double error = 1.0;
    
    while (step++ < 10 && error > 1.0e-3) {
      auto previous = p1.solutionAsVector();
      
      cout << "Assembling data for iteration " << step << "... ";
        //p1.assemblePenalties (eps);
        //printmatrix(std::cout, p1.P, "Penalty matrix","");
      
      cout << "Done.\n" << "Solving... "; std::flush(cout);  // Yuk!
      clock_t t = std::clock();
      p1.solve();

      bool ok = true;
      for (const auto& it : p1.solution())
        for (int i = 0; i < dim; ++i)
          ok = ok && (it[i] != NAN);
      
      if (!ok) {
        cout << "NaN!\n" << "What a calamity, I choose to quit.\n";
        exit (1);
      }

        // poor man's benchmarking
      double lapse = (1.0*std::clock() - 1.0*t)/(1.0*CLOCKS_PER_SEC);
      
      auto oldPrecision = cout.precision();
      cout << " done in " << std::setprecision(3) << lapse << " seconds.\n";
      cout.precision(oldPrecision);
      
      error = computeError (p1.solutionAsVector(), previous);
      cout << "New solution diverged by: " << error << "\n";

      cout << "Outputting... ";
      std::ostringstream oss;
      oss << "/tmp/SignoriniFEM" << dim << "d-"
          << std::setfill('0') << std::setw(3) << step;
      VTKWriter<grid_t::LeafGridView> vtkwriter (grid.leafView());
      vtkwriter.addVertexData (p1.solutionAsVector(), "u", dim);
      vtkwriter.write (oss.str(), VTKOptions::binaryappended);
      cout << "Done.\n";
    }
    
    return 0;
  } catch (MatrixBlockError& e) {
    cout << "DEAD! " << e.what() << "\n";// << p1.A[e.r][e.c] << "\n";
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
