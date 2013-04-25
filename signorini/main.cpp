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

#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

  //// Project includes

#include "benchmark.hpp"
#include "functors.hpp"
#include "shapefunctions.hpp"
#include "penaltymethod.hpp"
#include "activeset.hpp"
#include "twobodies.hpp"
#include "physicalgroupdescriptors.hpp"

  //// stdlib includes

//#include <string>
//using std::string;

int main (int argc, char** argv)
{
  const int          dim = 3;

  const double         E = 8.0e9;
  const double        E2 = 8.0e9;
  const double        nu = 0.3;
  const double       nu2 = 0.3;
  
  //const double         E = 5.0e9;       // Young's modulus (in Pa) [FV05, p.35]
  //const double       eps = 1.0e-5 / E;  // See [KO88, p.140]
  //const double       eps = 1.0e-14 / E;  // See [KO88, p.140]
  //const double tolerance = 1.0e-5;      // For the iterative penalty method
  //const int     maxsteps = 10;

    //// Grid setup

//  typedef SGrid<dim, dim>          grid_t;
  typedef ALUCubeGrid<dim, dim>    grid_t;
  typedef GridFactory<grid_t>   factory_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;

  FieldVector<int, dim> N(1);
  coord_t       origin (0.0);
  coord_t     topright (1.0);
  /*
  grid_t gridSlave (N, origin, topright);
  gridSlave.globalRefine (3);
  
    // Careful changing this! Check the supports of the boundary functors!!
  origin[1] = -1;
  topright[1] = 0;
  grid_t gridMaster (N, origin, topright);
  gridMaster.globalRefine (3);
*/
  
  string path ("/Users/miguel/Devel/Signorini/meshes/");
  factory_t factories[2];
  grid_t*       grids[2];
  GmshReader<grid_t> gmshReader;
  std::vector<int> boundary_id_to_physical_entity[2];
  std::vector<int> element_index_to_physical_entity[2];
  gmshReader.read (factories[MASTER],
                   path + string("cylinder master.msh"),
                   boundary_id_to_physical_entity[MASTER],
                   element_index_to_physical_entity[MASTER]);
  gmshReader.read (factories[SLAVE],
                   path + string("cylinder slave.msh"),
                   boundary_id_to_physical_entity[SLAVE],
                   element_index_to_physical_entity[SLAVE]);

  grids[MASTER] = factories[MASTER].createGrid();
  grids[SLAVE]  = factories[SLAVE].createGrid();
  
    //// Setup contact surfaces
  
  typedef typename grid_t::LevelGridView LevelGV;
  typedef Codim1Extractor<LevelGV> SurfaceExtractor;
  typedef PSurfaceMerge<dim-1, dim, double> SurfaceMergeImpl;
  typedef ::GridGlue<SurfaceExtractor, SurfaceExtractor> GlueType;
  
  std::set<int> groups[2];
  groups[MASTER] << 3 << 4;  // FIXME: what are the right identifiers?
  groups[SLAVE]  << 5 << 6;  // FIXME: what are the right identifiers?

  PhysicalFaceDescriptor<LevelGV, factory_t> masterDescriptor (factories[MASTER],
                                                               boundary_id_to_physical_entity[MASTER],
                                                               groups[MASTER]);
  PhysicalFaceDescriptor<LevelGV, factory_t> slaveDescriptor (factories[SLAVE],
                                                              boundary_id_to_physical_entity[SLAVE],
                                                              groups[SLAVE]);

  SurfaceExtractor masterExtractor (grids[MASTER]->levelView(0), masterDescriptor);
  SurfaceExtractor slaveExtractor (grids[SLAVE]->levelView(0), slaveDescriptor);
  
  SurfaceMergeImpl merger;
  GlueType glue (slaveExtractor, masterExtractor, &merger);   // FIXME: careful with the order
  
  glue.build();

  assert (glue.size() > 0);
  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  
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
  VolumeF   f (coord3 (0.0, 0.0, 0.0));
  VolumeF   f2 (coord3 (0.0, 0.0, 0.0));
  Dirichlet d (coord3 (0.0, -0.07, 0.0));
  Dirichlet d2 (coord3 (0.0, 0.0, 0.0));
  BoundaryF p (coord3 (-3e6, -7e6, 0.0));
  BoundaryF p2 (coord3 (3e6, -4e6, 0.0));
//  Gap       g (0.0, 0.0);
//  Gap       g2 (0.0, 0.0);
  Gap g (0.0, 0.05);

//  VolumeF   f (coord2 (0.0, 0.0));
//  VolumeF   f2 (coord2 (0.0, 0.0));
//  Dirichlet d (coord2 (0.0, -0.07));
//  Dirichlet d2 (coord2 (0.0, 0.0));
//  BoundaryF p (coord2 (-3e6, -7e6));
//  BoundaryF p2 (coord2 (3e6, -4e6));
//    //  Gap       g (0.0, 0.0);
//    //  Gap       g2 (0.0, 0.0);
//  Gap g (0.0, 0.05);
  
  
    //PMSolver  fem (gridSlave.leafView(), a, f, p, g, d, eps);
  IASolver fem2 (grids[SLAVE]->leafView(), a, f, p, g);
//  TwoSolver twoFem (gridMaster.leafView(), gridSlave.leafView(),
//                    a, a2, f, f2, d, d2, p, p2, g, g2, 4);

    //// Solution
  
  try {   // Pokemon Exception Handling
    cout << "Gremlin population: " << grids[SLAVE]->size(dim) << "\n";
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
