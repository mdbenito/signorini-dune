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
solZ = pdi.GetPointData().GetArray('u[2]')
#sum  = pdi.GetPointData().GetArray('Sum') # in case we use the calculator
 
newPoints = vtk.vtkPoints()
for i in range(0, pdi.GetNumberOfPoints()):
  coord = pdi.GetPoint(i)
  x, y, z = coord[:3]
  x = x + solX.GetValue(i)
  y = y + solY.GetValue(i)
  z = z + solZ.GetValue(i)
  newPoints.InsertPoint(i, x, y, z)
 
pdo.SetPoints(newPoints)

 */

#include "config.h"             // (initially) constructed by ./configure script

  //// Dune includes

#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

  //// Project includes

#include "tests.hpp"
#include "benchmark.hpp"
#include "functors.hpp"
#include "shapefunctions.hpp"
#include "penaltymethod.hpp"
#include "activeset.hpp"
#include "twobodies.hpp"
#include "physicalgroupdescriptors.hpp"

  //// Ok, this sucks big time...

#define PLATE    0
#define PRISM    1
#define CYLINDER 2

static const std::string meshes_path ("/Users/miguel/Devel/Signorini/meshes/");
static const std::string problem_meshes[2][3] = {
  {"plate master.msh", "prism master.msh", "cylinder master.msh"},
  {"plate slave.msh",  "prism slave.msh",  "cylinder slave.msh"}
};

int main (int argc, char** argv)
{
  const int problem = PRISM;
  const int     dim = 3;
  const bool  tests = false;
  const double    E = 8.0e9;
  const double   E2 = 8.0e9;
  const double   nu = 0.3;
  const double  nu2 = 0.3;

    //// Grid setup

  typedef UGGrid<dim>              grid_t;
  typedef GridFactory<grid_t>   factory_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;
  typedef typename grid_t::LevelGridView LevelGV;
  typedef Codim1Extractor<LevelGV> SurfaceExtractor;
  typedef PSurfaceMerge<dim-1, dim, double> SurfaceMergeImpl;
  typedef ::GridGlue<SurfaceExtractor, SurfaceExtractor> GlueType;
  typedef PhysicalFaceDescriptor<LevelGV, factory_t> FaceDescriptor;

  factory_t factories[2];
  grid_t*       grids[2];
  GmshReader<grid_t> gmshReader[2];
  std::vector<int> boundary_id_to_physical_entity[2];
  std::vector<int> element_index_to_physical_entity[2];
  
    //// Problem data setup

  typedef HookeTensor<ctype, dim>                                   HookeT;
  typedef Constraint<coord_t, dim>                          ConstraintHack;
  typedef ConstantEvaluation<ctype, dim, coord_t>               VectorEval;
  typedef GmshVolumeFunctor<ctype, dim, factory_t, VectorEval>     VolumeF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval, ConstraintHack> BoundaryF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval, ConstraintHack> Dirichlet;
  
//  typedef CylinderHackGapEvaluation<ctype, dim>                    GapHack;
//  typedef PlateHackGapEvaluation<ctype, dim>                       GapHack;
//  typedef PrismHackGapEvaluation<ctype, dim>                       GapHack;
//  typedef GmshBoundaryFunctor<ctype, dim, factory_t, GapHack, ConstraintHack> Gap;
//  typedef BetterLinearShapeFunction<ctype, dim, 1, 0> ShapeF;
//  typedef MLinearShapeFunction<ctype, dim>         ShapeF;
//  typedef Q1ShapeFunctionSet<ctype, dim, ShapeF> ShapeSet;
//  typedef BetterLinearShapeFunction<ctype, dim, 3, -1> LSShapeF;
//  typedef LagrangeSpaceShapeFunction<ctype, dim>     LSShapeF;
//  typedef Q1ShapeFunctionSet<ctype, dim, LSShapeF> LSShapeSet;

  typedef P1ShapeFunctionSet<ctype, dim, 1, 0> ShapeSet;
  typedef P1ShapeFunctionSet<ctype, dim, 3, -1> LSShapeSet;
  
//  typedef SignoriniFEPenalty<GV, HookeT, VolumeF, BoundaryF, Gap, Dirichlet, ShapeSet>
//          PMSolver;
//  typedef SignoriniIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, Gap, ShapeSet, LSShapeSet>
//          IASolver;
  typedef TwoBodiesIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, GlueType, ShapeSet, LSShapeSet>
          TwoSolver;
  
    // TODO: patch the GmshParser to read Physical Group names and use those
  
  std::set<int> volumeGroups[2], contactGroups[2], dirichletGroups[2], neumannGroups[2];
  switch (problem) {
    case PLATE:
      volumeGroups[MASTER] << 1;
      volumeGroups[SLAVE]  << 1;
      contactGroups[MASTER] << 3;  // plate
      contactGroups[SLAVE]  << 1;  // plate
      dirichletGroups[MASTER] << 1;  // plate
      dirichletGroups[SLAVE]  << 3;  // plate
      neumannGroups[MASTER] << 2 << 4;  // plate
      neumannGroups[SLAVE]  << 2 << 4;  // plate
      break;
    case PRISM:
      volumeGroups[MASTER] << 1;
      volumeGroups[SLAVE]  << 1;
      contactGroups[MASTER] << 1;       // prism
      contactGroups[SLAVE]  << 1;       // prism
      dirichletGroups[MASTER] << 6;       // prism
      dirichletGroups[SLAVE]  << 6;       // prism
      neumannGroups[MASTER] << 2 << 3 << 4 << 5;  // prism
      neumannGroups[SLAVE]  << 2 << 3 << 4 << 5;  // prism
      break;
    case CYLINDER:
      volumeGroups[MASTER] << 1;
      volumeGroups[SLAVE]  << 1;
      contactGroups[MASTER] << 3 << 4;  // cylinder
      contactGroups[SLAVE]  << 5 << 6;  // cylinder
      dirichletGroups[MASTER] << 1 << 2;  // cylinder
      dirichletGroups[SLAVE]  << 1 << 2;  // cylinder
      neumannGroups[MASTER] << 5 << 6;  // cylinder
      neumannGroups[SLAVE]  << 3 << 4;  // cylinder
      break;
    default:
      assert ("huh?" && false);
  }

  HookeT a (E, nu);
  HookeT a2 (E2, nu2);

  VolumeF f (factories[MASTER],
             element_index_to_physical_entity[MASTER],
             volumeGroups[MASTER],
             new VectorEval (coord3 (0.0, 4e8, 0.0)));
//              new VectorEval (coord2 (0.0, 4e8)));
  VolumeF f2 (factories[SLAVE],
              element_index_to_physical_entity[SLAVE],
              volumeGroups[SLAVE],
              new VectorEval (coord3 (0.0, -4e8, 0.0)));
//              new VectorEval (coord2 (0.0, -4e8)));
  Dirichlet d (factories[MASTER],
               boundary_id_to_physical_entity[MASTER],
               dirichletGroups[MASTER],
               new VectorEval (coord3 (0.0, 0.01, 0.0)));
//                new VectorEval (coord2 (0.0, 0.01)));
  Dirichlet d2 (factories[SLAVE],
                boundary_id_to_physical_entity[SLAVE],
                dirichletGroups[SLAVE],
                new VectorEval (coord3 (0.0, -0.01, 0.0)));
//                new VectorEval (coord2 (0.0, -0.01)));
  BoundaryF p (factories[MASTER],
               boundary_id_to_physical_entity[MASTER],
               neumannGroups[MASTER],
               new VectorEval (coord3 (3e6, 0.0, 0.0)));
//                new VectorEval (coord2 (3e6, 0.0)));
  BoundaryF p2 (factories[SLAVE],
                boundary_id_to_physical_entity[SLAVE],
                neumannGroups[SLAVE],
                new VectorEval (coord3 (0.0, -3e6, 0.0)));
//                new VectorEval (coord2 (0.0, -3e6)));

//  Gap g (factories[MASTER],
//         boundary_id_to_physical_entity[MASTER],
//         contactGroups[MASTER],
//         new GapHack ());
//  Gap g2 (factories[SLAVE],
//          boundary_id_to_physical_entity[SLAVE],
//          contactGroups[SLAVE],
//          new GapHack ());
//

    //// Setup grids
  
  DOBOTH(body) {
    gmshReader[body].read (factories[body],
                           meshes_path + problem_meshes[body][problem],
                           boundary_id_to_physical_entity[body],
                           element_index_to_physical_entity[body]);
    grids[body] = factories[body].createGrid();
  }

    //// Setup contact surfaces
  
  SurfaceMergeImpl merger;
  FaceDescriptor descriptors[2] = {
    FaceDescriptor (factories[MASTER], boundary_id_to_physical_entity[MASTER],contactGroups[MASTER]),
    FaceDescriptor (factories[SLAVE], boundary_id_to_physical_entity[SLAVE], contactGroups[SLAVE]) };
  SurfaceExtractor extractors[2] = {
    SurfaceExtractor(grids[MASTER]->levelView (0), descriptors[MASTER]),
    SurfaceExtractor(grids[SLAVE]->levelView (0), descriptors[SLAVE]) };
  
  GlueType glue (extractors[MASTER], extractors[SLAVE], &merger);     // FIXME: I've switched the names/roles!!

  glue.build();
  assert (glue.size() > 0);
  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  
  if (tests) {
    testShapes<ctype, dim, ShapeSet> ("/tmp/basis");
    testShapes<ctype, dim, LSShapeSet> ("/tmp/multbasis");
    test_basisOfPlaneNormalTo<ctype, dim>();
    
      //  TwoRefs<Gap> gaps (g, g2);
    TwoRefs<BoundaryF> neumann (p, p2);
    TwoRefs<Dirichlet> dirichlet (d, d2);
    TwoRefs<GV> gv (grids[MASTER]->leafView(), grids[SLAVE]->leafView());
    
    DOBOTH (body) {
        //    testGmshBoundaryFunctor (gv[body], gaps[body], string ("/tmp/test-gap-") + body);
      testGmshBoundaryFunctor (gv[body], neumann[body], string ("/tmp/test-neumann-") + body);
      testGmshBoundaryFunctor (gv[body], dirichlet[body], string ("/tmp/test-dirichlet-") + body);
    }
    
    testContactSurfaces<GlueType, MASTER> (glue, "/tmp/test-gap");
    testContactSurfaces<GlueType, SLAVE> (glue, "/tmp/test-gap");
    exit (1);
  }

//  PMSolver fem (grids[MASTER]->leafView(), a, f, p, g, d, 1.0e-14 / E);
//  IASolver fem (grids[MASTER]->leafView(), a, f, d, p, g);
  TwoSolver fem (grids[MASTER]->leafView(), grids[SLAVE]->leafView(),
                 a, a2, f, f2, d, d2, p, p2, glue, 3);

    //// Solution
  
  try {   // Pokemon Exception Handling
//    cout << "Gremlin population: " << grids[MASTER]->size(dim) + grids[SLAVE]->size(dim) << LF;
    
//      // Penalty method, one body
//    fem.initialize();
//    fem.solve (10, 1e-6);  // args: maxsteps, tolerance
//      // Active / inactive
    fem.solve();
    return 0;
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
