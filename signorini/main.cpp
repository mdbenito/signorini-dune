/******************************************************************************
 * main.cpp                                                                   *
 * TODO (among many other things)                                             *
 *   - patch the GmshParser to read Physical Group names and use those        *
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

typedef enum { PLATE=0, PRISM=1, CYLINDER=2 } ProblemType;
static const std::string meshPath ("/Users/miguel/Devel/Signorini/meshes/");
static const std::string meshNames[2][3] = {
  {"plate master.msh", "prism master.msh", "cylinder master.msh"},
  {"plate slave.msh",  "prism slave.msh",  "cylinder slave.msh"}
};

int main (int argc, char** argv)
{
  ProblemType problem = PRISM;
  const int       dim = 3;
  const bool    tests = true;
  const double   E[2] = { 8.0e9, 8.0e9 };
  const double  nu[2] = { 0.3,     0.3 };

  typedef UGGrid<dim>              grid_t;
  typedef GridFactory<grid_t>   factory_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef grid_t::LeafGridView         GV;
  typedef Codim1Extractor<GV>                    SurfaceExtractor;
  typedef PSurfaceMerge<dim-1, dim, double>      SurfaceMergeImpl;
  typedef ::GridGlue<SurfaceExtractor, SurfaceExtractor> GlueType;
  typedef PhysicalFaceDescriptor<GV, factory_t>    FaceDescriptor;
  typedef HookeTensor<ctype, dim>                          HookeT;
  typedef Constraint<coord_t, dim>                 ConstraintHack;
  typedef ConstantEvaluation<ctype, dim, coord_t>      VectorEval;
  typedef GmshVolumeFunctor<ctype, dim, factory_t, VectorEval>                     VolumeF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval, ConstraintHack> BoundaryF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval, ConstraintHack> Dirichlet;
  typedef P1ShapeFunctionSet<ctype, dim, 1, 0> ShapeSet;
  typedef P1ShapeFunctionSet<ctype, dim, 3, -1> LSShapeSet;
  typedef TwoBodiesIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, GlueType, ShapeSet, LSShapeSet>
          TwoSolver;
  /*
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
   //  typedef SignoriniFEPenalty<GV, HookeT, VolumeF, BoundaryF, Gap, Dirichlet, ShapeSet>
   //          PMSolver;
   //  typedef SignoriniIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, Gap, ShapeSet, LSShapeSet>
   //          IASolver;
   */

  factory_t        factories[2];
  grid_t*              grids[2];
  GmshReader<grid_t> readers[2];
  std::vector<int>     bi2pe[2]; // Boundary ids to gmsh physical entities
  std::vector<int>     ei2pe[2]; // Element indices to gmsh physical entities

  std::set<int> volumeGroups[2], contactGroups[2], dirichletGroups[2], neumannGroups[2];
  HookeT* a[2]; VolumeF* f[2];  Dirichlet* d[2];  BoundaryF* p[2]; //Gap* g[2];
  VectorEval* fEval[2], *dEval[2], *pEval[2];
  FaceDescriptor*  descriptors[2];
  SurfaceExtractor* extractors[2];

  switch (problem) {
    case PLATE:
      volumeGroups   [MASTER] << 1;       volumeGroups   [SLAVE] << 1;
      contactGroups  [MASTER] << 3;       contactGroups  [SLAVE] << 1;
      dirichletGroups[MASTER] << 1;       dirichletGroups[SLAVE] << 3;
      neumannGroups  [MASTER] << 2 << 4;  neumannGroups  [SLAVE] << 2 << 4;
/*      fEval[MASTER] = new VectorEval (coord2 (0.0,  4e8));
      fEval[SLAVE]  = new VectorEval (coord2 (0.0, -4e8));
      dEval[MASTER] = new VectorEval (coord2 (0.0, 0.01));
      dEval[SLAVE]  = new VectorEval (coord2 (0.0, -0.01));
      pEval[MASTER] = new VectorEval (coord2 (3e6, 0.0));
      pEval[SLAVE]  = new VectorEval (coord2 (0.0, -3e6));*/
      break;
    case PRISM:
      volumeGroups   [MASTER] << 1;      volumeGroups   [SLAVE] << 1;
      contactGroups  [MASTER] << 1;      contactGroups  [SLAVE] << 1;
      dirichletGroups[MASTER] << 6;      dirichletGroups[SLAVE] << 6;
      neumannGroups  [MASTER] << 2 << 3 << 4 << 5;
      neumannGroups  [SLAVE]  << 2 << 3 << 4 << 5;
      fEval[MASTER] = new VectorEval (coord3 (0.0,  4e8, 0.0));
      fEval[SLAVE]  = new VectorEval (coord3 (0.0, -4e8, 0.0));
      dEval[MASTER] = new VectorEval (coord3 (0.0, 0.01, 0.0));
      dEval[SLAVE]  = new VectorEval (coord3 (0.0, -0.01, 0.0));
      pEval[MASTER] = new VectorEval (coord3 (3e6, 0.0, 0.0));
      pEval[SLAVE]  = new VectorEval (coord3 (0.0, -3e6, 0.0));
      break;
    case CYLINDER:
      volumeGroups   [MASTER] << 1;       volumeGroups   [SLAVE] << 1;
      contactGroups  [MASTER] << 3 << 4;  contactGroups  [SLAVE] << 5 << 6;
      dirichletGroups[MASTER] << 1 << 2;  dirichletGroups[SLAVE] << 1 << 2;
      neumannGroups  [MASTER] << 5 << 6;  neumannGroups  [SLAVE] << 3 << 4;
      fEval[MASTER] = new VectorEval (coord3 (0.0,  4e8, 0.0));
      fEval[SLAVE]  = new VectorEval (coord3 (0.0, -4e8, 0.0));
      dEval[MASTER] = new VectorEval (coord3 (0.0, 0.01, 0.0));
      dEval[SLAVE]  = new VectorEval (coord3 (0.0, -0.01, 0.0));
      pEval[MASTER] = new VectorEval (coord3 (3e6, 0.0, 0.0));
      pEval[SLAVE]  = new VectorEval (coord3 (0.0, -3e6, 0.0));
      break;
    default:
      std::cerr << "huh?" << LF;
      exit (42);
  }

  DOBOTH (body) {
    readers[body].read (factories[body], meshPath + meshNames[body][problem], bi2pe[body], ei2pe[body]);
    grids[body] = factories[body].createGrid();

    descriptors[body] = new FaceDescriptor (factories[body], bi2pe[body], contactGroups[body]);
    
      ///FIXME!!! instantiating here the extractors[] results in both referring to the same GridView:
      // the glue has then both views pointing to the master.
//    extractors[body]  = new SurfaceExtractor (grids[body]->leafView(), *descriptors[body]);
    
    a[body] = new HookeT (E[body], nu[body]);
    f[body] = new VolumeF (factories[body], ei2pe[body], volumeGroups[body], fEval[body]);
    d[body] = new Dirichlet (factories[body], bi2pe[body], dirichletGroups[body], dEval[body]);
    p[body] = new BoundaryF (factories[body], bi2pe[body], neumannGroups[body], pEval[body]);
//    g[body] = new Gap (factories[body], bi2pe[body], contactGroups[body], new GapHack ());
  }
    // HACK: see above why we don't do this in the loop.
  extractors[MASTER]  = new SurfaceExtractor (grids[MASTER]->leafView(), *descriptors[MASTER]);
  extractors[SLAVE]  = new SurfaceExtractor (grids[SLAVE]->leafView(), *descriptors[SLAVE]);
  
  SurfaceMergeImpl merger;
  GlueType glue (*extractors[MASTER], *extractors[SLAVE], &merger); // FIXME: I've switched the names/roles!!
  glue.build ();
  assert (glue.size() > 0);
  
  TwoRefs<GV> gv (grids[MASTER]->leafView(), grids[SLAVE]->leafView());

  if (tests) {
    testShapes<ctype, dim, ShapeSet> ("/tmp/basis");
    testShapes<ctype, dim, LSShapeSet> ("/tmp/multbasis");
    test_basisOfPlaneNormalTo<ctype, dim>();

    DOBOTH (body) {
      testGmshBoundaryFunctor (gv[body], *p[body], string ("/tmp/test-neumann-") + body);
      testGmshBoundaryFunctor (gv[body], *d[body], string ("/tmp/test-dirichlet-") + body);
//    testGmshBoundaryFunctor (gv[body], *g[body], string ("/tmp/test-gap-") + body);
    }
    testContactSurfaces<GlueType, MASTER> (glue, "/tmp/test-gap");
    testContactSurfaces<GlueType, SLAVE> (glue, "/tmp/test-gap");
//    exit (1);
  }

//  PMSolver fem (gv[MASTER], a[MASTER], f[MASTER], p[MASTER], g[MASTER], d[MASTER], 1.0e-14 / E, 1e-6, 10);
//  IASolver fem (gv[MASTER], a[MASTER], f[MASTER], d[MASTER], p[MASTER], g[MASTER]);
  TwoSolver fem (gv[MASTER], gv[SLAVE], *a[MASTER], *a[SLAVE], *f[MASTER], *f[SLAVE], *d[MASTER], *d[SLAVE], *p[MASTER], *p[SLAVE], glue, 3);
  
  try {   // Pokemon Exception Handling
    fem.solve();
    return 0;
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
