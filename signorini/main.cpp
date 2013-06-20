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

typedef enum { PLATE=0, PRISM=1, CYLINDER=2, PEG=3, INDENT=4, COGS=5, TOOTH=6 } ProblemType;
static const std::string meshPath ("/Users/miguel/Devel/Signorini/meshes/");
static const std::string meshNames[2][7] = {
  {"plate master.msh", "prism master.msh", "cylinder master.msh", "prism master.msh", "plate master.msh", "cog master.msh", "tooth master.msh"},
  {"plate slave.msh",  "prism slave.msh",  "cylinder slave.msh", "cylinder vertical.msh", "indent slave.msh", "cog slave.msh", "tooth slave.msh"}
};

int main (int argc, char** argv)
{
  std::srand (clock ());   // Needed somewhere...
  
  #define         DIM   2                 // HACK because of VectorEval...
  ProblemType problem = COGS;
  const bool    tests = false;
  const int       dim = DIM;
  const double   E[2] = { 1.0e9, 1.0e9 };
  const double  nu[2] = { 0.3,     0.3 };

  typedef UGGrid<dim>              grid_t;
  typedef GridFactory<grid_t>   factory_t;
  typedef grid_t::ctype             ctype;
  typedef FieldVector<ctype, dim> coord_t;
  typedef FieldVector<ctype, 1>  scalar_t;
  typedef grid_t::LeafGridView         GV;
  typedef Codim1Extractor<GV>                    SurfaceExtractor;
  typedef PSurfaceMerge<dim-1, dim, double>      SurfaceMergeImpl;
  typedef ::GridGlue<SurfaceExtractor, SurfaceExtractor> GlueType;
  typedef PhysicalFaceDescriptor<GV, factory_t>    FaceDescriptor;
  typedef HookeTensor<ctype, dim>                          HookeT;
  typedef ConstantEvaluation<ctype, dim, coord_t>      VectorEval;
  typedef ConstantEvaluation<ctype, dim, scalar_t>     ScalarEval;
  typedef std::map<int, shared_ptr<VectorEval> >   VectorEvalsMap;
  typedef std::map<int, shared_ptr<ScalarEval> >   ScalarEvalsMap;
  typedef GmshVolumeFunctor<ctype, dim, factory_t, VectorEval>      VolumeF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval>   NeumannF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, VectorEval> DirichletF;
  typedef GmshBoundaryFunctor<ctype, dim, factory_t, ScalarEval>   ContactF;
  
  /*
    //////HACK HACK HACK: what I want is a linked list of Constraints, but each has its
    // own type, so for now I'm including the type of the next element in the list
    // yep, this defeats its whole purpose...
  typedef Constraint<coord_t, dim, DirichletF, DummyConstraint> EndDirichlet;
  typedef Constraint<coord_t, dim, ContactF, EndDirichlet>           Contact;
  typedef Constraint<coord_t, dim, DirichletF, Contact>            Dirichlet;
  typedef Constraint<coord_t, dim, NeumannF, Dirichlet>              Neumann;
   */
  
  typedef P1ShapeFunctionSet<ctype, dim, 1, 0> ShapeSet;
  typedef LagrangeShapeFunctionSet<ctype, dim> LSShapeSet;
  typedef TwoBodiesIASet<GV, HookeT, VolumeF, DirichletF, NeumannF, GlueType, ShapeSet, LSShapeSet> TwoSolver;
  
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
   //  typedef SignoriniFEPenalty<GV, HookeT, VolumeF, BoundaryF, Gap, Dirichlet, ShapeSet> PMSolver;
   //  typedef SignoriniIASet<GV, HookeT, VolumeF, Dirichlet, BoundaryF, Gap, ShapeSet, LSShapeSet> IASolver;
   */

  factory_t        factories[2];
  grid_t*              grids[2];
  GmshReader<grid_t> readers[2];
  std::vector<int>     bi2pe[2]; // Boundary ids to gmsh physical entities
  std::vector<int>     ei2pe[2]; // Element indices to gmsh physical entities

  std::set<int> volumeGroups[2], contactGroups[2], dirichletGroups[2], neumannGroups[2];
  HookeT* a[2]; VolumeF* f[2];  DirichletF* d[2];  NeumannF* p[2];  ContactF* g[2];
  VectorEvalsMap fEvals[2], dEvals[2], pEvals[2];
  ScalarEvalsMap cEvals[2];
  FaceDescriptor*  descriptors[2];
  SurfaceExtractor* extractors[2];

  switch (problem) {
    case PLATE:
    case INDENT:
      contactGroups  [MASTER] << 3;       contactGroups  [SLAVE] << 1;
#if DIM == 2
      fEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 1e8)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, -1e8)));
      dEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      dEvals[SLAVE][3]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, -0.01)));
      pEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      pEvals[MASTER][4] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      pEvals[SLAVE][2]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      pEvals[SLAVE][4] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      cEvals[MASTER][3] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][1]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
      break;
    case COGS:
      contactGroups  [MASTER] << 1;       contactGroups  [SLAVE] << 1;
#if DIM == 2
      fEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 2e8)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, -2e8)));
      dEvals[MASTER][3] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      dEvals[SLAVE][3]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      pEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord2 (-3e6, 0.0)));
      pEvals[SLAVE][2] = shared_ptr<VectorEval> (new VectorEval (coord2 (3e6, 0.0)));
      cEvals[MASTER][1] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][1]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
      break;
    case TOOTH:
      contactGroups  [MASTER] << 1;       contactGroups  [SLAVE] << 1;
#if DIM == 2
      fEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      dEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      dEvals[SLAVE][2]  = shared_ptr<VectorEval> (new VectorEval (coord2 (0.0, 0.0)));
      pEvals[MASTER][3] = shared_ptr<VectorEval> (new VectorEval (coord2 (-3e6, 0.0)));
      pEvals[SLAVE][3] = shared_ptr<VectorEval> (new VectorEval (coord2 (3e6, 0.0)));
      cEvals[MASTER][1] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][1]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
      break;
    case PRISM:
      contactGroups  [MASTER] << 1;      contactGroups  [SLAVE] << 1;
#if DIM == 3
      fEvals[MASTER][1] = shared_ptr<VectorEval>(new VectorEval (coord3 (0.0,  0.0, 0.0)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      dEvals[MASTER][6] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.01, 0.0)));
      dEvals[SLAVE][6]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, -0.01, 0.0)));
      pEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[MASTER][3] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[MASTER][4] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[MASTER][5] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][2]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][3]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][4]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][5]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      cEvals[MASTER][1] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][1]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
      break;
    case CYLINDER:
      contactGroups  [MASTER] << 3 << 4;  contactGroups  [SLAVE] << 4;
#if DIM == 3
      fEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0,  4e7, 0.0)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, -4e7, 0.0)));
      dEvals[MASTER][1] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      dEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      dEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      dEvals[SLAVE][2]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[MASTER][5] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[MASTER][6] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][3]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][5]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      cEvals[MASTER][3] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[MASTER][4] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][4]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
      break;
    case PEG:
      contactGroups  [MASTER] << 1;  contactGroups  [SLAVE] << 1;
#if DIM == 3
      fEvals[MASTER][1] = shared_ptr<VectorEval>(new VectorEval (coord3 (0.0,  4e8, 0.0)));
      fEvals[SLAVE][1]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, -4e7, 0.0)));
      dEvals[MASTER][6] = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.01, 0.0)));
      dEvals[SLAVE][2]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, -0.01, 0.0)));
      pEvals[MASTER][2] = shared_ptr<VectorEval> (new VectorEval (coord3 (3e6, 0.0, 0.0)));
      pEvals[MASTER][3] = shared_ptr<VectorEval> (new VectorEval (coord3 (3e6, 0.0, 0.0)));
      pEvals[MASTER][4] = shared_ptr<VectorEval> (new VectorEval (coord3 (-3e6, 0.0, 0.0)));
      pEvals[MASTER][5] = shared_ptr<VectorEval> (new VectorEval (coord3 (-3e6, 0.0, 0.0)));
      pEvals[SLAVE][3]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][4]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][5]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      pEvals[SLAVE][6]  = shared_ptr<VectorEval> (new VectorEval (coord3 (0.0, 0.0, 0.0)));
      cEvals[MASTER][1] = shared_ptr<ScalarEval> (new ScalarEval (1.0));
      cEvals[SLAVE][1]  = shared_ptr<ScalarEval> (new ScalarEval (1.0));
#endif
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
      // the glue has then both views pointing to the master... WTF?!?!?
//    extractors[body]  = new SurfaceExtractor (grids[body]->leafView(), *descriptors[body]);
    
    a[body] = new HookeT (E[body], nu[body]);
    f[body] = new VolumeF (factories[body], ei2pe[body], &fEvals[body]);
    d[body] = new DirichletF (factories[body], bi2pe[body], &dEvals[body]);
    p[body] = new NeumannF (factories[body], bi2pe[body], &pEvals[body]);
    g[body] = new ContactF (factories[body], bi2pe[body], &cEvals[body]);

    /*
      //// Build constraints, taking the interfaces into account!
      // Nodes marked as Dirichlet have "priority"
    dd[body] = new Dirichlet (*d[body], true);
      // Neumann nodes have the lowest priority (are overrriden by any other type)
    pp[body] = new Neumann (*p[body], true, new Dirichlet (*d[body], false, new Contact (*g[body], false)));
      // Contact nodes are overriden only by Dirichlet nodes
    gg[body] = new Contact (*g[body], true, new EndDirichlet (*d[body], false));
     */
  }
    // HACK: see above why we don't do this in the loop.
  extractors[MASTER] = new SurfaceExtractor (grids[MASTER]->leafView(), *descriptors[MASTER]);
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
      testGmshBoundaryFunctor (gv[body], *g[body], string ("/tmp/test-gap-") + body);
    }
    testContactSurfaces<GlueType, MASTER> (glue, "/tmp/test-glue");
    testContactSurfaces<GlueType, SLAVE> (glue, "/tmp/test-glue");
    exit (1);
  }

//  PMSolver fem (gv[MASTER], a[MASTER], f[MASTER], p[MASTER], g[MASTER], d[MASTER], 1.0e-14 / E, 1e-6, 10);
//  IASolver fem (gv[MASTER], a[MASTER], f[MASTER], d[MASTER], p[MASTER], g[MASTER]);
  TwoSolver fem (gv[MASTER], gv[SLAVE], *a[MASTER], *a[SLAVE], *f[MASTER], *f[SLAVE], *d[MASTER], *d[SLAVE], *p[MASTER], *p[SLAVE], glue);
  
  try {   // Pokemon Exception Handling
    fem.solve();
    return 0;
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
  }
  return 1;
}
