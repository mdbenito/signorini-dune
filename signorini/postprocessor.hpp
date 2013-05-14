/******************************************************************************
 * postprocessor.hpp                                                          *
 *                                                                            *
 * Common postprocessing object. Computes Von Mises stress, error and outputs *
 * to VTK files.                                                              *
 ******************************************************************************/

#ifndef SIGNORINI_POSTPROCESSOR_HPP
#define SIGNORINI_POSTPROCESSOR_HPP

#include "config.h"
#include <iostream>
#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "shapefunctions.hpp"
#include "utils.hpp"
#include "benchmark.hpp"

using namespace Dune;

/*! 
 
 Template type names:
 TGV: TGridView
 TET: TElasticityTensor
 TMP: TMapper
 TSS: TShapeFunctionSet
 
 TODO: Should only take one template parameter: the FEM class, and from it
 read types for GridView, Elasticity tensor, shape function set, etc.
 */
template<class TGV, class TET, class TMP, class TSS>
class PostProcessor
{
public:
  static const int dim = TGV::dimension;
  
  typedef typename TGV::ctype                       ctype;
  typedef FieldVector<ctype, dim>                 coord_t;
  typedef FieldMatrix<ctype, dim, dim>            block_t;
  typedef BlockVector<coord_t>                CoordVector;
  typedef std::vector<ctype>                   FlatVector;

  typedef LeafMultipleCodimMultipleGeomTypeMapper
          <typename TGV::Grid, MCMGVertexLayout> VertexMapper;
  
private:
  const TGV& gv;      //!< Grid view
  const TMP& mapper;  //!< Index mapper
  const TET& a;       //!< Elasticity tensor
  
  
  CoordVector u;      //!< Solution. Uses VertexMapper ordering.
  FlatVector vm;      //!< Von Mises stress of solution. Uses VertexMapper ordering.

public:
  PostProcessor (const TGV& _gv, const TMP& _m, const TET& _a);

  long double computeError (const CoordVector& v);
  void   computeVonMisesSquared ();
  
  std::string writeVTKFile (std::string base, int step) const;
  
protected:
  void check (const CoordVector& v) const;
  void setSolution (const CoordVector& v);
};

template<class TGV, class TET, class TMP, class TSS>
PostProcessor<TGV, TET, TMP, TSS>::PostProcessor (const TGV& _gv,
                                                  const TMP& _m,
                                                  const TET& _a)
  : gv (_gv), mapper(_m), a (_a)
{
  u.resize (gv.size (dim), 0.0);
  vm.resize (gv.size (dim), 0.0);
}

/*! Copies the solution vector.
 
 Only those entries of the input vector which correspond to nodes of the mesh
 are copied, using the mapper provided. In particular, any additional entries
 are disregarded (e.g. lagrange multipliers at the end of the CoordVector, etc.)
 */
template<class TGV, class TET, class TMP, class TSS>
void PostProcessor<TGV, TET, TMP, TSS>::setSolution (const CoordVector& v)
{
  VertexMapper defaultMapper (gv.grid ());
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int from = mapper.map (*it), to = defaultMapper.map (*it);
    u[to] = v[from];
  }
  bench().report ("Postprocessing", "New solution copied");
}

/*! Computes the error wrt. to the previous solution.
 
 Sets the current solution to the argument and then computes the error as
 
       new - old / |new|
 
 Creates a copy of the data. The old solution data is deleted after computation.
 
 In case no previous solution was set, returns 1, but still copies the data.
 */
template<class TGV, class TET, class TMP, class TSS>
long double PostProcessor<TGV, TET, TMP, TSS>::computeError (const CoordVector& v)
{
  check (v);
  long double r = 0.0;
  VertexMapper defaultMapper (gv.grid ());
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int from = mapper.map (*it);
    int   to = defaultMapper.map (*it);
    r += coord_t (v[from] - u[to]).one_norm();
  }
  setSolution (v); // u will hold the new solution values
  r /= u.one_norm();
  bench().report ("Postprocessing", string ("New solution diverged by: ") + r);

  return r;  // Should be 1 on first run.
}

/*! Returns the *squared* von Mises stress. (brute force)
 
 FIXME: the discontinuity of the gradients of the basis functions at the nodes
 leads to very weird results. How can I fix it?

 Also: shouldn't I normalize the contribution of each vertex?
 */
template<class TGV, class TET, class TMP, class TSS>
void PostProcessor<TGV, TET, TMP, TSS>::computeVonMisesSquared ()
{
  bench().report ("Postprocessing", "Computing von Mises stress...", false);
  const auto& basis = TSS::instance ();
  VertexMapper defaultMapper (gv.grid());

  std::fill (vm.begin(), vm.end(), 0.0);
  std::vector<int> count (gv.size (dim), 1);
  assert (defaultMapper.size() == gv.size (dim));
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    const int   vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i) {
      auto ii = defaultMapper.map (*it, i, dim);
      ++count.at(ii);
      auto iipos = ref.position (i, dim);
      block_t  s = a.stress (u[ii], basis[i].evaluateGradient (iipos));
      double   t = trace(s);
      double r = -0.5*t*t + 1.5*trace (s.rightmultiplyany (s));
//      cout << "stressing: " << ii << " at " << iipos << LF;
      /*
      r = 0.0;
      for (int k=0; k<dim; ++k) {
        for (int l=0; l<dim; ++l) {
          r += (s[k][l] + (k==l ? -0.5*t : 0))*(s[k][l] + (k==l ? -0.5*t : 0));
        }
      }
       */
      vm.at(ii) += r;
    }
  }
  for (int i = 0; i < vm.size(); ++i)
    vm.at(i) /= count.at(i);

  bench().report ("Postprocessing", " done.");
}


template<class TGV, class TET, class TMP, class TSS>
std::string PostProcessor<TGV, TET, TMP, TSS>::writeVTKFile (std::string base, int step) const
{
  const int numVertices = gv.size (dim);
  std::ostringstream oss;
  oss << base << dim << "d-" << std::setfill('0') << std::setw(3) << step;
  VTKWriter<typename TGV::Grid::LeafGridView> vtkwriter (gv.grid().leafView());

  VertexMapper defaultMapper (gv.grid());
  
  // debugging indices
  /*
  std::vector<int> indices (numVertices);
  std::vector<int> mapped (numVertices);

  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int from = mapper.map (*it);
    int to = defaultMapper.map (*it);
    mapped.at(to) = from;
    indices.at(to) = to;
  }
  cout << "Adding index data" << LF;
  vtkwriter.addVertexData (indices, "idx", 1);
  vtkwriter.addVertexData (mapped, "map", 1);
   */
  FlatVector uu (numVertices * CoordVector::block_type::dimension);
  FlatVector vvmm (numVertices);
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int ii = defaultMapper.map (*it);
    for (int c = 0; c < CoordVector::block_type::dimension; ++c) {
      uu.at (ii*dim+c) = u[ii][c];
      vvmm.at (ii) = vm.at(ii);
    }
  }
  cout << "Adding vertex data" << LF;
  vtkwriter.addVertexData (uu, "u", dim);
  vtkwriter.addVertexData (vvmm, "vm", 1);
  
  cout << "Writing file" << LF;
  vtkwriter.write (oss.str(), VTK::ascii);
  bench().report ("Postprocessing", string ("Output written to ").append (oss.str()));
  return oss.str();
}


/*! Old sanity check.
 Wrong setup of the stiffness matrix would lead to NaN results in some places.
 This check should no longer be necessary, but it doesn't hurt either.
 */
template<class TGV, class TET, class TMP, class TSS>
void PostProcessor<TGV, TET, TMP, TSS>::check (const CoordVector& v) const {
  bool ok = true;
  for (const auto& it : v)
    for (int i = 0; i < dim; ++i)
      ok = ok && (it[i] != NAN);
  
  if (!ok) {
    std::cout << "*** NaN! *** \n" << "\tWhat a calamity, I choose to quit.\n";
    exit (1);
  }
}

#endif