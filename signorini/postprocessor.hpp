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
 
 FIXME: Should only take one template parameter: the FEM class, and from it
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
  typedef BlockVector<FieldVector<ctype,1> > ScalarVector;
  typedef std::vector<ctype>                   FlatVector;

  typedef LeafMultipleCodimMultipleGeomTypeMapper
          <typename TGV::Grid, MCMGVertexLayout> VertexMapper;
  
private:
  const TGV& gv;      //!< Grid view
  const TMP& mapper;  //!< Index mapper
  const TET& a;       //!< Elasticity tensor
  
  
  CoordVector* u;      //!< Solution
  ScalarVector vm;    //!< Von Mises stress of solution

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
  : gv (_gv), mapper(_m), a (_a), u (NULL), vm (NULL)
{
  VertexMapper defaultMapper (gv.grid());
  const auto totalVertices = defaultMapper.size ();
  vm.resize (totalVertices);
}

template<class TGV, class TET, class TMP, class TSS>
void PostProcessor<TGV, TET, TMP, TSS>::setSolution (const CoordVector& v)
{
  delete u;
  u = new CoordVector (v);
  bench().report ("Postprocessing", string("New solution has size: ") + u->size());
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
  
  if (u == NULL) {
    r = 1.0;
  } else {
    long double n = 0.0;
    auto uit = u->begin();
    auto vit = v.begin();
    
    for (; uit != u->end() && vit != v.end(); ++uit, ++vit) {
      for (int i=0; i < dim; ++i) {
        r += std::abs ((*uit)[i] - (*vit)[i]);
        n += std::abs ((*uit)[i]);
      }
    }
    r = r / n;
    bench().report ("Postprocessing", string ("New solution diverged by: ") + r);
  }
  
  setSolution (v);
  return r;
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
  
  const auto totalVertices = defaultMapper.size ();
  
  if (u == NULL || u->size() < totalVertices)
    DUNE_THROW (Exception, "call PostProcessor::computeError() first");

  vm = 0;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType typ = it->type ();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (typ);
    const int   vnum = ref.size (dim);

    double r;
    for (int i = 0; i < vnum; ++i) {
      auto ii = defaultMapper.map (*it, i, dim);
      auto iipos = ref.position (i, dim);
      block_t  s = a.stress ((*u)[ii], basis[i].evaluateGradient (iipos));
      double   t = trace(s);
      r = -0.5*t*t + 1.5*trace (s.rightmultiplyany (s));
        //cout << "stressing: " << ii << " ";
      /*
      r = 0.0;
      for (int k=0; k<dim; ++k) {
        for (int l=0; l<dim; ++l) {
          r += (s[k][l] + (k==l ? -0.5*t : 0))*(s[k][l] + (k==l ? -0.5*t : 0));
        }
      }
       */
      vm[ii] += r / vnum;
    }
  }
  
  bench().report ("Postprocessing", " done.");
}


  /// This is about the ugliest code I've written in a while...
  /// UPD: scrap that. I just wrote something far worse. :(((((

template<class TGV, class TET, class TMP, class TSS>
std::string PostProcessor<TGV, TET, TMP, TSS>::writeVTKFile (std::string base, int step) const
{
  std::ostringstream oss;
  oss << base << dim << "d-" << std::setfill('0') << std::setw(3) << step;
  VTKWriter<typename TGV::Grid::LeafGridView> vtkwriter (gv.grid().leafView());
  
  std::vector<int> indices (gv.size(dim));
  std::vector<int> mapped (gv.size(dim));
  VertexMapper defaultMapper (gv.grid());
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int from = mapper.map (*it);
    int to = defaultMapper.map (*it);
    mapped[to] = from;
    indices[to] = to;
  }
  cout << "Adding vertex data" << LF;
  vtkwriter.addVertexData (indices, "idx", 1);
  vtkwriter.addVertexData (mapped, "map", 1);

  FlatVector uu (gv.size(dim) * CoordVector::block_type::dimension);
  FlatVector vvmm (gv.size(dim));
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    int from = mapper.map (*it);
    int to = defaultMapper.map (*it);
    for (int c = 0; c < CoordVector::block_type::dimension; ++c) {
      uu[to*dim+c] = (*u)[from][c];
      vvmm[to] = vm[from];
    }
  }
  cout << "Adding more vertex data" << LF;
  vtkwriter.addVertexData (uu, "u", dim);
  vtkwriter.addVertexData (vvmm, "vm", 1);
  
  cout << "Writing vertex data" << LF;
  vtkwriter.write (oss.str(), VTK::appendedraw);
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