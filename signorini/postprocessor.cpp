
/******************************************************************************
 * postprocessor.cpp                                                          *
 *                                                                            *
 * Common postprocessing object. Computes Von Mises stress, error and outputs *
 * to VTK files.                                                              *
 ******************************************************************************/
#include "config.h"
#include <iostream>
#include <vector>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "shapefunctions.hpp"
#include "utils.hpp"

using namespace Dune;
using std::cout;

/*! 
 
 Template type names:
 TGV: TGridView
 TET: TElasticityTensor
 TSS: TShapeFunctionSet
 
 FIXME: Should only take one template parameter: the FEM class, and from it
 read types for GridView, Elasticity tensor, shape function set, etc.
 */
template<class TGV, class TET, class TSS>
class PostProcessor
{
public:
  static const int dim = TGV::dimension;
  
  typedef            typename TGV::ctype ctype;
  typedef      FieldVector<ctype, dim> coord_t;
  typedef FieldMatrix<ctype, dim, dim> block_t;
  typedef     BlockVector<coord_t> CoordVector;
        // FIXME: not the right type!!!!
  typedef BlockVector<FieldVector<ctype,1> > ScalarVector;
  
private:
  const TGV& gv;  //!< Grid view
  const TET& a;   //!< Elasticity tensor
  
  CoordVector* u;      //!< Solution
  ScalarVector* vm;    //!< Von Mises stress of solution

public:
  PostProcessor (const TGV& _gv,  const TET& _a);

  void setSolution (const CoordVector& v);
  double computeError (const CoordVector& v);
  void computeVonMises ();
  
  void writeVTKFile (std::string base, int step) const;
  void check (const CoordVector& v) const;
  
protected:
  template<class V>
  std::vector<ctype> asVector(const V* u) const;};

template<class TGV, class TET, class TSS>
PostProcessor<TGV, TET, TSS>::PostProcessor (const TGV& _gv,  const TET& _a)
  : gv (_gv), a (_a), u (NULL), vm (NULL)
{

}

template<class TGV, class TET, class TSS>
void PostProcessor<TGV, TET, TSS>::setSolution (const CoordVector& v)
{
  if (u != NULL)
    delete u;
  u = new CoordVector (v);
}

/*! Computes the error wrt. to the previous solution.
 
 Sets the current solution to the argument and then computes the error as
 
       new - old / |new|
 
 Creates a copy of the data. The old solution data is deleted after computation.
 
 In case no previous solution was set, returns 1, but still copies the data.
 */
template<class TGV, class TET, class TSS>
double PostProcessor<TGV, TET, TSS>::computeError (const CoordVector& v)
{
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
  }
  setSolution (v);
  return r;
}

/*! Brute force! 
 
 Don't think, just code...
 
 Also: shouldn't I normalize the contribution of each vertex?
 */
template<class TGV, class TET, class TSS>
void PostProcessor<TGV, TET, TSS>::computeVonMises ()
{
  const auto& basis = TSS::instance ();
  const auto&  iset = gv.indexSet ();

  if (vm != NULL) delete vm;
  
  vm = new ScalarVector (gv.size (dim), gv.size (dim));

  for (auto& p : *vm)
    p = FieldVector<ctype,1>(0.0);
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType typ = it->type ();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (typ);
    const int   vnum = ref.size (dim);
    
    for (int i = 0; i < vnum; ++i) {
      auto ii = iset.subIndex (*it, i, dim);
      auto iipos = ref.position (i, dim);

      block_t s = a.stress ((*u)[ii], basis[i].evaluateGradient (iipos));
      double  r = trace(s);      
      r = r*r - 1.5*(r*r - trace (s.rightmultiplyany (s)));
      r = std::sqrt (r);
      
      (*vm)[ii] += r;
    }
  }
}

template<class TGV, class TET, class TSS>
void PostProcessor<TGV, TET, TSS>::writeVTKFile (std::string base, int step) const
{
  std::ostringstream oss;
  oss << base << dim << "d-" << std::setfill('0') << std::setw(3) << step;
  VTKWriter<typename TGV::Grid::LeafGridView> vtkwriter (gv.grid().leafView());
  if (u != NULL)  vtkwriter.addVertexData (asVector (u), "u", dim);
  if (vm != NULL) vtkwriter.addVertexData (asVector (vm), "vm", 1);
  vtkwriter.write (oss.str(), VTKOptions::binaryappended);
}

// FIXME: lots of copying! Choose some memory policy and fix this.
template<class TGV, class TET, class TSS>
template <class V>
std::vector<typename TGV::ctype>
PostProcessor<TGV, TET, TSS>::asVector(const V* v) const {
  std::vector<ctype> ret;
  for (auto& p : *v)
    for (auto& c : p)
      ret.push_back (c);
  return ret;
}

template<class TGV, class TET, class TSS>
void PostProcessor<TGV, TET, TSS>::check (const CoordVector& v) const {
  bool ok = true;
  for (const auto& it : v)
    for (int i = 0; i < dim; ++i)
      ok = ok && (it[i] != NAN);
  
  if (!ok) {
    cout << "*** NaN! *** \n" << "\tWhat a calamity, I choose to quit.\n";
    exit (1);
  }
}