/******************************************************************************
 * tests.hpp                                                                  *
 ******************************************************************************/

#ifndef SIGNORINI_TESTS_HPP
#define SIGNORINI_TESTS_HPP

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

using namespace Dune;

template <class ctype, int dim, class ShapeSet>
bool testShapes (const std::string& filename)
{
  const auto& basis = ShapeSet::instance ();
  GeometryType gt (basis.basicType, dim);
  const auto& element = GenericReferenceElements<ctype, dim>::general (gt);
  typedef FieldVector<ctype, dim> coord_t;
  coord_t x;
  
  cout << "Testing " << basis.basicType << " shapes for dimension " << dim << ":\n";
  cout << "Testing vertices:\n";
  for (int i=0; i < basis.size(); ++i) {
    for (int v = 0; v < element.size (dim); ++v) {
      x = element.position (v, dim);
      cout << "   Basis[" << i << "](" << x << ") = "
           << basis[i].evaluateFunction (x) << "\n";
    }
  }
  
  cout << "Testing grid:\n";
  FieldMatrix<ctype, 41, 41> r, gx, gy;
  
  for (int i=0; i < basis.size(); ++i) {
    for (int x = -20; x <= 20; ++x) {
      for (int y = -20; y <= 20; ++y) {
        coord_t v;
        if (dim == 2)      v <<= x*0.1 , y*0.1;
        else if (dim == 3) v <<= x*0.1 , y*0.1, 0.5;
        else DUNE_THROW (Exception, "Invalid dimension for testShapes");
        
        r[x+20][y+20]  = basis[i].evaluateFunction (v);
        gx[x+20][y+20] = basis[i].evaluateGradient (v)[0];
        gy[x+20][y+20] = basis[i].evaluateGradient (v)[1];
          //        cout << "   Basis[" << i << "](" << x << ", " << y << ") = "
          //            << basis[i].evaluateFunction (v) << "\n";
      }
    }
    writeMatrixToMatlab (r, filename + i);
    writeMatrixToMatlab (gx, filename + std::string("-gradx") + i);
    writeMatrixToMatlab (gy, filename + std::string("-grady") + i);
  }
  
  cout << "Testing gradient:\n";
  for (int i=0; i < basis.size(); ++i) {
    for (int v = 0; v < element.size (dim); ++v) {
      x = element.position (v, dim);
      cout << "   Basis[" << i << "](" << x << ") = "
           << basis[i].evaluateGradient (x) << "\n";
    }
  }
  
  cout << "Testing quadrature points of order two:\n";
  for (int i=0; i < basis.size(); ++i) {
    for (auto& x : QuadratureRules<ctype, dim>::rule (gt, 2)) {
      cout << "   Basis[" << i << "](" << x.position() << ") = "
           << basis[i].evaluateFunction (x.position()) << "\n";
    }
  }

  return true;
}


/*
TGV: TGridView
TFN: TFunctor
*/
template<class TGV, class TFN>
void testGmshBoundaryFunctor (const TGV& gv, const TFN& func, std::string base)
{
  static const int     dim = TGV::dimension;
  static const int ret_dim = TFN::return_dim;
  
  typedef typename TGV::ctype     ctype;
  typedef typename TGV::Grid     grid_t;
  typedef LeafMultipleCodimMultipleGeomTypeMapper <grid_t, MCMGVertexLayout> VertexMapper;
  
  VertexMapper mapper (gv.grid());
  const auto& gids = gv.grid().globalIdSet();
  
  VTKWriter<typename grid_t::LeafGridView> vtkwriter (gv.grid().leafView());
 
  std::vector<ctype> values (gv.size(dim) * ret_dim, 0.0);
  std::vector<ctype> support (gv.size (dim), 0.0);
//  std::vector<int> indices (gv.size (dim), 0);

  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is)
      if (func.isSupported (*is)) {
        const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
        for (int i = 0; i < ivnum; ++i) {
          int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          auto  idx = mapper.map (*it, subi, dim);
          const auto& in = is->inside();
          support[idx] = 1.0;
//          indices[idx] = idx;

          auto global = it->geometry().global (ref.position (subi, dim));
          auto  ret = func (global);
          for (int c = 0; c < ret_dim; ++c)
            values[idx*ret_dim+c] = ret[c];
        }
      }

  vtkwriter.addVertexData (support, "support", 1);
  vtkwriter.addVertexData (values, "value", ret_dim);
//  vtkwriter.addVertexData (indices, "index", 1);
  base = base + std::string("-") + dim + std::string("d");
  cout << "Writing file " << base  << LF;
  vtkwriter.write (base, VTK::appendedraw);
}

#endif  // SIGNORINI_TESTS_HPP
