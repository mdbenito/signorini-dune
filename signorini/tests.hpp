/******************************************************************************
 * tests.hpp                                                                  *
 ******************************************************************************/

#ifndef SIGNORINI_TESTS_HPP
#define SIGNORINI_TESTS_HPP

#include "config.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "shapefunctions.hpp"
#include "utils.hpp"

using namespace Dune;

template <class ctype, int dim>
bool test_basisOfPlaneNormalTo () {
  typedef FieldVector<ctype, dim> coord_t;
  static const double tolerance = 1.0e-15;
  
  std::cout << "Testing basis of plane normal to a vector n (tolerance= " << tolerance << "):\n";
  bool ret = true;
  for (int c = 0; c < 5; ++c) {
    coord_t n;
    for (int d = 0; d < dim; ++d)
      n[d] = 2.0 * (-0.5 + rand() / (RAND_MAX + 1.0));  // n[d] ∈ (-1,1)
    
    std::cout << "\tTest #" << c <<": n= " << n << LF;
    auto t = basisOfPlaneNormalTo (n);
    for (int i = 0; i < dim - 1; ++i) {
      ctype prod = n * t[i];
      bool orth = std::abs(prod) <= tolerance;
      ret = ret && orth;
      std::cout << "\t\tt[" << i << "]= " << t[i] << ", norm= " << t[i].two_norm();
      if (orth) std::cout << " is orthogonal to n.\n";
      else      std::cout << " is NOT orthogonal to n (" << prod << ")." << LF;
    }
    for (int i = 0; i < dim - 1; ++i) {
      for (int j = i+1; j < dim - 1; ++j) {
        ctype prod = t[i] * t[j];
        bool orth = std::abs(prod) <= tolerance;
        ret = ret && orth;
        std::cout << "\t\tt[" << i << "] and t[" << j << "]";
        if (orth) std::cout << " are orthogonal to each other.\n";
        else      std::cout << " are NOT orthogonal to each other (" << prod << ")." << LF;
      }
    }
  }
  return ret;
}


template <class ctype, int dim, class ShapeSet>
bool testShapes (const std::string& filename)
{
  typedef FieldVector<ctype, dim> coord_t;

  const auto& basis = ShapeSet::instance ();
  GeometryType gt (basis.basicType, dim);
  const auto& element = GenericReferenceElements<ctype, dim>::general (gt);
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
  
  cout << "Testing quadrature points of order two (TODO):\n";
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
  typedef typename TGV::ctype     ctype;
  typedef typename TGV::Grid     grid_t;
  typedef LeafMultipleCodimMultipleGeomTypeMapper <grid_t, MCMGVertexLayout> VertexMapper;
 
  const int     dim = TGV::dimension;
  const int ret_dim = TFN::return_dim;
  const int numVertices = gv.size(dim);

  VertexMapper mapper (gv.grid());
  VTKWriter<typename grid_t::LeafGridView> vtkwriter (gv.grid().leafView());
 
  std::vector<ctype> values (numVertices * ret_dim, 0.0);
  std::vector<ctype> support (numVertices, 0.0);

  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (func.isSupported (*is)) {
        const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
        for (int i = 0; i < ivnum; ++i) {
          int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          auto  idx = mapper.map (*it, subi, dim);
          support.at(idx) = 1.0;
          
          auto global = it->geometry().global (ref.position (subi, dim));
          auto  ret = func (*is, global);
          for (int c = 0; c < ret_dim; ++c)
            values.at (idx*ret_dim+c) = ret[c];
        }
      }
    }
  }

  vtkwriter.addVertexData (support, "support", 1);
  vtkwriter.addVertexData (values, "value", ret_dim);
  base = base + std::string("-") + dim + std::string("d");
  cout << "Writing file " << base  << LF;
  vtkwriter.write (base, VTK_OUTPUT_MODE);
}


template <class Glue, int body>
void testContactSurfaces (const std::vector<const Glue*>& gluedPieces, std::string base)
{
  dune_static_assert((body==0 || body==1), "'body' can only be 0 or 1");
  
  typedef typename SelectType <(body==0), typename Glue::Grid0View, typename Glue::Grid1View>::Type GV;
  typedef typename GV::ctype ctype;
  typedef LeafMultipleCodimMultipleGeomTypeMapper <typename GV::Grid, MCMGVertexLayout> VertexMapper;
  
  const int         dim = GV::dimension;
  const int     domdimw = GV::dimensionworld;
  const auto&        gv = gluedPieces.at(0)->template gridView<body>();
  const int numVertices = gv.size(dim);
  
  VertexMapper mapper (gv.grid());
  VTKWriter<GV> vtkwriter (gv);
  std::vector<ctype> support (numVertices, 0.0);
  std::vector<ctype> projection (numVertices * dim, 0.0);
  std::vector<ctype> normal (numVertices * dim, 0.0);
  for (auto glue : gluedPieces) {
    for (auto is = glue->template ibegin<body>(); is != glue->template iend<body>(); ++is) {
      if (is->self() && is->neighbor()) {
        const auto& ref = GenericReferenceElements<ctype, dim>::general (is->inside()->type());
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
        for (int i = 0; i < ivnum; ++i) {
          int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          auto  idx = mapper.map (*(is->inside()), subi, dim);
          support.at(idx) = 1.0;
          auto pos = is->geometryInInside().local (ref.position (i, dim));
          auto nr = is->geometryOutside().global (pos) - is->geometry().global (pos);
          nr /= nr.two_norm ();
          for (int c = 0; c < dim; ++c)
            projection.at (idx*dim + c) = nr[c];
          
          auto slaveVertex = is->geometry().corner (subi);
          auto localSlave = is->geometry().local (slaveVertex);
          nr = is->unitOuterNormal (localSlave);
          nr *= (body == 0) ? 1.0 : -1.0;
          for (int c = 0; c < dim; ++c)
            normal.at (idx*dim + c) = nr[c];
        }
      }
    }
  }
  vtkwriter.addVertexData (support, "contact", 1);
  vtkwriter.addVertexData (projection, "projection", dim);
  vtkwriter.addVertexData (normal, "normal", dim);

  base = base + std::string("-") + body + std::string("-") + dim + std::string("d");
  cout << "Writing file " << base  << LF;
  vtkwriter.write (base, VTK_OUTPUT_MODE);
}

#endif  // SIGNORINI_TESTS_HPP
