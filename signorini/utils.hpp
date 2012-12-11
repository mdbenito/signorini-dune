/******************************************************************************
 * utils.hpp                                                                  *
 ******************************************************************************/

#ifndef SIGNORINI_UTILS_HPP
#define SIGNORINI_UTILS_HPP

#include <cmath>
#include "shapefunctions.hpp"

using std::cout;

template<class K, int R, int C>
FieldMatrix<K, R, C> operator* (const FieldMatrix<K, R, C>& m, const K& d)
{
  FieldMatrix<K, R, C> ret;
  for (int i=0; i < R; ++i)
    for (int j=0; j < C; ++j)
      ret[i][j] = m[i][j] * d;
  
  return ret;
}

template<class K, int N>
FieldVector<K, N> operator* (const FieldVector<K, N>& v, const K& d)
{
  FieldVector<K, N> ret;
  for (int i=0; i < N; ++i)
    ret[i] = v[i] * d;
  
  return ret;
}

/*! Should fully print a sparse matrix in Matlab syntax. It doesn't.
 
 NOTE: Dune already provides a printSparseMatrix()
 */
template<class M>
void printSparseMatrix (const M& m)
{
  cout << "M = [ ";
  for (const auto& r : m) {
    for (const auto& c : r)
      cout << c << "  ";
    cout << ";\n      ";
  }
  cout << " ]\n";
}


template <class VectorType>
void writeVectorToFile (const VectorType& vector,
                        const std::string& filename, int outputPrecision = 18)
{
  std::ofstream outs (filename.c_str ());
  outs.precision (outputPrecision);
  
  for (const auto& c : vector)
    outs << c << " ";
  outs << "\n";
  
  outs.close();
}

template <class T>
long double computeError (const std::vector<T>& uu, const std::vector<T>& vv)
{
  long double r = 0.0;
  long double n = 0.0;
  auto u = uu.begin();
  auto v = vv.begin();
  
  for (; u != uu.end() && v != vv.end(); ++u, ++v) {
    r += std::abs (*u - *v);
    n += std::abs (*u);
  }
  return r / n;
}

template <class ctype, int dim, template <class, int> class ShapeSet>
bool testShapes()
{
  const auto& basis = ShapeSet<ctype, dim>::instance ();
  GeometryType gt (basis.basicType, dim);
  const auto& element = GenericReferenceElements<ctype, dim>::general (gt);
  FieldVector<ctype, dim> x;
  
  cout << "Testing " << basis.basicType << " shapes for dimension " << dim << ":\n";
  cout << "Testing vertices:\n";
  for (int i=0; i < basis.N; ++i) {
    for (int v = 0; v < element.size (dim); ++v) {
        x = element.position (v, dim);
        cout << "   Basis[" << i << "](" << x << ") = "
             << basis[i].evaluateFunction (x) << "\n";
    }
  }
  
  cout << "Testing quadrature of order two:\n";
  for (int i=0; i < basis.N; ++i) {
    for (auto& x : QuadratureRules<ctype, dim>::rule (gt, 2)) {
        cout << "   Basis[" << i << "](" << x.position() << ") = "
             << basis[i].evaluateFunction (x.position()) << "\n";
    }
  }

  return true;
}


#endif  // SIGNORINI_UTILS_HPP
