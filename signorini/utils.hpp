/******************************************************************************
 * utils.hpp                                                                  *
 ******************************************************************************/

#ifndef SIGNORINI_UTILS_HPP
#define SIGNORINI_UTILS_HPP

#include "config.h"
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <set>
#include <vector>

using namespace Dune;
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

template<class K, int R, int C>
K trace (const FieldMatrix<K, R, C>& m)
{
  K ret(0);
  for (int i=0; i < R; ++i)
    ret += m[i][i];

  return ret;
}

template<class K, int R, int C>
FieldVector<K, C> operator* (const FieldVector<K, C>& v, const FieldMatrix<K, R, C>& m)
{
  FieldVector<K, C> ret (0.0);
  for (int i=0; i < C; ++i)
    for (int j=0; j < R; ++j)
      ret[i] += v[j] * m[j][i];
  
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

  // convenience funcs
template <class T>
std::vector<T>& operator<< (std::vector<T>& v, const T& d)
{
  v.push_back (d);
  return v;
}

template <class T>
std::set<T>& operator<< (std::set<T>& s, const T& d)
{
  s.insert (d);
  return s;
}

struct clear_vector { };

template <class T>
std::vector<T>& operator<< (std::vector<T>& v, const clear_vector& c)
{
  (void) c;
  v.clear();
  return v;
}
  ////// Shitty, wasteful, slow string manipulation:


std::string operator+ (const std::string& a, const std::string& b);

  // Use with integers, etc.
template <class T>
std::string operator+ (const std::string& a, T b)
{
  std::ostringstream oss;
  oss << a << b;
  return oss.str();
}

/*! Wrap a BCRSMatrix class and offer a window into it.
 
 Simply wraps operator[] and adds offsets.
 */
template <class T>
class MatrixWindow : public T
{
  typedef typename T::row_type row_type;
  typedef typename T::size_type size_type;

  class RowWindow : public row_type
  {
  public:
    int offset;
    row_type* r;

    RowWindow () : offset(0) { }
    T& operator[] (size_type i) { return (*r)[i+offset]; }
    const T& operator[] (size_type i) const { return (*r)[i+offset]; }
  };

  const T& M;
  const int roffset;
  const int coffset;
  
  mutable RowWindow row;

public:
  MatrixWindow (const T& m, int r, int c) : M(m), roffset(r), coffset(c)
  {
    row.offset = coffset;
  }
  
  RowWindow& operator[] (size_type i)
  {
    row.r = &M[i + roffset];
    return *row;
  }
  
  const RowWindow& operator[] (size_type i) const
  {
    row.r = &M[i + roffset];
    return *row;
  }
};

template <class ctype, int R, int C>
class MatrixFlattener : BCRSMatrix<ctype>
{
  typedef BCRSMatrix<FieldMatrix<ctype, R, C> > BaseMatrix;
  const BaseMatrix& M;

public:
  MatrixFlattener (const BaseMatrix& from) : M(from)
  {
    setSize (M.rows() * R, M.columns() * C);
      // 1. Set sparsity pattern from M
    
      // 2. Fill from M
  }
};


#endif  // SIGNORINI_UTILS_HPP
