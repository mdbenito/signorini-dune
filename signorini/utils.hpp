/******************************************************************************
 * utils.hpp                                                                  *
 ******************************************************************************/

#ifndef SIGNORINI_UTILS_HPP
#define SIGNORINI_UTILS_HPP

#include "config.h"
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/common/exceptions.hh>

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <set>
#include <vector>

#define LF std::endl

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

template<class M, class X, class Y, int l=1>
class DummyPreconditioner : public Preconditioner<X,Y> {
public:
  typedef typename Dune::remove_const<M>::type matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;
  enum { category=SolverCategory::sequential };
  DummyPreconditioner () { }
  virtual void pre (X& x, Y& b) { }
  virtual void apply (X& v, const Y& d) { v = d; }
  virtual void post (X& x) { }
};

/* Now, THIS is ugly... */
template <class T>
class TwoRefs {
  const T& a;
  const T& b;
  
public:
  TwoRefs (const T& _a, const T& _b) : a(_a), b(_b) { }
    //TwoRefs (T& _a, T& _b) : a(_a), b(_b) { }
  TwoRefs (const TwoRefs<T>& other) : a(other.a), b(other.b) { }
  TwoRefs& operator= (const T& other) { a = other.a; b = other.b; return *this;}

  const T& operator[] (int idx) const throw (Exception) {
    if (idx == 0) return a;
    if (idx == 1) return b;
    throw (new Exception ());
  }
/*
  T& operator[] (int idx) throw (Exception) {
    if (idx == 0) return a;
    if (idx == 1) return b;
    throw (new Exception ());
  }
*/
};

  // Flatten a BCRSMatrix<B> to a BCRSMatrix<B::field_type>
template <class B, class BB>
void flattenMatrix (const BCRSMatrix<B>& A, BCRSMatrix<BB>& dest)
{
  int roff = 0;
  int coff = 0;

    // Flatten the adjacency pattern
  std::map<int, std::set<int> > adjacencyPattern;
  
  for (int r = 0; r < A.N(); ++r) {
//    cout << "At row: " << r << LF;
    int rtmp = 0;
    coff = 0;
    for (int c = 0; c < A.M(); ++c) {
//      cout << "   At column: " << c << LF;
      if (! A.exists(r, c)) {    // Empty entry. Traverse the column to find its width
        bool done=false;
        int l = 0;
        for (; l < A.N(); ++l) if (A.exists (l,c)) { done = true; break; }
        if (done) {
          coff += A[l][c].M();
//          cout << "   EMPTY entry. Colum offset is now: " << coff << LF;
        } else {
          cout << "WTF?! empty column in A" << LF;
          exit (1);
        }
      } else {
        rtmp = A[r][c].N();
        for (int rr = 0; rr < A[r][c].N(); ++rr) {
//          cout << "      At sub-row: " << rr << LF;
          for (int cc = 0; cc < A[r][c].M(); ++cc) {
//            cout << "         At sub-column: " << cc << LF;
            if (A[r][c].exists(rr, cc)) {
//              cout << "         Inserting (" << roff+rr << ", " << coff+cc << ")" << LF;
              adjacencyPattern[roff+rr].insert((int)coff+cc);
            }
          }
        }   // End sub-block in current row
        coff += A[r][c].M();
//        cout << "   Column offset is now: " << coff << LF;
      }
    } // End row
    roff += rtmp;
//    cout << "Row offset is now: " << roff << LF;
  }
  
//  cout << "Size: (" << roff << ", " << coff << ")\n";
  dest.setBuildMode (BCRSMatrix<BB>::row_wise);
  dest.setSize (roff, coff);
  
  for (auto row = dest.createbegin(); row != dest.createend(); ++row)
    for (const auto& col : adjacencyPattern[row.index()])
      row.insert (col);

    // Copy A

  roff = 0;
  coff = 0;
  for (int r = 0; r < A.N(); ++r) {
//    cout << "At row: " << r << LF;
    int rtmp = 0;
    coff = 0;
    for (int c = 0; c < A.M(); ++c) {
//      cout << "   At column: " << c << LF;
      if (! A.exists(r, c)) {    // Empty entry. Traverse the column to find its width
        bool done=false;
        int l = 0;
        for (; l < A.N(); ++l) if (A.exists (l,c)) { done = true; break; }
        if (done) {
          coff += A[l][c].M();
//          cout << "   EMPTY entry. Colum offset is now: " << coff << LF;
        } else {
          cout << "WTF?! empty column in A" << LF;
          exit (1);
        }
      } else {
        rtmp = A[r][c].N();
        for (int rr = 0; rr < A[r][c].N(); ++rr) {
//          cout << "      At sub-row: " << rr << LF;
          for (int cc = 0; cc < A[r][c].M(); ++cc) {
//            cout << "         At sub-column: " << cc << LF;
            if (A[r][c].exists(rr, cc)) {
//              cout << "         Copying (" << roff+rr << ", " << coff+cc << ")" << LF;
              dest[roff+rr][coff+cc] = A[r][c][rr][cc];
            }
          }
        }   // End sub-block in current row
        coff += A[r][c].M();
//        cout << "   Column offset is now: " << coff << LF;
      }
    } // End row
    roff += rtmp;
//    cout << "Row offset is now: " << roff << LF;
  }
}


template <class B>
void subMatrix (const BCRSMatrix<B>& A, BCRSMatrix<B>& dest, int r, int c, int n, int m)
{
  dest.setBuildMode (BCRSMatrix<B>::random);
  dest.setSize (n, m);
  
  for (int i = 0; i < n; ++i) {
    int nz = 0;
    for (int j = 0; j < m; ++j) {
      if (A.exists (r+i, c+j))
        ++nz;
    }
    dest.setrowsize (i, nz);
  }
  dest.endrowsizes();
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      if (A.exists (r+i, c+j))
        dest.addindex (i, j);
  dest.endindices();
  
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      if (A.exists (r+i, c+j))
        dest[i][j] = A[r+i][c+j];
}

template <class B>
void transpose (BCRSMatrix<B>& dest, const BCRSMatrix<B>& A) {
  dest.setBuildMode (BCRSMatrix<B>::random);
  dest.setSize (A.M(), A.N());

  for (int col = 0; col < A.M(); ++col) {
    int nz = 0;
    for (int row = 0; row < A.N(); ++row) {
      if (A.exists (row, col))
        ++nz;
    }
    dest.setrowsize (col, nz);
  }
  dest.endrowsizes();

  for (int col = 0; col < A.M(); ++col)
    for (int row = 0; row < A.N(); ++row)
      if (A.exists (row, col))
        dest.addindex (col, row);
  dest.endindices();
  
  for (int col = 0; col < A.M(); ++col)
    for (int row = 0; row < A.N(); ++row)
      if (A.exists (row, col))
        dest[col][row] = A[row][col];
}

#endif  // SIGNORINI_UTILS_HPP
