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
#include <cstdlib>
#include <ctime>

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

template<class K, int C>
bool operator== (const FieldVector<K, C>& v, const K& m)
{
  for (int i=0; i < C; ++i)
//    if (std::abs(v[i] - m) <= std::numeric_limits<ctype>::epsilon())
    if (v[i] != m)
      return false;

  return true;
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
  virtual void pre (X& x, Y& b) { (void) x; (void) b;}
  virtual void apply (X& v, const Y& d) { v = d; }
  virtual void post (X& x) { (void) x; }
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

/* Gram-Schmidt (sort of)
 Returns a vector of v_1, ..., v_{dim-1} of orthonormal vectors determining the
 plane whose normal is this functions argument.
 
 FIXME: avoid the idiotic copy and removal of the argument.
 */
template <class ctype, int dim>
std::vector<FieldVector<ctype, dim> > basisOfPlaneNormalTo (const FieldVector<ctype, dim>& vec) {
  typedef FieldVector<ctype, dim> coord_t;

  std::vector<coord_t> v (dim);
  std::srand (static_cast<unsigned int> (std::time (0)));
  v[0] = vec;
  for (int i = 1; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      v[i][j] = (double)(std::rand() % 1000) * 1e-3;
    }
  }
  
  for (int i = 0; i < dim; ++i) {
    v[i] /= v[i].two_norm();
    for (int j = i+1; j < dim; ++j) {
      v[j] -= v[i] * (v[j] * v[i]); // v[i] * <v[j], v[i]>
      if (v[j] == 0.0)
        return basisOfPlaneNormalTo (vec);
    }
  }
  
  v.erase (v.begin ());
  return v;
}

  // FINISH THIS! (check orthonormality)
template <class ctype, int dim>
bool test_basisOfPlaneNormalTo () {
  typedef FieldVector<ctype, dim> coord_t;
  coord_t n(0); n[0] = n[1] = 1.0;
  auto r = basisOfPlaneNormalTo (n);
  cout << "n= " << n << LF;

  for (int i=0; i<dim-1; ++i)
    cout << "\t\tr[" << i << "]= " << r[i] << ", norm= " << r[i].two_norm() << LF;
  
  return true;
}

  // More hacks...
template <class K>
FieldVector<K, 2> coord2 (K x, K y) {
  FieldVector<K, 2> ret;
  ret[0] = x; ret[1] = y;
  return ret;
}

template <class K>
FieldVector<K, 3> coord3 (K x, K y, K z) {
  FieldVector<K, 3> ret;
  ret[0] = x; ret[1] = y; ret[2] = z;
  return ret;
}

/*** Poor man's polynomials ***/

template <class ctype, int dim>
class Monomial {
  typedef FieldVector<ctype, dim> coord_t;
  
  int indices[dim];  // array {0 1 1} means x_2*x_3, etc.
  ctype coeff;
  
  void clear() {
    for (int i=0; i < dim; ++i) indices[i] = 0;
  }
public:
  typedef std::vector<Monomial<ctype, dim> > Polynomial;
  
  Monomial (ctype _coeff=0) : coeff(_coeff) { clear(); }
  Monomial (std::vector<int> ilist, ctype _coeff) : coeff (_coeff) {
    clear();
    for (auto& x : ilist) indices[x] = 1;
  };
  
  ctype operator() (const coord_t& local) const
  {
    ctype r = coeff;
    for (int i=0; i < dim; ++i)
      if (indices[i] > 0)
        r *= local[i];
    
    return r;
  }
  
  void differentiate (int which)
  {
    assert (which >=0 && which < dim);
    if (indices[which] == 0) {
      clear();
      coeff = 0;
    } else {
      indices[which] = 0;
    }
  }
  
  /* TODO: translate this into C++!! and delete the rest!!!!
   
   (define dim 3)
   
   (define (prepend x)
   (lambda (y) (append x y)))
   
   (define (monomials ord from prog)
   (cond ((>= from dim) '())
   ((== ord 1) (cons `(,from) (monomials ord (+ 1 from) prog)))
   (else
   (with np (append prog `(,from))
   (append
   (map (prepend np) (monomials (- ord 1) (+ 1 from) '()))
   (monomials ord (+ 1 from) prog))))))
   */
  static Polynomial monomialsOfOrder (int order)
  {
    if (dim > 3)
      DUNE_THROW (Exception, "general monomialsOfOrder not implemented in c++");
    
    /* MEGA-HACK!!
     
     This is just plain stupid, but it's too late now for recursion in c++
     without sensible reference counting.
     */
    
    Polynomial p;
    std::vector<int> ids;
    
    if (order == 1) {
      for (int i=0; i < dim; ++i)
        p.push_back (Monomial (ids << clear_vector() << i, 1));
    } else if (order == 2) {
      if (dim == 2) {
        p.push_back (Monomial (ids << clear_vector() << 0 << 1, 1));
      } else if (dim == 3) {
        p.push_back (Monomial (ids << clear_vector() << 0 << 1, 1));
        p.push_back (Monomial (ids << clear_vector() << 0 << 2, 1));
        p.push_back (Monomial (ids << clear_vector() << 1 << 2, 1));
      }
    } else if (order == 3) {
      p.push_back (Monomial (ids << clear_vector() << 0 << 1 << 2, 1));
    }
    return p;
  }
  
  
  friend std::ostream& operator<< (std::ostream& os, const Monomial& m)
  {
    os << m.coeff;
    for (int i=0; i < dim; ++i)
      if (m.indices[i] > 0)
        os << "x[" << i << "]";
    return os;
  }
  
  friend Monomial& operator* (Monomial& m, ctype c)
  {
    m.coeff *= c;
    return m;
  }
  
  friend std::ostream& operator<< (std::ostream& os, const Polynomial& p)
  {
    for (auto& m : p)
      os << m << " + ";
    return os;
  }
  
  friend Polynomial& operator* (Polynomial& p, ctype c)
  {
    for (auto& x : p)
      x = x*c;
    return p;
  }
};

#endif  // SIGNORINI_UTILS_HPP
