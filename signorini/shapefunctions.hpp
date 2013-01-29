/******************************************************************************
 * shapefunctions.hpp                                                         *
 *                                                                            *
 * TODO:                                                                      *
 *   - Replace *Q1ShapeFunctionSet for a template class, taking the shape as    *
 *     template parameter (needs changing the definition of P1 functions).    *
 ******************************************************************************/

#ifndef SIGNORINI_SHAPEFUNCTIONS_HPP
#define SIGNORINI_SHAPEFUNCTIONS_HPP

#include <iostream>
//#include <bitset>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include "utils.hpp"

using namespace Dune;
using std::cout;

/*! LinearShapeFunction:
 
 Represents a shape function and provides methods to evaluate the function
 and its gradient on a reference element.
 */
template<class ctype, int dim>
class LinearShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;
  
public:
  enum { N = dim };
  
  LinearShapeFunction () : coeff0(0.0), coeff1(0.0) {}
  
  LinearShapeFunction (ctype coeff0_, const coord_t& coeff1_)
  : coeff0(coeff0_), coeff1(coeff1_) {}
  
  ctype evaluateFunction (const coord_t& local) const
  {
    ctype result = coeff0;
    for (int i = 0; i < dim; ++i)
      result += coeff1[i] * local[i];
    return result;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const { return coeff1; }
  
private:
  ctype   coeff0;
  coord_t coeff1;
  
  template<class TC, int d> friend class P1ShapeFunctionSet;
};

/*! P1ShapeFunctionSet
 
 Singleton collection of LinearShapeFunction for P1 elements.
 */
template<class ctype, int dim>
class P1ShapeFunctionSet
{
  typedef FieldVector<ctype, dim> coord_t;
  
public:
  enum { N = dim + 1 };
  const GeometryType::BasicType basicType = GeometryType::simplex;

  typedef LinearShapeFunction<ctype, dim> ShapeFunction;
  
  static const P1ShapeFunctionSet& instance()
  {
    static const P1ShapeFunctionSet sfs;
    return sfs;
  }
  
  inline const ShapeFunction& operator[] (int i) const
  {
    if (i < 0 || i >= N)
      DUNE_THROW (Exception, "Index out of bounds for shape function.");
    return (i == 0) ? f0 : f1[i - 1];
  }
  
private:

  P1ShapeFunctionSet()
  {
    coord_t e(-1.0);
    f0.coeff0 = 1.0;
    f0.coeff1 = e;
    for (int i = 0; i < dim; ++i) {
      coord_t e(0.0);
      e[i] = 1.0;
      f1[i].coeff0 = 0.0;
      f1[i].coeff1 = e;
    }
  }
  
  P1ShapeFunctionSet(const P1ShapeFunctionSet& other) {}
  
  ShapeFunction f0;
  ShapeFunction f1[dim];
};


/*! MLinearShapeFunction:
 
 Represents a shape function and provides methods to evaluate the function
 and its gradient on a cubic reference element.
 
 Any Lagrange shape function is identified by the coefficients of the vertex
 were it's 1. On the plane, vertex (0,0) has function (1-x)(1-y), vertex (0,1)
 function (1-x)y, vertex (1,1) function xy, and so on.

 With this scheme (product of factors being x_i or (1-x_i)) the functions are
 linear on the restriction to each edge originating at the vertex where they are
 one, since they are polynomes of order N with N-1 components fixed.
 
 ...
 
 The mask data member is actually just the coordinates of the vertex where this
 shape function is 1, i.e. we use the bits:
 
    mask = 5 = 0...0101 => vertex = (...,0,1,0,1)

 TODO:
  - Test unit!!
  - Memoization seems unnecessary since so little time is comparatively spent
    evaluating the basis functions
 */
template<class ctype, int dim>
class MLinearShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;
  
public:
    //enum { D = dim };
  
  MLinearShapeFunction (unsigned long _mask=0) : mask(_mask)
  {
      //std::cout << "Creating shape function with mask " << mask << "\n";
  }

  ctype evaluateFunction (const coord_t& local) const
  {
    ctype r = 1.0;
    for (int i = 0; i < dim; ++i)
      r *= (mask & (1u<<i)) ? local[i] : (1-local[i]);
    
      //std::cout << "Evaluating shape #" << mask << " on " << local
      //          << " returns " << r << "\n";

    return r;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
    coord_t r;
    for (int i = 0; i < dim; ++i) {
      ctype t = 1.0;
        // evaluate partial derivative in i direction
        // For a node (1,0,1) the function is x(1-y)z and the gradient
        // ((1-y)z,x(-1)z,x(1-y))
      for (int j = 0; j < dim; ++j) {
        if (i == j) {
          t *= (mask & (1u<<j)) ? 1.0 : -1.0;
        } else
          t *= (mask & (1u<<j)) ? local[j] : (1-local[j]);
      }
      r[i] = t;
    }
    
      //cout << "Evaluating GRADIENT of shape #" << mask << " on " << local
      //     << " returns " << r << "\n";
    return r;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > ulong bitsize
};


/*! For the discrete space of Lagrange multipliers.
 
 These functions take value 2 on their associated vertex and -1 on the rest.
 
 FIXME!! Lots of mindless copying. I should also memoize. And this just sucks...
 */
  // FIXME: Why doesn't this work as I want it to?
template <class T>
std::vector<T>& operator<< (std::vector<T>& v, T d)
{
  v.push_back (d);
  return v;
}

template<class ctype, int dim>
class LagrangeSpaceShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;

  class Monome {
    int indices[dim];  // array (0 1 1) means x_2*x_3, etc.
    void clear() {
      for (int i=0; i < dim; ++i) indices[i] = 0;
    }
  public:
    ctype coeff;
    Monome (ctype _coeff=0) : coeff(_coeff) { clear(); }
    Monome (std::vector<int> ilist, ctype _coeff) : coeff (_coeff) {
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
      if (indices[which] == 0) {
        clear();
        coeff = 0;
      } else {
        indices[which] = 0;
      }
    }
    
    friend std::ostream& operator<< (std::ostream& os, const Monome& m)
    {
      os << m.coeff;
      for (int i=0; i < dim; ++i)
        if (m.indices[i] > 0)
          os << "x[" << i << "]";
      return os;
    }
  };
  
  typedef std::vector<Monome> Polynome;
  
  friend std::ostream& operator<< (std::ostream& os, const Polynome& p)
  {
    for (auto& m : p)
      os << m << " + ";
    return os;
  }
  
  friend Monome& operator* (Monome& m, ctype c)
  {
    m.coeff = c;
    return m;
  }
  
  friend Polynome& operator* (Polynome& p, ctype c)
  {
    for (auto& x : p)
      x = x*c;

    return p;
  }

  Polynome monomesOfOrder (int order)
  {
    /* TODO: translate this!! and delete the rest!!!!
     
     (define dim 3)

     (define (prepend x)
       (lambda (y) (append x y)))
     
     (define (monomes ord from prog)
       (cond ((>= from dim) '())
             ((== ord 1) (cons `(,from) (monomes ord (+ 1 from) prog)))
             (else
                (with np (append prog `(,from))
                  (append
                    (map (prepend np) (monomes (- ord 1) (+ 1 from) '()))
                    (monomes ord (+ 1 from) prog))))))
     */
    
    if (dim > 3)
      DUNE_THROW (Exception, "general monomesOfOrder not implemented in c++");

    /*!  MEGA-HACK!!
     
     This is just plain stupid, but it's too late now for recursion in c++
     without reasonable reference counting.
     */
    
    Polynome p;
    Monome m;
    std::vector<int> indices;
    if (order == 1) {
      for (int i=0; i < dim; ++i) {
        indices.push_back (i);
        m = Monome (indices, 1);
        p.push_back (m);
        indices.clear();
      }
    } else if (order == 2) {
      if (dim == 2) {
        indices << 0 << 1;
        m = Monome (indices, 1);
        p.push_back (m);
      } else if (dim == 3) {
        m = Monome (indices << 0 << 1, 1);
        p.push_back (m);
        indices.clear();
        m = Monome (indices << 0 << 2, 1);
        p.push_back (m);
        indices.clear();
        m = Monome (indices << 1 << 2, 1);
        p.push_back (m);
      }
    } else if (order == 3) {
      indices.clear();
      m = Monome (indices << 0 << 1 << 2, 1);
      p.push_back (m);

    }
    return p;
  }

public:
  Polynome poly;
  Polynome gradient[dim];
  
  LagrangeSpaceShapeFunction (unsigned long _mask=0) : mask(_mask)
  {
    std::cout << "Creating shape function with mask " << mask << "\n";

    Polynome p;
    for (int i=0; i < dim; ++i) {
      p = monomesOfOrder (i+1);
      int sign = (i%2 == 0) ? -1 : 1;
      p = p * (sign * 3.0);
      poly.insert (poly.end(), p.begin(), p.end());
    }
    
    Monome m(2);
    poly.push_back(m);
    
    for (int i=0; i < dim; ++i) {
      for (const auto& monome : poly) {
        m = Monome (monome);
        m.differentiate (i);
        gradient[i].push_back (m);
      }
    }
    
    cout << "   Shape has poly: " << poly << "\n" << "       and grad: ";
    for (int i=0; i < dim; ++i)
       cout << gradient[i] << ", ";
    cout << "\n";
    
  }
  
  /* The basis function at (0 1 1) is
   
     -3*(x_1 + (1-x_2) + (1-x_3)) + 3*(x_1*(1-x_2) + x_1*(1-x_3) + (1-x_2)*(1-x_3))
     -3*(x_1*(1-x_2)*(1-x_3)) + 2
   
   */  
  ctype evaluateFunction (const coord_t& local) const
  {
    ctype r = 0.0;
    coord_t idiotic_copy (local);
    for (int i = 0; i < dim; ++i)
      idiotic_copy[i] = (mask & (1u<<i)) ? 1-local[i] : local[i];
    
    for (auto& monome : poly)
      r += monome (idiotic_copy);

      //std::cout << "Evaluating shape #" << mask << " on " << local
      //      << " returns " << r << "\n";
    
    return r;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
    coord_t r;
    coord_t idiotic_copy (local);
    for (int i = 0; i < dim; ++i)
      idiotic_copy[i] = (mask & (1u<<i)) ? 1-local[i] : local[i];
    
    for (int i = 0; i < dim; ++i)
      for (auto& monome : gradient[i])
        r[i] += monome (idiotic_copy);

          //cout << "Evaluating GRADIENT of shape #" << mask << " on " << local
      //     << " returns " << r << "\n";
    return r;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > ulong bitsize
};


/*! Q1ShapeFunctionSet.
 
 Singleton collection of *linear* ShapeFunction
 */
template<class ctype, int dim, class ShapeFunction>
class Q1ShapeFunctionSet
{
  typedef FieldVector<ctype, dim> coord_t;

public:
  const GeometryType::BasicType basicType = GeometryType::cube;
  enum { N = 1 << dim};

  static const Q1ShapeFunctionSet& instance()
  {
    if (! _instance)
      _instance = new Q1ShapeFunctionSet();
    return *_instance;
  }

  inline const ShapeFunction& operator[] (int i) const
  {
    if (i < 0 || i >= N)
      DUNE_THROW (Exception, "Index out of bounds for shape function.");
    
    return *(f[i]);
  }
  
private:
  static Q1ShapeFunctionSet* _instance;
  ShapeFunction* f[N];
  
  /*! Convert from DUNE vertex indices in the reference element to our encoding.
   
   This provides the mapping between vertex numbers in DUNE:
   
   http://www.dune-project.org/doc/doxygen/dune-geometry-html/group___geometry_generic_reference_elements.html
   
   and our system which assigns to the vertex with coordinates (1,1,0), i.e.
   number 3 in DUNE the index 6 = 110b. We basically do a poor man's bit order
   reversal (lsb->msb) of dim bits.
   FIXME? breaks at dim = 16
   */
  inline unsigned int mapDuneIndex (int i)
  {
    unsigned int msb_i = 0;
    i = i & 0xFFFF;
    for (int b = 0; b < dim; ++b)
      msb_i |= ((i>>b)&1)<<(dim-1-b);
    
      //std::bitset<8> lsb(i);
      //std::bitset<8> msb(msb_i);
      //std::cout << "i= " << lsb << ",\t\tmsb_i= " << msb << "("
      //          << (int)msb_i << ")\n";
    return msb_i;
  }
  
  Q1ShapeFunctionSet & operator = (const Q1ShapeFunctionSet &);
  Q1ShapeFunctionSet (const Q1ShapeFunctionSet& other) { }
  Q1ShapeFunctionSet ()
  {
    for (int i=0; i < N; ++i)
      f[i] = new ShapeFunction (mapDuneIndex (i));
    atexit (this->atExit);
  }

  ~Q1ShapeFunctionSet ()
  {
    for (auto i=0; i < N; ++i)
      delete f[i];
  }
  
  static void atExit ()
  {
    delete _instance;
    _instance = 0;
  }
};

template <class C, int D, class TS>
Q1ShapeFunctionSet<C, D, TS>* Q1ShapeFunctionSet<C, D, TS>::_instance = 0;

#endif  // defined (SIGNORINI_SHAPEFUNCTIONS_HPP)
