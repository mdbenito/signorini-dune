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
#include <dune/geometry/quadraturerules.hh>

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
  static const size_t N = dim;
  static const GeometryType::BasicType basicType = GeometryType::simplex;

  LinearShapeFunction () : indep(0.0), coeff(0.0) { }
  LinearShapeFunction (const coord_t& _coeff, const ctype _indep)
  : coeff (_coeff), indep (_indep) { }
  
  ctype evaluateFunction (const coord_t& local) const
  {
    ctype result = indep;
    for (int i = 0; i < dim; ++i)
      result += coeff[i] * local[i];
    return result;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
    return coeff;
  }
  
    // FIXME: HACK and WRONG
  inline bool isSupported (const coord_t& local) const
  {
    ctype s = 0.0;
    for (int i = 0; i < dim; ++i) {
      if (local[i] < 0 || local[i] > 1) return false;
      s += local[i];
    }
    return (s <= 1);
  }

  
private:
  coord_t coeff;
  ctype   indep;

  template<class TC, int d, int f, int i> friend class P1ShapeFunctionSet;
};


/*! P1ShapeFunctionSet
 
 Singleton collection of LinearShapeFunction for P1 elements. Scaled by "factor"
 and "indep".
 */
template<class ctype, int dim, int factor, int indep>
class P1ShapeFunctionSet
{
  typedef FieldVector<ctype, dim> coord_t;
  static const size_t N = dim + 1;
  
public:
  static const GeometryType::BasicType basicType = GeometryType::simplex;

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
    return *(f[i]);
  }
  
  static size_t size()
  {
    return N;
  }

private:

  P1ShapeFunctionSet()
  {
    f[0] = new ShapeFunction (coord_t(-1.0*factor), factor+indep);
    coord_t e;
    for (int i = 0; i < dim; ++i) {
      e = 0.0;
      e[i] = factor;
      f[i+1] = new ShapeFunction (e, indep);
    }
  }
  
  P1ShapeFunctionSet (const P1ShapeFunctionSet& other) {}
  
  ShapeFunction* f[dim+1];
};


/* 
 This is a complete "hut" basis function: the support extends to the four
 cuadrilaterals adjacent to the vertex where the value of the function is 1.
 
 FIXME! This is wrong! Plot the files output by testShapes() to see how the
 values and gradients are different of those of MLinearShapeFunction
 */
template<class ctype, int dim, int factor, int indep>
class BetterLinearShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;

    // Scalar function to construct the n-dimensional one
  ctype b0 (ctype x) const
  {
    if (-1 <= x && x <= 0)    return 1.0 + x;
    else if (0 < x && x <= 1) return 1.0 - x;
    else                      return     0.0;
  }

      // Derivative of the scalar function to construct the gradient
  ctype g0 (ctype x) const
  {
    if (-1 <= x && x < 0)      return  1.0;
    else if (0 <= x && x <= 1) return -1.0;
    else                       return  0.0;
  }
  
    /*! Returns true if this shape function has its apex at a node whose i-th
     coordinate is 1. */
  bool atCoord (int i) const
  {
    return (mask & (1u<<i));
  }
  
public:
  static const GeometryType::BasicType basicType = GeometryType::cube;
  
  BetterLinearShapeFunction (unsigned long _mask = 0) : mask (_mask)
  {
    cout << "BetterLinearShapeFunction needs FIXING! Aborting.\n";
    assert (false);
  }
  
  ctype evaluateFunction (const coord_t& local) const
  {
    if (! isSupported (local)) return 0.0;
    ctype r = factor;
    for (int i = 0; i < dim; ++i) {
      if (atCoord (i)) r *= b0 (local[i]-1);
      else             r *= b0 (local[i]);
    }
    r += indep;
      //std::cout << "Shape[" << mask << "] (" << local << ")= " << r << "\n";    
    return r;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
//    if (! isSupported (local)) return coord_t (0.0);
    coord_t r (factor);
    
    for (int i = 0; i < dim; ++i) {
      if (atCoord (i)) r[i] *= g0 (local[i]-1);
      else             r[i] *= g0 (local[i]);
      for (int j = 0; j < dim; ++j) {
        if (j == i) continue;
        if (atCoord (j)) r[i] *= b0 (local[j]-1);
        else             r[i] *= b0 (local[j]);
      }
    }
      //std::cout << "GRAD Shape[" << mask << "] (" << local << ")= " << r << "\n";
    return r;
  }
  
  inline bool isSupported (const coord_t& local) const
  {
    for (int i = 0; i < dim; ++i)
      if (std::abs ((atCoord (i) ? 1.0 : 0) - local[i]) > 1)
        return false;
    return true;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > bitsize of ulong
};

/*
 Basis elements for simplices are identified by their mask and an additional
 "artificial" index. If i,j,k are the 1st, 2nd and 3rd bits of the mask, then
 we define this index as l=1-i-j-k
 
template<class ctype, int dim, int factor, int indep, GeometryType::BasicType GeometryType::simplex>
class BetterLinearShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;
     
public:
  static const GeometryType::BasicType basicType = GeometryType::simplex;
  
  BetterLinearShapeFunction (unsigned long _mask = 0)
  : mask (_mask) { }
  
  ctype evaluateFunction (const coord_t& local) const
  {
    if (! isSupported (local)) return 0.0;
    ctype r = factor;

      // TODO
    
    r += indep;
      //std::cout << "Shape[" << mask << "] (" << local << ")= " << r << "\n";
    return r;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
    if (! isSupported (local)) return coord_t (0.0);
    coord_t r (factor);

      // TODO
    
      //std::cout << "GRAD Shape[" << mask << "] (" << local << ")= " << r << "\n";
    return r;
  }
  
  inline bool isSupported (const coord_t& local) const
  {
    for (int i = 0; i < dim; ++i) {
      ctype t = (mask & (1u<<i)) ? 1.0 : 0.0;
      if (t < 0 || t > 1)
        return false;
    }
    return true;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > bitsize of ulong
};

*/

/*! Q1ShapeFunctionSet.
 
 Singleton collection of linear ShapeFunctions
 */
template<class ctype, int dim, class ShapeFunction>
class Q1ShapeFunctionSet
{
  typedef FieldVector<ctype, dim> coord_t;

public:
  static const size_t N = (1 << dim);
  static const GeometryType::BasicType basicType = ShapeFunction::basicType;

  static const Q1ShapeFunctionSet& instance ()
  {
    if (! _instance)
      _instance = new Q1ShapeFunctionSet ();
    return *_instance;
  }

  inline const ShapeFunction& operator[] (int i) const
  {
    if (i < 0 || i >= N)
      DUNE_THROW (Exception, "Index out of bounds for shape function.");
    
    return *(f[i]);
  }
  
  static size_t size()
  {
    return N;
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
   */
  inline unsigned int mapDuneIndex (int i)
  {
    assert (dim < 16);
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
    assert (basicType == GeometryType::cube);
    std::cerr << "WARNING: *not* using mapDuneIndex() in Q1ShapeFunctionSet\n";
    for (int i=0; i < N; ++i) {
        /// WTF??!??!?!??! mapDuneIndex() is no longer necessary? was it ever?
        // when did I change the shape functions' evaluation?
      unsigned int mask = i; //mapDuneIndex(i);
//          cout << "Creating shape function " << i << " with mask " << mask << "\n";
      f[i] = new ShapeFunction (mask);
    }
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


/******************************************************************************
 *                                                                            *
 *                                    OLD CODE!                               *
 *                                                                            *
 ******************************************************************************/



/*! MLinearShapeFunction:
 
 Represents a shape function and provides methods to evaluate the function
 and its gradient on a cubic reference element.
 
 Any Lagrange shape function is identified by the coefficients of the vertex
 were it's 1. On the plane, vertex (0,0) has function (1-x)(1-y), vertex (0,1)
 function (1-x)y, vertex (1,1) function xy, and so on.
 
 With this scheme (product of factors being x_i or (1-x_i)) the functions are
 linear on the restriction to each edge originating at the vertex where they take
 the value one, since they are polynomes of order N with N-1 components fixed.
 
 The mask data member is actually just the coordinates of the vertex where this
 shape function is 1, i.e. we use the bits:
 
 mask = 5 = 0...0101 => vertex = (...,0,1,0,1)
 
 We extend the basis functions out of the reference element as "huts", i.e.
 flipping the slope:
 
   ^
  /|\
 / | \
/  |  \
 --*---*-
 ^   ^
 |   |
 element nodes
 
 This extension was previously used as a quick (and bogus) replacement for a proper
 projection map between the contact boundaries of the bodies and is now unused.
 
 
 TODO:
 - Test unit!!
 - Implement extension of the gradient
 - Memoization seems unnecessary since so little time is comparatively spent
 evaluating the basis functions
 */
template<class ctype, int dim>
class MLinearShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;

public:
  static const GeometryType::BasicType basicType = GeometryType::cube;
  
  MLinearShapeFunction (unsigned long _mask = 0) : mask (_mask) { }
  
  ctype evaluateFunction (const coord_t& local) const
  {
    if (!isSupported(local)) return 0.0;
    ctype r = 1.0;
    for (int i = 0; i < dim; ++i) {
      if (local[i]<0)
        r *= (mask & (1u<<i)) ? local[i] : (1.0 + local[i]);
      else if (local[i]>1)
        r *= 2.0 - local[i];
      else
        r *= (mask & (1u<<i)) ? local[i] : (1.0 - local[i]);
    }
    
      //std::cout << "Evaluating shape #" << mask << " on " << local
      //          << " returns " << r << "\n";
    
    return r;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
    coord_t r;

      // evaluate partial derivative in i direction
      // For a node (1,0,1) the function is x(1-y)z and the gradient
      // ((1-y)z,x(-1)z,x(1-y))
    for (int i = 0; i < dim; ++i) {
      ctype t = 1.0;
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
  
  inline bool isSupported (const coord_t& local) const
  {
    for (int i = 0; i < dim; ++i)
      if (std::abs(((mask & (1u<<i)) ? 1.0 : 0) - local[i]) > 1)
        return false;
    return true;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > bitsize of ulong
};


/*! For the discrete space of Lagrange multipliers.
 
 This works both for simplices and tetrahedra in 2 and 3 dimensions.
 
 These functions take value 2 on their associated vertex and -1 on the rest.
 
 Example: the basis function at (0 1 1) is
 
 -3*(x_1 + (1-x_2) + (1-x_3)) + 3*(x_1*(1-x_2) + x_1*(1-x_3) + (1-x_2)*(1-x_3))
 -3*(x_1*(1-x_2)*(1-x_3)) + 2
 
 FIXME!! Lots of mindless copying. I should also memoize. And this just sucks...
 */

template<class ctype, int dim>
class LagrangeSpaceShapeFunction
{
  typedef FieldVector<ctype, dim> coord_t;
  typedef Monomial<ctype, dim> Monomial;
  typedef typename Monomial::Polynomial Polynomial;
  
  Polynomial poly;
  Polynomial gradient[dim];
  
public:
  static const GeometryType::BasicType basicType = GeometryType::cube;
  
  LagrangeSpaceShapeFunction (unsigned long _mask = 0) : mask (_mask)
  {
    Polynomial p;
    for (int i = 1; i < dim; ++i) {
      p = Monomial::monomialsOfOrder (i);
      int sign = (i%2 == 0) ? -1 : 1;
      p = p * (sign * 3.0);
      poly.insert (poly.end(), p.begin(), p.end());
    }
    
    Monomial m(2);
    poly.push_back (m);
    
    for (int i=0; i < dim; ++i) {
      for (const auto& mon : poly) {
        m = Monomial (mon);
        m.differentiate (i);
        gradient[i].push_back (m);
      }
    }
    /*
     cout << "   Shape has poly: " << poly << "\n" << "       and grad: ";
     for (int i=0; i < dim; ++i)
     cout << gradient[i] << ", ";
     cout << "\n";
     */
  }
  
  ctype evaluateFunction (const coord_t& local) const
  {
    if (!isSupported(local)) return 0.0;
    
    coord_t c (local);
    /*
     ctype r=0.0;
     for (int i = 0; i < dim; ++i)
     copy[i] = (mask & (1u<<i)) ? (1 - local[i]) : local[i];
     
     for (auto& mon : poly)
     r += mon (copy);
     */
    
    switch (mask) {
      case 0:  // node x=0, y=0 (, z=0)
        if (c[0]<0) c[0] *= -1.0;
        if (c[1]<0) c[1] *= -1.0;
        if (dim == 3 && c[2]<0) c[2] *= -1.0;   // FIXME: guessing
        break;
      case 1:  // node x=1, y=0 (, z=0)
        if (c[0]>1) c[0] -= 1.0;
        else        c[0] = 1.0 - c[0];
        if (c[1]<0) c[1] *= -1.0;
        if (dim == 3 && c[2]<0) c[2] *= -1.0;   // FIXME: guessing
        break;
      case 2:  // node x=0, y=1 (, z=0)
        if (c[0]<0) c[0] *= -1.0;
        if (c[1]>1) c[1] -= 1.0;
        else        c[1] = 1.0 - c[1];
        if (dim == 3 && c[2]<0) c[2] *= -1.0;   // FIXME: guessing
        break;
      case 3:  // node x=1, y=1 (, z=0)
        if (c[0]>1) c[0] -= 1.0;
        else        c[0] = 1.0 - c[0];
        if (c[1]>1) c[1] -= 1.0;
        else        c[1] = 1.0 - c[1];
        if (dim == 3 && c[2]<0) c[2] *= -1.0;   // FIXME: guessing
        break;
      case 4:  // node x=0, y=0, z=1
        if (c[0]<0) c[0] *= -1.0;
        if (c[1]<0) c[1] *= -1.0;
        if (dim == 3) {    // FIXME: guessing
          if (c[2]>1) c[2] -= 1.0;
          else        c[2] = 1.0 - c[2];
        }
        break;
      case 5:  // node x=1, y=0, z=1
        if (c[0]>1) c[0] -= 1.0;
        else        c[0] = 1.0 - c[0];
        if (c[1]<0) c[1] *= -1.0;
        if (dim == 3) {    // FIXME: guessing
          if (c[2]>1) c[2] -= 1.0;
          else        c[2] = 1.0 - c[2];
        }
        break;
      case 6:  // node x=0, y=1, z=1
        if (c[0]<0) c[0] *= -1.0;
        if (c[1]>1) c[1] -= 1.0;
        else        c[1] = 1.0 - c[1];
        if (dim == 3) {    // FIXME: guessing
          if (c[2]>1) c[2] -= 1.0;
          else        c[2] = 1.0 - c[2];
        }
        break;
      case 7:  // node x=1, y=1, z=1
        if (c[0]>1) c[0] -= 1.0;
        else        c[0] = 1.0 - c[0];
        if (c[1]>1) c[1] -= 1.0;
        else        c[1] = 1.0 - c[1];
        if (dim == 3) {    // FIXME: guessing
          if (c[2]>1) c[2] -= 1.0;
          else        c[2] = 1.0 - c[2];
        }
        break;
      default:
        DUNE_THROW (Exception, "Invalid index for evaluateFunction");
    }
    if (dim==2)
      return -3.0 * c[0] - 3.0 * c[1] + 3.0 * c[0] * c[1] + 2.0;
    else if (dim==3)
      return - 3.0 * (c[0] + c[1] + c[2])
      + 3.0 * (c[0] * c[1] + c[1] * c[2] + c[0] * c[2])
      - 3.0 * (c[0] * c[1] * c[2])
      + 2.0;
  }
  
  inline coord_t evaluateGradient (const coord_t& local) const
  {
      //    if (!isSupported(local)) return coord_t(0.0);
    coord_t r (0.0);
    coord_t copy (local);
    for (int i = 0; i < dim; ++i)
      copy[i] = (mask & (1u<<i)) ? 1-local[i] : local[i];
    
    for (int i = 0; i < dim; ++i)
      for (auto& mon : gradient[i])
        r[i] += mon (copy);
    
    return r;
  }
  
  inline bool isSupported (const coord_t& local) const
  {
    for (int i = 0; i < dim; ++i)
      if (std::abs(((mask & (1u<<i)) ? 1.0 : 0) - local[i]) > 1)
        return false;
    return true;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > ulong bitsize
};
#endif