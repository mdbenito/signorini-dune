#ifndef SHAPEFUNCTIONS_HPP
#define SHAPEFUNCTIONS_HPP

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

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
  enum { D = dim };
  
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
    // private constructor prevents additional instances
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
  enum { D = dim };
  
  MLinearShapeFunction (unsigned long _mask=0) : mask(_mask)
  {
      //std::cout << "Creating shape function with mask " << mask << "\n";
  }

  ctype evaluateFunction (const coord_t& local) const
  {
    
    ctype r = 1.0;
    for (int i = 0; i < dim; ++i)
      r *= (mask & (1u<<i)) ? local[i] : (1-local[i]);
    
      //std::cout << "Evaluating shape #" << mask << " on " << local << " returns " << r << "\n";

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
    
      //cout << "Evaluating GRADIENT of shape #" << mask << " on " << local << " returns " << r << "\n";
    return r;
  }
  
private:
  unsigned long mask;  //!< Would break for dim > ulong bitsize
};



/*! Q1ShapeFunctionSet
 
 Singleton collection of LinearShapeFunction for Q1 elements.
 */
template<class ctype, int dim>
class Q1ShapeFunctionSet
{
  typedef FieldVector<ctype, dim> coord_t;

public:
  enum { N = 1 << dim};
  
  typedef MLinearShapeFunction<ctype, dim> ShapeFunction;
  
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
  
  Q1ShapeFunctionSet & operator = (const Q1ShapeFunctionSet &);
  Q1ShapeFunctionSet (const Q1ShapeFunctionSet& other) { }
  Q1ShapeFunctionSet ()
  {
    cout << "Initialiting Q1ShapeFunctionSet of size " << N << "\n";
    for (unsigned long i=0; i < N; ++i)  // overkill type
      f[i] = new ShapeFunction (i);
      //FIXME!
      //atexit (&Q1ShapeFunctionSet::atExit);
  }

  ~Q1ShapeFunctionSet ()
  {
    for (auto i=0; i < N; ++i)  // overkill type
      delete f[i];
  }
  
  void atExit ()
  {
    delete _instance;
    _instance = 0;
  }
};

template <class C, int D> Q1ShapeFunctionSet<C,D>*
Q1ShapeFunctionSet<C,D>::_instance = 0;


#endif  //SHAPEFUNCTIONS_HPP
