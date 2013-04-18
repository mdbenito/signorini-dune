/******************************************************************************
 * functors.cpp                                                               *
 *                                                                            *
 * Collection of functors needed for a specific problem.                      *
 * Functors expect global coordinates for evaluation. They implement a method *
 * isSupported() accepting Dune::Geometry as argument and another a coord_t   *
 *                                                                            *
 * TODO:                                                                      *
 *   - Use data associated to a grid to define the support of any functor.    *
 *   - Create a ProblemData template class to gather all functors, their type *
 *     information and the reading of the grid to define the support.         *
 *   - Tests                                                                  *
 *   - For fun: use policies in a Functor template class to specify how the   *
 *     support is defined (grid or manually), whether coordinates are local   *
 *     or global, ...)                                                        *
 ******************************************************************************/

#ifndef SIGNORINI_FUNCTORS_HPP
#define SIGNORINI_FUNCTORS_HPP

#include <dune/common/fvector.hh>
#include <dune/common/fassign.hh>

#include "utils.hpp"
#include <cmath>

/*! Kronecker's delta. */
template <typename T>
inline double kron (const T& i, const T& j)
{
  return (i==j) ? 1.0 : 0.0;
}

/*! Compute the action of the Hooke tensor on two vectors.
 
 This models a hyperelastic, isotropic and homogeneous material with Young's
 modulus and Poisson's ratio specified in the constructor.
 */
template <typename ctype, int dim>
class HookeTensor {
  typedef         FieldVector<ctype, dim> coord_t;
  typedef    FieldMatrix<ctype, dim, dim> block_t;
  typedef FieldMatrix<block_t, dim, dim> tensor_t;
  
  tensor_t a;

public:

  /*! Compute Lamé coefficients from Young's modulus and Poisson's ratio
   and initialize the coefficients of the tensor.
   */
  HookeTensor (double E=1.0, double nu=1.0)
  {
    double lambda = (E*nu)/((1+nu)*(1-2*nu));
    double     mu = E/(2*(1+nu));
    
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l) {
            a[i][j][k][l] = lambda * kron(i, j) * kron(k, l) +
                            mu * (kron(i,k)*kron(j,l) + kron(i,l)*kron(j,k));
          }
      //printmatrix(std::cout, a, "Hooke Tensor", "row", 12, 10);
  }
  
  block_t operator() (const coord_t& f, const coord_t& g) const
  {
    block_t r (0.0);
    
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            r[i][k] += a[i][j][k][l] * f[l] * g[j];
    
    return r;
  }
  
  /* Yes, we could return a vector here because the stress tensor
   is symmetric.
   
   In order to calculate ∂u_k/∂x_l we use that u is a linear combination of the
   basis functions:
   
      u = u^a * phi^a  ==> ∂u_k/∂x_l = u^a_k * ∂phi^a/∂x_l
   
   FIXME: Is this right?
   */
  block_t stress (const coord_t& u, const coord_t& phi) const
  {
    block_t s (0.0);
    
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            s[i][j] += a[i][j][k][l] * u[k] * phi[l];

    return s;
  }
};

/*! A functor to model constant volume forces. */
template <typename ctype, int dim>
class VolumeLoad {
  typedef FieldVector<ctype, dim> coord_t;
  coord_t ret;
public:
  VolumeLoad (const coord_t& _ret) : ret(_ret) { }

  coord_t operator() (const coord_t& x) const
  {
    return ret;
    
      //ret[0] = 9.6154; ret[1] = -4.8077; return ret * 1.0e7;  // [FV05, p.36]
      //ret <<= zero; return ret;                               // [HW05]
      //ret[0] = 0; ret[1] = -10; return ret * 1.0e7;           // DATA3
  }
};


/*! A functor to model Dirichlet data
 
 This just makes no sense. Some stuff hardcoded, some not...
 */
template <typename ctype, int dim>
class DirichletFunctor {
  typedef FieldVector<ctype, dim> coord_t;
  ctype cx, cy;
  coord_t ret;
public:
  DirichletFunctor (const coord_t& _ret) : ret(_ret) { }
  
  coord_t operator() (const coord_t& x) const
  {
    return ret;
  }
  
  template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
  bool isSupported (const class Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geo) const
  {
    for (int i = 0; i < geo.corners(); ++i)
      if (! isSupported (geo.corner (i)))
        return false;
    
    return true;
  }
  
  inline bool isSupported (const coord_t& x) const
  {
    return /* (x[0]>0 && x[0] < 1) && */ (x[1]==1.0 || x[1]<=-1.0);
  }

};


/*! A boundary vectorial functor. Implements isSupported()
 Hackish hackish
 */
template <typename ctype, int dim>
class Tractions {
  typedef FieldVector<ctype, dim> coord_t;
  coord_t ret;
public:
  Tractions (const coord_t& _ret) : ret(_ret) { }

  inline coord_t operator() (const coord_t& x) const
  {
      // Values from [FV05, p.36]
    /*
    if (x[1] == 1.0) {
      ret[0] =  1.9231*x[0] - 3.8462;  ret[1] = -13.462*x[0] + 2.8846;
    } else if (x[0] == 1.0) {
      ret[0] =  6.7308*x[1] - 5.7692;  ret[1] = -3.8462*x[1] + 3.9231;
    } else {
        ret <<= zero;
    }
    
    return ret*1.0e7;
    */
      //ret[0] = -2; ret[1] = -12; return ret * 1.0e7;       // DATA3
      //ret[0] = 0; ret[1] = -12; return ret * 1.0e7;       // DATA4
      //ret[0] = (0.5-x[0])*10; ret[1] = -12; return ret * 1.0e7;   // DATA5
    /* // [HW05]
    ret[0] = cx *(0.5-x[0]);
    ret[1] = cy;
     */
    
      // HACK
    if ((ret[0]>0 && x[0] == 1) || (ret[0]<0 && x[0] == 0))
      return coord_t(0.0);

    return ret;
  }

  template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
  bool isSupported (const class Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geo) const
  {
    for (int i = 0; i < geo.corners(); ++i)
      if (! isSupported (geo.corner (i)))
        return false;

    return true;
  }
  
  inline bool isSupported (const coord_t& x) const
  {
      //return (x[0] > 1 - x[1]);       // [FV05]
      //return x[0] > 0 && x[0] < 1;    // For the upper and lower tractions
      //return x[0] > 0 && x[0] < 1 && x[1] > 0;  // For the upper tractions
    return x[0] == 0 || x[0] == 1;   // [HW05]
  }
};


/*! A boundary scalar functor.
 
 WARNING! This is the *normalized* normal gap (d'oh!). It's not clear to me what
 the correct order of magnitude is or how to relate it with the size of the grid,
 etc.
 */
template <typename ctype, int dim>
class NormalGap {  
  typedef FieldVector<ctype, dim> coord_t;
  ctype ySupport;
  ctype ret;
public:
  NormalGap (ctype y=0.0, ctype r=0.05) : ySupport(y), ret(r) { }
  
    // Careful! remember that it must be g(x) > 0
  inline ctype operator() (const coord_t& x) const
  {
    return ret;                       // [HW05] uses 0.05
      //return sin (x[0]*6*M_PI) / 50.0;   // DATA3,4
      //return std::abs (sin (x[0]*6*M_PI) / 20.0);   // DATA5
  }
  
  template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
  bool isSupported (const class Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geo) const
  {
    for (int i = 0; i < geo.corners(); ++i)
      if (! isSupported (geo.corner(i)))
        return false;
    return true;
  }

  /* Why do I need the check with epsilon()? ctype is a double and literals are
   double by default, so if I construct with NormalGap(-0.05) the check here
   should work. But it does not!!!
   */
  inline bool isSupported (const coord_t& x) const
  {
      //return false;   // [FV05]
    bool ret = (std::abs(x[1] - ySupport) <= std::numeric_limits<ctype>::epsilon());
//    cout << "   isSupported(" << x << ") with ySupport= " << ySupport << ": "
//         << (ret ? "YES" : "NO") << LF;
    return ret;// && (x[0] != 0.0) && (x[0] != 1.0);
  }
};


/*! A boundary scalar functor for the active / inactive set strategy.
 */
template <class ctype, int dim, class ScalarVector>
class ActiveSetFunctor {
  typedef FieldVector<ctype, dim> coord_t;
  
  const ScalarVector& gap;
  const ScalarVector& solution;
  const ScalarVector& multipliers;
  ctype c;
  ctype d;
  
public:
    /// Parameter d achieves nothing, has no theoretical justification I know of
    /// and should be left alone as 1.0
  ActiveSetFunctor (const ScalarVector& _gap, const ScalarVector& _sol,
                    const ScalarVector& _mul, ctype _c=1.0, ctype _d=1.0)
  : gap (_gap), solution (_sol), multipliers (_mul), c (_c), d (_d) { }
  
  inline ctype operator() (int i) const
  {
    return d*multipliers[i]+c*(solution[i]-gap[i]);
  }
  
  inline bool isSupported (int i) const
  {
    /*
    cout << "ActiveSetFunctor::isSupported(" << i << ")= d * "
         << multipliers[i] << " + c * (" << solution[i] << " - " << gap[i]
         << ") = " << (*this)(i) << "\n";
     */
    return (d*multipliers[i]+c*(solution[i]-gap[i])) > 0;
  }
};



#endif // FUNCTORS_HPP
