/******************************************************************************
 * functors.cpp                                                               *
 *                                                                            *
 * Collection of functors needed for a specific problem.                      *
 * Functors expect global coordinates for evaluation. They implement a method *
 * isSupported() accepting Dune::Geometry as argument and another a coord_t   *
 *                                                                            *
 * TODO:                                                                      *
 *   - Use data associated to a grid to define the support of any functor.    *
 *   - Tests                                                                  *
 *   - For fun: use policies in a Functor template class to specify how the   *
 *     support is defined (grid or manually), whether coordinates are local   *
 *     or global, ...)                                                        *
 *                                                                            *
 ******************************************************************************/
#ifndef FUNCTORS_HPP
#define FUNCTORS_HPP

#include <dune/common/fvector.hh>
#include <dune/common/fassign.hh>

#include "utils.hpp"

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
   
   In order to approximate ∂u_k/∂x_l we use that u a linear combination of the
   basis functions:
   
      u = u^a * phi^a  ==> ∂u_k/∂x_l = u^a_k * ∂phi^a/∂x_l
   */
  block_t stress (const coord_t& u, const coord_t& phi) const
  {
    block_t s (0.0);
    
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            s[k][l] += a[i][j][k][l] * u[k] * phi[l];

    return s;
  }
};

/*! A functor to model constant volume forces. */
template <typename ctype, int dim>
class VolumeLoad {
  typedef FieldVector<ctype, dim> coord_t;

public:
  coord_t operator() (const coord_t& x) const
  {
    coord_t ret;
    
      // Values from [FV05, p.36]
      //ret[0] = 9.6154; ret[1] = -4.8077;
      //return ret * 1.0e7;
    
    ret[0] = 0; ret[1] = -10;
    return ret * 1.0e8;
  }
};

/*! A boundary vectorial functor. Implements isSupported() */
template <typename ctype, int dim>
class Tractions {
  typedef FieldVector<ctype, dim> coord_t;

public:
  inline coord_t operator() (const coord_t& x) const
  {
    coord_t ret;

      // Values from [FV05, p.36]
    /*
    if (x[1] == 1.0) {
      ret[0] =  1.9231*x[0] - 3.8462;
      ret[1] = -13.462*x[0] + 2.8846;
    } else if (x[0] == 1.0) {
      ret[0] =  6.7308*x[1] - 5.7692;
      ret[1] = -3.8462*x[1] + 3.9231;
    } else {
        ret <<= zero;
    }
    
    return ret*1.0e7;
    */
    
    /* The values in Hüber&Wohlmuth2005
    if (isSupported(x)) ret <<= (0.5-x[0])*30, 6.5;  // 2D
    else                ret <<= zero;
     
    return ret;
     */
    
      // upper and lower tractions
    ret[0] = -2; ret[1] = -12;
    return ret * 1.0e7;
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
      //return (x[0] > 1 - x[1]);  // For the values from [FV05, p.36]
      //return x[0] > 0 && x[0] < 1;   // For the upper and lower tractions
    return x[0] > 0 && x[0] < 1 && x[1] > 0;  // For the upper tractions
  }
};


/*! A boundary scalar functor. */
template <typename ctype, int dim>
class NormalGap {  
  typedef FieldVector<ctype, dim> coord_t;

public:

  inline ctype operator() (const coord_t& x) const
  {
      // Careful! remember that it must be g(x) > 0
    return isSupported (x) ? sin (x[0]*M_PI) / 200.0 : 0;
  }
  
  template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
  bool isSupported (const class Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geo) const
  {
    for (int i = 0; i < geo.corners(); ++i)
      if (! isSupported (geo.corner(i))) {
          //cout << "Not supported at " << geo.corner(i) << "\n";
        return false;
      }

    return true;
  }
  
  inline bool isSupported (const coord_t& x) const
  {
    return (x[1] == 0 && x[0] > 0 && x[0] < 1);
  }
};

#endif // FUNCTORS_HPP
