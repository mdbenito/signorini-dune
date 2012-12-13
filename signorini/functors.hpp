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

/*! Kronecker's delta. */
template <typename T>
inline int kron (const T& i, const T& j)
{
  return (i==j) ? 1.0 : 0.0;
}

/*! Compute the action of the Hooke tensor on two vectors.
 
 This models a hyperelastic, isotropic and homogeneous material via Young's
 modulus and Poisson's ratio.
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
    block_t r(0.0);
    
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        for (int k = 0; k < dim; ++k)
          for (int l = 0; l < dim; ++l)
            r[i][k] += a[i][j][k][l] * f[l] * g[j];
    
    return r;
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
    ret <<= 9.6154e7, -4.8077e7;   // 2D
    return ret;
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
/* Tests: the commented lines are as in the paper by Figueiredo & Viaño.
    if (x[1] > 0.98 && x[0] > 0)
        //ret <<= (1.9231*x[0] - 3.8462)*1.0e7, (-13.462*x[0] + 2.8846)*1.0e7;
      ret <<= (-1.9231*x[0] + 3.8462)*1.0e7, (-13.462*x[0] + 2.8846)*1.0e7;
    else if (x[0] > 0.98)
        //ret <<= (6.7308*x[1] - 5.7692)*1.0e7, (1.9231 - 3.8462*x[1])*1.0e7;
      ret <<= (-1.9231*x[1] + 3.8462)*1.0e7, (-13.462*x[1] + 2.8846)*1.0e7;
    else
      ret <<= zero;
*/
    if (x[0] >= 1-x[1])
      ret <<= (-1.9231*x[1] + 3.8462)*1.0e7, (-13.462*x[1] + 2.8846)*1.0e7;
    else
      ret <<=  zero;
    
    /* The values in Hüber&Wohlmuth2005
    if (isSupported(x)) ret <<= (0.5-x[0])*30, 6.5;  // 2D
    else                ret <<= zero;
     */
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
    return (x[0] > 1 - x[1]);  // TEMPORARY
  }
};


/*! A boundary scalar functor. */
template <typename ctype, int dim>
class NormalGap {  
  typedef FieldVector<ctype, dim> coord_t;

public:
  inline ctype operator() (const coord_t& x) const
  {
    return isSupported(x) ? 0.005 : 0;  // Temporary. FIXME
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
    return (x[1] < 0.01 && x[0] != 0); // Temporary. FIXME
  }
};

#endif // FUNCTORS_HPP
