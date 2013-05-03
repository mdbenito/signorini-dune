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
#include "traverse.hpp"
#include <cmath>
#include <vector>
#include <set>

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


/*! A boundary scalar functor for the active / inactive set strategy.
 */
template <class ctype, int dim, class ScalarVector>
class ActiveSetFunctor {
  typedef FieldVector<ctype, dim> coord_t;
  
  const ScalarVector& gap;
  const ScalarVector& solution;
  const ScalarVector& multipliers;
  ctype c;
  
public:
  ActiveSetFunctor (const ScalarVector& _gap, const ScalarVector& _sol,
                    const ScalarVector& _mul, ctype _c=1.0)
  : gap (_gap), solution (_sol), multipliers (_mul), c (_c) { }
  
  inline ctype operator() (int i) const
  {
    return multipliers[i]+c*(solution[i]-gap[i]);
  }
  
  inline bool isSupported (int i) const
  {
    /*
     cout << "ActiveSetFunctor::isSupported(" << i << ")= "
     << multipliers[i] << " + c * (" << solution[i] << " - " << gap[i]
     << ") = " << (*this)(i) << "\n";
     */
    return (multipliers[i]+c*(solution[i]-gap[i])) > 0;
  }
};


  /// REMOVE ME!!! Use PSurface!
template <typename ctype, int dim>
class CylinderHackGapEvaluation {
  typedef FieldVector<ctype, dim> coord_t;
  typedef FieldVector<ctype, 1> scalar_t;
public:
  typedef scalar_t return_t;
  static const int return_dim = 1;
  
  CylinderHackGapEvaluation () { assert (dim == 3); }
  
    // Careful! remember that it must be g(x) >= 0
  return_t operator() (const coord_t& global) const
  {
    if (std::abs (global[0]) > 1.0 || std::abs (global[2]) > 1.0)
      return 100000.0;
    return 1.0 - std::sqrt (1.0 - global[0]*global[0]) +
           2.0 - std::sqrt (1.0 - global[2]*global[2]) - 1.0;
  }
};


/*! A quick hack to specify constraints to apply to GmshFunctors.
 
 I only need this until I understand what's wrong with the GmshReader and the
 Physical Entities.
 */
template <typename coord_t, int dim>
class Constraint {
  typedef bool (*Predicate)(const coord_t&);
  Predicate      f;
  bool      negate;
  Constraint* next;

public:
    // use the argument "affirmative" to negate the Predicate
  Constraint (Predicate _f=0, bool affirmative=true, Constraint* _next=0)
    : f (_f), negate (!affirmative), next (_next) { };
  
  bool at (const coord_t& global) const
  {
    if (f && (negate == f (global))) return false;  // == is a xor for booleans!
    if (next)            return next->at (global);
    return true;
  }
};


/*! An evaluation to model constant scalar data
 */
template <typename ctype, int dim, class ret_t>
class ConstantEvaluation {
  typedef FieldVector<ctype, dim> coord_t;
public:
  typedef ret_t return_t;
  static const int return_dim = return_t::dimension;
  
  ConstantEvaluation (const return_t& _ret) : ret (_ret) { }
  
  return_t operator() (const coord_t& global) const
  {
    (void) global;  // ignored!
    return ret;
  }

private:
  return_t ret;
};


/*! A functor which has support on some Physical Entity defined in Gmsh.
 
 class Evaluation must define operator() and export return_t.
 class Constraint must define at(). REMOVE THIS! it's just a hack!
 
 We take ownership of the Evaluation and Constraints objects and delete them
 when necessary.
 */
template <typename ctype, int dim, class TGF, class Evaluation, class Constraint>
class GmshBoundaryFunctor {
  typedef FieldVector<ctype, dim> coord_t;
  
  const TGF&          gf;   //!< Grid factory
  std::vector<int> bi2pe;
  std::set<int>   groups;   //!< List of physical groups managed by this functor
  Evaluation*       eval;
  Constraint* constraint;
  std::string       name;

public:
  typedef typename Evaluation::return_t return_t;
  static const int return_dim = Evaluation::return_dim;  
  static const int codim = 1;
  
  GmshBoundaryFunctor (const TGF& _gf,
                       const std::vector<int>& boundary_id_to_physical_entity,
                       const std::set<int>& _groups,
                       Evaluation* _eval,
                       Constraint* _cons=0,
                       std::string _name="")
  : gf (_gf), bi2pe (boundary_id_to_physical_entity), groups (_groups),
    eval (_eval), constraint (_cons), name (_name)
  { }
  
  ~GmshBoundaryFunctor () { delete eval; delete constraint; }
  
  return_t operator() (const coord_t& global) const
  {
    return (*eval) (global);
  }
  
  template <class Intersection>
  bool isSupported (const Intersection& is) const
  {
    if (is.neighbor() || !is.boundary())
           // FIXME! adding: || !gf.wasInserted (is))  evaluates to false always?!??! (cube grids only??!)
      return false;

      //    auto idx = gf.insertionIndex (is);  // not implemented for UGGrid
    auto idx = is.boundarySegmentIndex();       // but it's ok since indices are not reordered ("if not load balancing"?)
    if (!(idx > 0 && idx < bi2pe.size() && groups.find (bi2pe[idx]) != groups.end()))
      return false;

    if (constraint) {
      const auto&   in = is.inside();
      const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
      const int   vnum = ref.size (is.indexInInside(), 1, dim);
      for (int i = 0; i < vnum; ++i) {
        int  subi = ref.subEntity (is.indexInInside (), 1, i, dim);
        auto global = in->geometry().global (ref.position (subi, dim));
        if (! constraint->at (global)) {
            //        cout << "ConstraintHack " << name << " not fulfilled! (at " << global << ")\n";
          return false;
        }
      }
    }
    return true;
  }
};


/*! A functor which has support on some Physical Entity defined in Gmsh.
 
 class Evaluation must define operator() and export return_t.
 
 We take ownership of the Evaluation object and delete it when necessary.
 */
template <typename ctype, int dim, class TGF, class Evaluation>
class GmshVolumeFunctor {
  typedef FieldVector<ctype, dim>                coord_t;
  typedef typename TGF::template Codim<0>::Entity Entity;

  const TGF&          gf;   //!< Grid factory
  std::vector<int> ei2pe;
  std::set<int>   groups;   //!< List of physical groups managed by this functor
  Evaluation*       eval;
  
public:
  typedef typename Evaluation::return_t return_t;
  static const int return_dim = Evaluation::return_dim;
  static const int codim = 0;
  
  GmshVolumeFunctor (const TGF& _gf,
                     const std::vector<int>& element_id_to_physical_entity,
                     const std::set<int>& _groups,
                     Evaluation* _eval)
  : gf (_gf), ei2pe (element_id_to_physical_entity), groups (_groups), eval (_eval)
  { }
  
  ~GmshVolumeFunctor () { delete eval; }
  
  return_t operator() (const coord_t& global) const
  {
    return (*eval) (global);
  }
  
  bool isSupported (const Entity& en) const
  {
    auto idx = gf.insertionIndex (en);
    return (idx > 0 && idx < ei2pe.size() && groups.find (ei2pe[idx]) != groups.end());
  }
};

#endif // FUNCTORS_HPP
