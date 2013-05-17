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
//     cout << "ActiveSetFunctor::isSupported(" << i << ")= "
//          << multipliers[i] << " + c * (" << solution[i] << " - " << gap[i]
//          << ") = " << (*this)(i) << "\n";
   
    return (multipliers[i]+c*(solution[i]-gap[i])) > 0;
  }
};

/*! A functor which has support on some Physical Entity defined in Gmsh.
 
 class Evaluation must define operator() and export return_t.
 */
template <typename ctype, int dim, class TGF, class Evaluation>
class GmshBoundaryFunctor
{
  typedef shared_ptr<Evaluation> EvaluationPtr;
  const TGF&          gf;   //!< Grid factory
  std::vector<int> bi2pe;
  std::map<int, EvaluationPtr>* evals;  //!< Physical group --> Evaluation function
  std::string       name;

public:
  typedef FieldVector<ctype, dim> coord_t;
  typedef typename Evaluation::return_t return_t;
  static const int return_dim = Evaluation::return_dim;  
  static const int codim = 1;
  
  GmshBoundaryFunctor (const TGF& _gf,
                       const std::vector<int>& boundary_id_to_physical_entity,
                       std::map<int, EvaluationPtr>* _evals,
                       std::string _name="")
  : gf (_gf), bi2pe (boundary_id_to_physical_entity), evals (_evals), name (_name)
  { }
  
  ~GmshBoundaryFunctor () { /*delete evals;*/ }

  /*! We need the intersection to know which Evaluation is needed, otherwise
   a costly search would be necessary per each evaluation.
   NOTE that we don't perform any checks here!
   */
  template <class Intersection>
  return_t operator() (const Intersection& is, const coord_t& global) const
  {
    auto ev = (*evals) [bi2pe [is.boundarySegmentIndex()]];
    return (*ev)(global);
  }
  
  template <class Intersection>
  bool isSupported (const Intersection& is) const
  {
      // FIXME! adding: || !gf.wasInserted (is))  evaluates to false always?!??! (cube grids only??!)
    if (is.neighbor() || !is.boundary())
      return false;
      //    auto idx = gf.insertionIndex (is);  // not implemented for UGGrid
    auto idx = is.boundarySegmentIndex();       // but it's ok since indices are not reordered ("if not load balancing"?)
    if (!(idx > 0 && idx < bi2pe.size() && evals->find (bi2pe[idx]) != evals->end()))
      return false;
    return true;
  }
};


/*! A functor which has support on some Physical Entity defined in Gmsh.
 
 class Evaluation must define operator() and export return_t.
  */
template <typename ctype, int dim, class TGF, class Evaluation>
class GmshVolumeFunctor
{
//  typedef typename TGF::GridType GridType;
  typedef UGGrid<dim> GridType;
  typedef shared_ptr<Evaluation> EvaluationPtr;
  const TGF&          gf;   //!< Grid factory
  std::vector<int> ei2pe;
  std::map<int, EvaluationPtr>* evals;  //!< Physical group --> Evaluation function
  
  template <class E>
  class ExportEntityImplementation : public E
  {
  public:
    typedef typename E::Implementation Implementation;
  private:
    ExportEntityImplementation (const E& e) : E(e) { }
  };

public:
  typedef typename TGF::template Codim<0>::Entity Entity;
  typedef FieldVector<ctype, dim>                coord_t;
  typedef typename Evaluation::return_t         return_t;

  static const int return_dim = Evaluation::return_dim;
  static const int codim = 0;
  
  GmshVolumeFunctor (const TGF& _gf,
                     const std::vector<int>& element_id_to_physical_entity,
                     std::map<int, EvaluationPtr>* _evals)
  : gf (_gf), ei2pe (element_id_to_physical_entity), evals (_evals)
  { }
  
  ~GmshVolumeFunctor () { delete evals; }

  return_t operator() (const Entity& en, const coord_t& global) const
  {
      // There's no insertionIndex() in UGGrid
//    auto ev = (*evals) [ei2pe [gf.insertionIndex (en)]];
    
      // Maybe something like:
//    grid.levelIndexSet(0).index (en);
    
    auto ev = (*evals) [1];  //HACK
    return (*ev) (global);
  }
  
  bool isSupported (const Entity& en) const
  {
    auto idx = gf.insertionIndex (en);
    return (idx > 0 && idx < ei2pe.size() && evals->find (ei2pe[idx]) != evals->end());
  }
};


/*! Build complex constraints using chains of functors linked with AND and NOT.
 
 Constraints build a single linked list. The head of the list deletes all items
 when it dies.
 
 HACK: I don't know how to create a pointer to a generic Constraint object since
 this is a template: if I use Constraint* then I get this particular instantiation
 So I need the extra template parameter TNextConstraint. This sucks big time.
 
 REMOVE ME. This won't work, since the interfaces of the subsets of the boundary
 are nodes and we are checking whole Intersections...

template <typename coord_t, int dim, class TFunctor, class TNextConstraint>
class Constraint
{
public:
  typedef typename TFunctor::return_t return_t;

  static const int return_dim = TFunctor::return_dim;
  static const int codim = TFunctor::codim;
  
    // use the argument "affirmative" to negate the action of the functor
  Constraint (const TFunctor& _f, bool affirmative=true, TNextConstraint* _next=0)
    : func (_f), negate (!affirmative), next (_next) { };
  ~Constraint() { delete next; }
  
  return_t operator() (const coord_t& global) const { return func(global); }
  
  void insert (TNextConstraint* _next)
  {
    if (next) next->insert (_next);
    else      next = _next;
  }
  
  template <typename Entity>
  bool isSupported (const Entity& is) const
  {
    if (negate == func.isSupported (is))  // == is a xor for booleans!
      return false;
    if (next)
      return next->isSupported (is);
    return true;
  }

private:
  const TFunctor&  func;
  bool           negate;
  TNextConstraint* next;
};


struct DummyConstraint {
  template <typename Entity>
  bool isSupported(const Entity& e) const { return true; }
};
*/

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

#endif // FUNCTORS_HPP
