/******************************************************************************
 * activeset.hpp                                                              *
 ******************************************************************************/

#ifndef SIGNORINI_ACTIVESET_HPP
#define SIGNORINI_ACTIVESET_HPP

#include "config.h"
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fassign.hh>     // operator <<= for vectors
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include "utils.hpp"
#include "benchmark.hpp"
#include "shapefunctions.hpp"
#include "functorsupportmapper.hpp"
#include "activeinactivemapper.hpp"

using namespace Dune;
using std::cout;

/*! Solve the one body Signorini problem with an active-inactive set strategy.
 
 This implements the algorithm described in [LW04].
 
 Template type names:
    TGV: TGridView
    TET: TElasticityTensor
    TFF: TForcesFunctor
    TTF: TTractionsFunctor
    TGF: TGapFunctor
    TSS: TShapeSet
    TLM: TLagrangeMultipliersShapeSet
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
class SignoriniIASet
{
public:
  static const int dim = TGV::dimension;
  
  typedef typename TGV::EntityPointer        EntityPointer;
  typedef std::vector<const EntityPointer>   EntityPointerVector;

  typedef typename TGV::ctype                        ctype;
  typedef FieldVector<ctype, dim>                  coord_t;
  typedef FieldMatrix<ctype, dim, dim>             block_t;
  typedef BCRSMatrix<block_t>                  BlockMatrix;
  typedef BCRSMatrix<coord_t>                 VectorMatrix;
  typedef BlockVector<coord_t>                 CoordVector;
  typedef BlockVector<FieldVector<ctype, 1> > ScalarVector;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector>    AIFunctor;
  typedef ActiveInactiveMapper<dim, TGV> AIMapper;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                  MCMGVertexLayout> VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGT> GapVertexMapper;

private:
  const TGV& gv;    //!< Grid view
  
  const TET& a;     //!< Elasticity tensor
  const TFT& f;     //!< Volume forces
  const TTT& p;     //!< Boundary forces
  const TGT& gap;   //!< Normal gap function (scalar)
  
  block_t      I;   //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  BlockMatrix  A;   //!< Stiffness matrix
  BlockMatrix  D;   //!< See [HW04, eq. (3.5)]
  VectorMatrix N;   //!< See [HW04, eq. (3.5)]
  VectorMatrix T;   //!< See [HW04, eq. (3.5)]
  CoordVector  b;   //!< RHS: volume forces and tractions
  CoordVector  u;   //!< Solution
  CoordVector  m;   //!< Lagrange multiplier
  ScalarVector g;   //!< Gap functor evaluated at the gap nodes
  ScalarVector n_u; //!< Normal component of the solution at gap nodes
  ScalarVector n_m; //!< Normal component of the lagrange multiplier at gap nodes
  
  int          quadratureOrder;
  EntityPointerVector boundary;
  EntityPointerVector   active;
  EntityPointerVector inactive;
  EntityPointerVector   others;
  
  const GapVertexMapper& gapMapper;
  AIMapper&               aiMapper;

public:
  SignoriniIASet (const TGV& _gv,  const TET& _a, const TFT& _f, const TTT& _p,
                  const TGT& _gap, int _quadratureOrder = 4);
  
  void setupMatrices ();
  void initialize ();
  void assemble ();
  void determineActive ();
  void solve ();
  
  const CoordVector& solution() const { return u; }
};


/******************************************************************************
 * Implementation                                                             *
 ******************************************************************************/

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::SignoriniIASet (const TGV& _gv,
                                                                   const TET& _a,
                                                                   const TFT& _f,
                                                                   const TTT& _p,
                                                                   const TGT& _gap,
                                                                   int _quadratureOrder)
: gv (_gv), a (_a), f (_f), p (_p), gap(_gap), quadratureOrder(_quadratureOrder),
  gapMapper(_gv, _gap)
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::setupMatrices ()
{
  const int total = gv.size (dim);
  const int ingap = gapMapper.size();  // *Should* be active.size() + inactive.size()

  std::vector<std::set<int> > adjacencyPattern (total);
  std::vector<std::set<int> > adjacencyPattern2 (ingap);
  
    //For each element we traverse all its vertices and set them as adjacent
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    int vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i)
      for (int j = 0; j < vnum; ++j)
        adjacencyPattern[aiMapper.map (*it, i, dim)] .
            insert (aiMapper.map (*it, j, dim));
  }
  
    // D is diagonal. FIXME: should I use another data type?
  for (auto& x : active) {
    int idx = gapMapper.map (x) - others.size();
    adjacencyPattern2[idx].insert (idx);
  }

    //// Initialize default values

  A.setSize (total, total);
  A.setBuildMode (BlockMatrix::random);
  
  D.setSize (ingap, ingap);
  D.setBuildMode (BlockMatrix::random);

  b.resize (total, false);
  u.resize (total, false);
  
  for (int i = 0; i < total; ++i)
    A.setrowsize (i, adjacencyPattern[i].size ());
  
  for (int i = 0; i < ingap; ++i)
      D.setrowsize (i, adjacencyPattern2[i].size ());
  
  A.endrowsizes ();
  D.endrowsizes ();
  
  for (int i = 0; i < total; ++i)
    for (const auto& it : adjacencyPattern[i])
      A.addindex (i, it);
  
  for (int i = 0; i < ingap; ++i)
    for (const auto& it : adjacencyPattern2[i])
      D.addindex (i, it);

  A.endindices ();
  D.endindices ();
  
  A = 0.0;
  D = 0.0;
  b = 0.0;
  u = 0.0;
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::initialize ()
{
  g.resize (gapMapper.size());
  n_u.resize (gapMapper.size());
  n_m.resize (gapMapper.size());
  
  /*
  for (const auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it); ++is) {
      if (is->boundary ()) {
        boundary << *it;
        break;
      }
    }
  }
   */
    //  Determine boundary entities (codim 1)
  boundary.clear();
  for (const auto it = gv.template begin<1>(); it != gv.template end<1>(); ++it)
    if (it->boundary ())
      boundary << *it;
  cout << "The boundary contains " << boundary.size() << " entities.\n";
  
  bench().start ("Active set initialization");
  determineActive();
  bench().stop ("Active set initialization");
  
  bench().start ("Adjacency computation", false);
  setupMatrices();
  bench().stop ("Adjacency computation");

}

/*! Assemble the system matrix from the matrices A, D, N, T.

 FIXME: I should simply reorder matrices D, N, T. Something like:
 
   get the old iaMapper indices, and for each permutation old<->new, swap the
   corresponding entries in the matrix.

 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::assemble ()
{
    //// Recalculate matrix D.

    // Recall that the integral is over the gap boundary, so we need not integrate
    // the basis functions of nodes outside it.
  
  const auto& multBasis = TLM::instance();
  const auto&     basis = TSS::instance();
  const int      offset = others.size();
  
  for (const auto& e : boundary) {
    EntityPointer in = e.inside();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (in.type());
    const int   vnum = ref.size (e.indexInInside (), 1, dim);
    const auto& igeo = e.geometry ();
    
    if (gap.isSupported (igeo)) {
      for (int i = 0 ; i < vnum; ++i) {
        int subi = ref.subEntity (e.indexInInside (), 1, i, dim);
        int   ii = aiMapper.map (e, subi, dim) - offset;  // FIXME: the offset is a hack
        
        for (auto& x : QuadratureRules<ctype, dim-1>::rule (e.type(), quadratureOrder)) {
          auto& global = igeo.global (x.position ());
          auto&  local = e.geometry().local (global);
          D[ii][ii] += basis[subi].evaluateFunction (local) *
                       multBasis[subi].evaluateFunction (local) *
                       x.weight () *
                       igeo.integrationElement (x.position ());
        }
      }
    }
  }
  
  /*
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    GeometryType typ = it->type ();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (typ);
    const auto&  geo = it->geometry();
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int  ivnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        const auto& ityp = is->type ();
        
        if (gap.isSupported (igeo)) {
          for (int i = 0 ; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            int   ii = aiMapper.map (*it, subi, dim);
            
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (ityp, quadratureOrder)) {
              auto& global = igeo.global (x.position ());
              auto&  local = it->geometry().local (global);
              D[ii][ii] += basis[subi].evaluateFunction (local) *
              multBasis[subi].evaluateFunction (local) *
              x.weight () *
              igeo.integrationElement (x.position ());
            }
          }
        }
      }
    }
  }
  */
  
  
    /// Recalculate matrix N
  const int K = active.size();
  N.setSize (K, K);
  N.setBuildMode (BlockMatrix::random);

  for (const auto& e : boundary) {
    EntityPointer in = e.inside();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (in.type());
    const int   vnum = ref.size (e.indexInInside (), 1, dim);
    const auto& igeo = e.geometry ();
    coord_t n = e.unitOuterNormal ();
    for (int i = 0 ; i < vnum; ++i) {
      int subi = ref.subEntity (e.indexInInside (), 1, i, dim);
      int   ii = aiMapper.map (e, subi, dim) - offset;  // FIXME: the offset is a hack
    }
  }
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::determineActive ()
{
  AIFunctor contact (g, n_u, n_m);  // 0 initial values *should* imply empty support
  
  active.clear();
  inactive.clear();
  others.clear();  // FIXME: no need to recalculate this!
  
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    if (gap.isSupported (it->geometry())) {
      if (contact.isSupported (it->geometry()))
        active << *it;
      else
        inactive << *it;
    } else {
      others << *it;
    }
  }
  
  aiMapper.update (active, inactive, others);
  
  cout << "Others:";
  for (auto& x : others)   cout << " " << aiMapper.map(x);
  cout << "\nInactive:";
  for (auto& x : inactive) cout << " " << aiMapper.map(x);
  cout << "\nActive:";
  for (auto& x : active)   cout << " " << aiMapper.map(x);
  cout << "\n";
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::solve ()
{
  
}

#endif /* defined (SIGNORINI_ACTIVESET_HPP) */
