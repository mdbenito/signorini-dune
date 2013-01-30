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
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
class SignoriniIASet
{
public:
  static const int dim = TGV::dimension;
  
  typedef            typename TGV::ctype ctype;
  typedef      FieldVector<ctype, dim> coord_t;
  typedef FieldMatrix<ctype, dim, dim> block_t;
  typedef      BCRSMatrix<block_t> BlockMatrix;
  typedef     BlockVector<coord_t> CoordVector;
  typedef BlockVector<FieldVector<ctype, 1> > ScalarVector;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                  MCMGVertexLayout> VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGT> GapVertexMapper;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector> AIFunctor;
  typedef ActiveInactiveMapper<dim, TGV, TGT, AIFunctor> AIMapper;

private:
  const TGV& gv;  //!< Grid view
  
  const TET& a;   //!< Elasticity tensor
  const TFT& f;   //!< Volume forces
  const TTT& p;   //!< Boundary forces
  const TGT& gap;   //!< Normal gap function (scalar)
  
  block_t      I;  //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  BlockMatrix  A;  //!< Stiffness matrix
  BlockMatrix  D;  //!< See [HW04, eq. (3.5)]
  CoordVector  b;  //!< RHS: volume forces and tractions
  CoordVector  u;  //!< Solution
  ScalarVector g;  //!< Gap functor evaluated at the gap nodes
  
  std::vector<std::set<int> > adjacencyPattern;
  int quadratureOrder;
  
public:
  SignoriniIASet (const TGV& _gv,  const TET& _a, const TFT& _f, const TTT& _p,
                  const TGT& _gap, int _quadratureOrder = 4);
  
  void setupMatrices ();
  void initialize ();
  void assembleMain ();
  void activateSet ();
  void solve ();
  
  const CoordVector& solution() const { return u; }
};


/******************************************************************************
 * Implementation                                                             *
 ******************************************************************************/

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::SignoriniIASet (const TGV& _gv,
                                                              const TET& _a,
                                                              const TFT& _f,
                                                              const TTT& _p,
                                                              const TGT& _gap,
                                                              int _quadratureOrder)
: gv (_gv), a (_a), f (_f), p (_p), gap(_gap), quadratureOrder(_quadratureOrder)
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::setupMatrices ()
{
  VertexMapper       mapper (gv.grid());
  GapVertexMapper gapMapper (gv, gap);
  
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    cout << "Vertex #" << mapper.map (*it)
         << " at " << it->geometry().center()
         << " has new index " << gapMapper.map (*it) << "\n";
  }
  
  g.resize (gapMapper.size());
  
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::initialize ()
{
  setupMatrices();
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::assembleMain ()
{
  
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::activateSet ()
{
  /*
   ScalarVector ngap (gapCount);
   ScalarVector nsol (gapCount);
   ScalarVector nmul (gapCount);
   AIFunctor activeset (ngap, nsol, nmul, c);
   
   */
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS>::solve ()
{
  
}

#endif /* defined (SIGNORINI_ACTIVESET_HPP) */
