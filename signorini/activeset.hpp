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
#include <dune/istl/matrixmatrix.hh>
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

  typedef typename TGV::Grid::GlobalIdSet      GlobalIdSet;
  typedef typename TGV::template Codim<dim>::Entity Vertex;
  typedef typename TGV::Grid::GlobalIdSet::IdType   IdType;
  typedef std::set<IdType>                           IdSet;
  typedef std::vector<const Vertex*>              IdVector;
  
  typedef typename TGV::ctype                        ctype;
  typedef FieldVector<ctype, dim>                  coord_t;
  typedef FieldMatrix<ctype, dim, dim>             block_t;
  typedef BCRSMatrix<block_t>                  BlockMatrix;
  typedef BCRSMatrix<coord_t>                 VectorMatrix;
  typedef BlockVector<coord_t>                 CoordVector;
  typedef BlockVector<FieldVector<ctype, 1> > ScalarVector;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector>  AIFunctor;
  typedef ActiveInactiveMapper<dim, TGV>               AIMapper;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                  MCMGVertexLayout> VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGT> GapVertexMapper;

private:
  const TGV& gv;    //!< Grid view
  const GlobalIdSet& gids; //

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
    //IdSet boundary;
  IdSet   active;
  IdSet inactive;
  IdSet   others;
    // tests
    //IdVector activeV;
    //IdVector inactiveV;
    //IdVector othersV;
  
  GapVertexMapper* gapMapper;  // FIXME: don't use gapMapper
  AIMapper*         aiMapper;

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
: gv (_gv), gids(_gv.grid().globalIdSet()), a (_a), f (_f), p (_p), gap(_gap),
  quadratureOrder(_quadratureOrder), gapMapper(new GapVertexMapper(_gv, _gap))
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
  
  g.resize (gapMapper->size());
  n_u.resize (gapMapper->size());
  n_m.resize (gapMapper->size());

  aiMapper = new AIMapper (gv, active, inactive, others);
}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::setupMatrices ()
{
  const int total = gv.size (dim);
  const int ingap = gapMapper->size();  // *Should* be active.size() + inactive.size()

  std::vector<std::set<int> > adjacencyPattern (total);
  
    //For each element we traverse all its vertices and set them as adjacent
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    int vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i) {
      for (int j = 0; j < vnum; ++j) {
        adjacencyPattern[aiMapper->map (*it, i, dim)].insert (aiMapper->map (*it, j, dim));
      }
    }
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
    D.setrowsize (i, 1);
  
  A.endrowsizes ();
  D.endrowsizes ();
  
  for (int i = 0; i < total; ++i)
    for (const auto& it : adjacencyPattern[i])
      A.addindex (i, it);
  
  for (int i = 0; i < ingap; ++i)
    D.addindex (i,i);
  

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
  bench().start ("Active set initialization");
  determineActive();
  bench().stop ("Active set initialization");
  
  bench().start ("Adjacency computation", false);
  setupMatrices();
  bench().stop ("Adjacency computation");

}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::determineActive ()
{
  AIFunctor contact (g, n_u, n_m, 1.0);  // 0 initial values *should* imply empty support
  
  active.clear();
  inactive.clear();
  others.clear();  // FIXME: no need to recalculate this!
  
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    auto id = gids.id (*it);
    if (gap.isSupported (it->geometry())) {
      if (contact.isSupported (gapMapper->map (*it))) {  // careful with the mapper!
        active << id;
      } else {
        inactive << id;
      }
    } else {
      others << id;
    }
  }
  
  aiMapper->update (active, inactive, others);
  
  cout << "Others:";
  for (auto x : others)   cout << " " << aiMapper->map(x);
  cout << "\nInactive:";
  for (auto x : inactive) cout << " " << aiMapper->map(x);
  cout << "\nActive:";
  for (auto x : active)   cout << " " << aiMapper->map(x);
  cout << "\n";
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
  int            offset = static_cast<int> (others.size());
  cout << "Offsetting by " << offset << "\n";
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        
        if (gap.isSupported (igeo)) {
          for (int i = 0 ; i < vnum; ++i) {
            int subi  = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids.subId (*it, subi, dim);
            int    ii = aiMapper->map (id) - offset;  // FIXME: the offset is a hack
            cout << "Setting index " << ii << "\n";
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (is->type(), quadratureOrder)) {
              const auto& global = igeo.global (x.position ());
              const auto&  local = it->geometry().local (global);
              D[ii][ii] = I * basis[subi].evaluateFunction (local) *
                          multBasis[subi].evaluateFunction (local) *
                          x.weight () *
                          igeo.integrationElement (x.position ());
            }
          }
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
            int   ii = aiMapper->map (*it, subi, dim);
            
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
  
  
    //// Recalculate matrices N, T
  
    // Ok, it's not exactly smart to use sparse matrices to store diagonal ones...

  const int K = (int)active.size();
  offset     += inactive.size();

  N.setSize (K, K);
  N.setBuildMode (VectorMatrix::random);
  T.setSize (K, K);
  T.setBuildMode (VectorMatrix::random);
  
  for (int i=0; i < K; ++i) {
    N.setrowsize (i, 1);
    T.setrowsize (i, 1);
  }
  
  N.endrowsizes();
  T.endrowsizes();
  
  AIFunctor contact (g, n_u, n_m, 1.0);  // 0 initial values *should* imply empty support
  
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    if (gap.isSupported (it->geometry())) {
      if (contact.isSupported (gapMapper->map (*it))) {  // careful with the mapper!
        int idx = aiMapper->map (*it);
        N.addindex (idx, idx);
        T.addindex (idx, idx);
      }
    }
  }
  
  N.endindices();
  T.endindices();
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        for (int i = 0 ; i < vnum; ++i) {
          int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          coord_t n = is->centerUnitOuterNormal ();  // ok?
          coord_t t; t[0] = -n[1]; t[1] = n[0];
          int   ii = aiMapper->map (*it, subi, dim) - offset;  // FIXME: the offset is a hack
          if (active.find (gids.subId (*it, subi, dim)) != active.end()) {  // ARGH!
            N[ii][ii] = n * D[ii][ii];
            T[ii][ii] = t * D[ii][ii];  // FIXMEMEMEMEEMME!!!!
          }
        }
      }
    }
  }
  
    //// Stiffness matrix:
  
  
  bench().report ("Assembly", "Feeding gremlins...", false);
  std::set<int> boundaryVisited;
  
  /*
   We compute the stiffness matrix in an unintuitive way: for each entry A_ij one
   must integrate (a function of) the gradients on all the domain (or rather on
   its very localized supports), but instead of traversing the whole grid for
   every entry, we traverse once and compute the contribution of each of the
   basis functions on each element to their corresponding matrix entry. This in
   turn uses the fact that we use P1 or Q1 elements with exactly as many nodes
   as vertices has an element (that's why we may use vnum)
   
   For each leaf we get a quadrature rule for the geometry type, then compute
   transformed gradients, obtain global indices of vertices i and j and update
   the associated matrix entry and right-hand side.
   */
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType typ = it->type ();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (typ);
    const auto&  geo = it->geometry();
    const int   vnum = ref.size (dim);
    
      /// Stiffness matrix
    
    for (auto& x : QuadratureRules<ctype, dim>::rule (typ, quadratureOrder)) {
      block_t jacInvTra = geo.jacobianInverseTransposed (x.position ());
      coord_t grad1, grad2;
      
      for (int i = 0; i < vnum; ++i) {
        jacInvTra.mv (basis[i].evaluateGradient (x.position ()), grad1);
        
        for (int j = 0; j < vnum; ++j) {
          jacInvTra.mv (basis[j].evaluateGradient (x.position ()), grad2);
          
          int ii = aiMapper->map (*it, i, dim);
          int jj = aiMapper->map (*it, j, dim);
            //cout << "ii= " << ii << ", jj= " << jj << "\n";
          try {
            A[ii][jj] += a (grad2, grad1) * x.weight () * geo.integrationElement (x.position ());
          } catch (ISTLError& e) {       // The adjacencyPattern does not match.
            cout << "FAILED setting data for A[" << ii << ", " << jj << "] = "
            << "(" << geo.corner(i) << ") x (" << geo.corner (j) << ")\n";
          }
        }
      }
        //// Integrand of the RHS
      
      for (int i = 0 ; i < vnum; ++i)
        b[aiMapper->map (*it, i, dim)] +=
            f (it->geometry ().global (x.position ())) *
            basis[i].evaluateFunction (x.position ()) *
            x.weight () *
            it->geometry ().integrationElement (x.position ());
      
    }
    
      //// Boundary conditions. See the discussion in the comments to the method
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int  ivnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        const auto& ityp = is->type ();
        
        if (p.isSupported (igeo)) {            // Neumann conditions.
                                               //cout << "Neumann'ing: "; printCorners (igeo);
          for (int i = 0 ; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            int   ii = aiMapper->map (*it, subi, dim);
              //auto   v = it->template subEntity<dim> (rsub)->geometry().center();
              //cout << "Neumann'ing node: " << ii << " at " << v << "\n";
            boundaryVisited.insert (ii);
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (ityp, quadratureOrder)) {
              auto global = igeo.global (x.position ());
              b[ii] += p (global) *
              basis[subi].evaluateFunction (it->geometry().local (global)) *
              x.weight () *
              igeo.integrationElement (x.position ());
            }
          }
        }
      }
    }
  }

  /* Dirichlet boundary conditions.
   
   For unvisited vertices on the boundary, replace the associated line of A
   and b with a trivial one.
   
   As stated before, this must be done in a separate step.
   */
  
  coord_t dirichlet;
  dirichlet <<= zero;
    //dirichlet[0] = 0; dirichlet[1] = -0.07;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
          //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
        
        for (int i = 0; i < ivnum; ++i) {
          auto subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            //auto v = it->template subEntity<dim> (rsub)->geometry().center();
          int ii = aiMapper->map (*it, subi , dim);
          if (boundaryVisited.count (ii) == 0) {
              //cout << "Dirichlet'ing node: " << ii << " at " << v << "\n";
            boundaryVisited.insert(ii);
            A[ii] = 0.0;
            A[ii][ii] = I;
            b[ii] = dirichlet;
          }
        }
      }
    }
  }
  
    //printmatrix (cout, A, "Stiffness matrix","");
  bench().report ("Assembly", "\tdone.");
  
}


template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGT, TSS, TLM>::solve ()
{
  bench().start ("Assembly");
  assemble();
  bench().stop ("Assembly");

    //// HACK of the HACKS...
  bench().start ("Gluing");
  const int n_N = (int)others.size();
  const int n_A = (int)active.size();
  const int n_I = (int)inactive.size();
  
  BlockMatrix B;
  B.setBuildMode (BlockMatrix::row_wise);
  B.setSize (n_N+n_I+2*n_A, n_N+n_I+n_A);
  CoordVector uu, c;
  c.resize (n_N+n_I+n_A, false);
  uu.resize (n_N+n_I+n_A, false);
  
  MatrixWindow<BlockMatrix> A_NN (A, 0, 0);
  MatrixWindow<BlockMatrix> A_NI (A, 0, n_N);
  MatrixWindow<BlockMatrix> A_NA (A, 0, n_N+n_I);
  MatrixWindow<BlockMatrix> A_IN (A, n_N, 0);
  MatrixWindow<BlockMatrix> A_II (A, n_N, n_N);
  MatrixWindow<BlockMatrix> A_IA (A, n_N, n_N+n_I);
  MatrixWindow<BlockMatrix> A_AN (A, n_N+n_I, 0);
  MatrixWindow<BlockMatrix> A_AI (A, n_N+n_I, n_N);
  MatrixWindow<BlockMatrix> A_AA (A, n_N+n_I, n_N+n_I);
  
  VectorMatrix T1(T);
    //matMultMat(T1, T, MatrixTransformer<BlockMatrix, VectorMatrix>(A_AN));

  VectorMatrix T2(T);
    //matMultMat(T2, T, MatrixTransformer<BlockMatrix, VectorMatrix>(A_AI));

  VectorMatrix T3(T);
    //matMultMat(T3, T, MatrixTransformer<BlockMatrix, VectorMatrix>(A_AA));
  
  for (auto row = B.createbegin(); row != B.createend(); ++row) {
    if (row.index() < n_N + n_I) {
      for (auto col = A[row.index()].begin(); col != A[row.index()].end(); ++col)
        row.insert (col.index());
    } else if (row.index() < n_N+n_I+n_A) {
      row.insert (row.index());
    } else if (row.index() < B.N()) {
      auto i = row.index() - n_N+n_I+n_A;
      for (auto col = T1[i].begin(); col != T1[i].end(); ++col)
        row.insert (col.index());
      for (auto col = T2[i].begin(); col != T2[i].end(); ++col)
        row.insert (col.index());
      for (auto col = T3[i].begin(); col != T3[i].end(); ++col)
        row.insert (col.index());
    }
  }

  for (auto row = B.begin(); row != B.end(); ++row) {
    if (row.index() < n_N+n_I) {
      for (auto col = A[row.index()].begin(); col != A[row.index()].end(); ++col)
        B[row.index()][col.index()] = A[row.index()][col.index()];
    } else if (row.index() < n_N+n_I+n_A) {
        //B[row.index()][row.index()] = N
    } else if (row.index() < B.N()) {
    }
  }
  
  bench().stop ("Gluing");
  
  bench().start ("Solving");
  InverseOperatorResult stats;
  MatrixAdapter<BlockMatrix, CoordVector, CoordVector> op (B);

   // initialize preconditioner and use bi-conjugate gradient solver
  SeqILUn<BlockMatrix, CoordVector, CoordVector> ilu1 (B, 1, 0.96);
  BiCGSTABSolver<CoordVector> bcgs (op, ilu1, 1e-15, 500, 0);
  bcgs.apply (uu, c, stats);
  bench().stop ("Solving");
}

#endif /* defined (SIGNORINI_ACTIVESET_HPP) */
