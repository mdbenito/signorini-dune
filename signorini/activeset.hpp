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
#include <dune/istl/bdmatrix.hh>
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
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
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
  typedef FieldVector<ctype, 1>                   scalar_t;
  typedef FieldVector<ctype, dim>                  coord_t;
  typedef FieldMatrix<ctype, dim, dim>             block_t;
  typedef BCRSMatrix<scalar_t>                ScalarMatrix;
  typedef BCRSMatrix<coord_t>                 VectorMatrix;
  typedef BCRSMatrix<block_t>                  BlockMatrix;
  typedef BlockVector<scalar_t>               ScalarVector;
  typedef BlockVector<coord_t>                 CoordVector;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector>  AIFunctor;
  typedef ActiveInactiveMapper<dim, TGV>               AIMapper;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                  MCMGVertexLayout> VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGF> GapVertexMapper;

private:
  const TGV& gv;    //!< Grid view
  const GlobalIdSet& gids; //

  const TET& a;     //!< Elasticity tensor
  const TFT& f;     //!< Volume forces
  const TTT& p;     //!< Boundary forces
  const TGF& gap;   //!< Normal gap function (scalar)
  
  block_t      I;   //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  BlockMatrix  A;   //!< Stiffness matrix
  BlockMatrix  D;   //!< See [HW05] (Should be Dune::BDMatrix, but cannot build it)
  CoordVector  b;   //!< RHS: volume forces and tractions
  CoordVector  u;   //!< Solution
  CoordVector  n_d;   //!< D * normal at the gap nodes
  ScalarVector g;   //!< Gap functor evaluated at the gap nodes
  ScalarVector n_u; //!< Normal component of the solution at gap nodes
  ScalarVector n_m; //!< Normal component of the lagrange multiplier at gap nodes
  
  int          quadratureOrder;
    //IdSet boundary;
  IdSet   active;
  IdSet inactive;
  IdSet    other;
    // tests
    //IdVector activeV;
    //IdVector inactiveV;
    //IdVector otherV;

  AIMapper*         aiMapper;

public:
  SignoriniIASet (const TGV& _gv,  const TET& _a, const TFT& _f, const TTT& _p,
                  const TGF& _gap, int _quadratureOrder = 4);
  
  void setupMatrices ();
  void assemble ();
  void determineActive ();
  void step ();
  void solve ();

  const CoordVector& solution() const { return u; }
};


/******************************************************************************
 * Implementation                                                             *
 ******************************************************************************/

template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::SignoriniIASet (const TGV& _gv,
                                                                   const TET& _a,
                                                                   const TFT& _f,
                                                                   const TTT& _p,
                                                                   const TGF& _gap,
                                                                   int _quadratureOrder)
: gv (_gv), gids(_gv.grid().globalIdSet()), a (_a), f (_f), p (_p), gap(_gap),
  quadratureOrder(_quadratureOrder)
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;

    //// Initialize the AIMapper: set all nodes in the gap to inactive:
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    auto id = gids.id (*it);
    if (gap.isSupported (it->geometry()))  inactive << id;
    else                                      other << id;
  }
  aiMapper = new AIMapper (gv, active, inactive, other);
  
    //// Other initializations
  g.resize (inactive.size());
  n_u.resize (inactive.size());
  n_m.resize (inactive.size());
  n_u = 0;
  n_m = 0;
  g   = 0;
}

/*! Initializes the stiffness matrix.
 
 At step k=0, the system matrix (without boundary conditions) is
 
  /                 \
  | A_NN   A_NS   0 |
  | A_SN   A_SS   D |
  |    0      0   I |
  \                 /

 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::setupMatrices ()
{
  const int total = gv.size (dim);
  const int ingap = active.size() + inactive.size();
  cout << "total= " << total << ", ingap= " << ingap << "\n";

  std::vector<std::set<int> > adjacencyPattern (total+2*ingap); //too big actually
  
    //For each element we traverse all its vertices and set them as adjacent
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    int vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i) {
      int ii = aiMapper->map (*it, i, dim);
      for (int j = 0; j < vnum; ++j) {
        int jj = aiMapper->map (*it, j, dim);
        adjacencyPattern[ii].insert (jj);
      }
      /* Next we add one entry corresponding to the diagonal matrix D in
                /                 \
                | A_NN   A_NS   0 |
                | A_SN   A_SS   D |
                |    0      0   I |
                \                 /
       */
      if (other.find (gids.subId (*it, i, dim)) == other.end())  // i.e.: vertex is in S
        adjacencyPattern[ii].insert (total + aiMapper->mapInBoundary (*it, i, dim));
    }
  }

    //// Initialize default values

  A.setSize (total+ingap, total+ingap);
  A.setBuildMode (BlockMatrix::random);
  D.setSize (ingap, ingap);
  D.setBuildMode (BlockMatrix::random);
  
  b.resize (total+ingap, false);
  
  for (int i = 0; i < total; ++i)
    A.setrowsize (i, adjacencyPattern[i].size ());
  for (int i = total; i < total+ingap; ++i)
    A.setrowsize (i, 1);
  A.endrowsizes ();
  
  for (int i = 0; i < total; ++i)
    for (const auto& it : adjacencyPattern[i])
      A.addindex (i, it);
  for (int i = total; i < total+ingap; ++i)
    A.addindex (i, i);
  A.endindices ();
  
  for (int i = 0; i < ingap; ++i) D.setrowsize (i, 1);
  D.endrowsizes();
  for (int i = 0; i < ingap; ++i) D.addindex (i, i);
  D.endindices();
    //for (int i = 0; i < ingap; ++i) D[i][i] = I;
  
  /* This crashes:
  D.setBuildMode (BlockMatrix::row_wise);
  for (auto row = D.createbegin(); row != D.createend(); ++row) {
    auto i = row.index();
    row.insert (i);
    D[i][i] = I;
  }
  */

  A = 0.0;
  D = 0.0;
  b = 0.0;
}


template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::determineActive ()
{
  AIFunctor contact (g, n_u, n_m, 1.0);  // 0 initial values *should* imply empty support
  
  active.clear();
  inactive.clear();

  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    auto id = gids.id (*it);
    if (gap.isSupported (it->geometry())) {
      if (contact.isSupported (aiMapper->map (*it))) active << id;
      else                                         inactive << id;
    }
  }
  
  aiMapper->update (active, inactive, other);
  
    //cout << "Others:";
    //for (auto x : other)   cout << " " << aiMapper->map(x);
  cout << "\nInactive:";
  for (auto x : inactive) cout << x << "->" << aiMapper->map(x) << " ";
  cout << "\nActive:";
  for (auto x : active)   cout << x << "->" << " " << aiMapper->map(x) << " ";
  cout << "\n";
}

/*! Assemble the system matrix from the matrices A, D, N, T.

 FIXME: I should simply reorder matrices D, N, T. Something like:
 
   get the old iaMapper indices, and for each permutation old<->new, swap the
   corresponding entries in the matrix.

 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::assemble ()
{
  const auto& multBasis = TLM::instance();
  const auto&     basis = TSS::instance();
  const int       total = gv.size (dim);
  
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
    
      //// Gap at boundary for the computation of the active index set

    for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it)
      if (gap.isSupported (it->geometry())) {
        g[aiMapper->mapInBoundary (gids.id (*it))] = gap (it->geometry().center());
        boundaryVisited.insert (aiMapper->map (*it));
      }

      //// Neumann Boundary conditions.
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int  ivnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        const auto& ityp = is->type ();
        
        if (p.isSupported (igeo)) {
            //cout << "Neumann'ing: "; printCorners (igeo);
          for (int i = 0 ; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            int   ii = aiMapper->map (*it, subi, dim);
              //auto   v = it->template subEntity<dim> (subi)->geometry().center();
              //cout << "Neumann'ing node: " << ii << " at " << v << "\n";
            boundaryVisited.insert (ii);
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (ityp, quadratureOrder)) {
                // Transform relative (dim-1)-dimensional coord. in local coord.
              const auto& global = igeo.global (x.position ());
              const auto& local  = it->geometry().local (global);
              b[ii] += p (global) *
                      basis[subi].evaluateFunction (local) *
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
   and b with a trivial one. This must be done in a separate step.
   */
  coord_t dirichlet;
  dirichlet[0] = 0; dirichlet[1] = -0.07;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
          //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
        
        for (int i = 0; i < ivnum; ++i) {
          auto subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            //auto v = it->template subEntity<dim> (subi)->geometry().center();
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

    //// Calculate submatrix D.
  
    // Recall that the integral is over the gap boundary, so we need not integrate
    // the basis functions of nodes outside it.
  
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
            int    ii = aiMapper->map (id);
            
              //auto   v = it->template subEntity<dim> (subi)->geometry().center();
              //cout << "Setting index for D: " << ii << " at " << v << "\n";
            
            int kk = aiMapper->mapInBoundary (id);
            int jj = total + kk;
              //cout << "     jj= " << jj << ", kk= " << kk << "\n";
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (is->type(), quadratureOrder)) {
                // Transform relative (dim-1)-dimensional coord. in local coord.
              const auto& global = igeo.global (x.position ());
              const auto&  local = it->geometry().local (global);
              /*
              cout << "        quadrature point= " << global << " (" << local
                   << ")\n" << "           basis eval[" << subi
                   << "]= " << basis[subi].evaluateFunction (local)
                   << "\n" << "       multbasis eval[" << subi
                   << "]= " << multBasis[subi].evaluateFunction (local) << "\n";
               */
              A[ii][jj] += I * basis[subi].evaluateFunction (local) *
                           multBasis[subi].evaluateFunction (local) *
                           x.weight () *
                           igeo.integrationElement (x.position ());
              D[kk][kk] += A[ii][jj];
            }
          }
        }
      }
    }
  }

  bench().report ("Assembly", "\tdone.");
    //printmatrix (cout, D, "D", "");
}


template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::step ()
{
    //// HACK of the HACKS...
  bench().report ("Stepping", "Gluing", false);
  const int total = gv.size (dim);
  const int n_N = (int)other.size();
  const int n_A = (int)active.size();
  const int n_I = (int)inactive.size();
  
  CoordVector uu, c;
  n_d.resize (n_I+n_A, false);
  n_d = 0.0;
  c.resize (total+n_I+2*n_A, false);
  u.resize (total+n_A+n_I, false);
  uu.resize (total+n_A+n_I, false);
  uu = 0.0;

    // Copy RHS vector
  for (int i = 0; i < total; ++i)
    c[i] = b[i];
  for (int i = total; i < total+n_I+2*n_A; ++i)
    c[i] = 0.0; // we initialize the last n_A elements to gap() later.

  /* We copy matrix A and append the last lines. B should be:
   
         /                                 \
         | A_NN   A_NI   A_NS      0       |
         | A_IN   A_II   A_IA    D_I     0 |
         | A_AN   A_AI   A_AA      0   D_A |
     B = |    0      0      0   Id_I     0 |
         |    0      0      0      0   T_A |
         |    0      0    N_A      0     0 |
         \                                /
   
   The line | 0 0 N_A 0 0 | has row indices in [total+n_I+n_A, total+n_I+2*n_A)
   
   */
  BlockMatrix B;
  B.setBuildMode (BlockMatrix::row_wise);
    //B.setSize (total+n_I+2*n_A, total+n_I+2*n_A);
    //cout << "B is " << total+n_I+2*n_A << " x " << total+n_I+2*n_A << "\n";
  B.setSize (total+n_I+2*n_A, total+n_I+n_A);
  cout << "B is " << total+n_I+2*n_A << " x " << total+n_I+n_A << "\n";
  
    // Build an adjacency pattern for the last block-line in B
  std::vector<std::set<int> > adjacencyPattern (total);
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    int vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i) {
      IdType id = gids.subId (*it, i, dim);
      if (active.find (id) != active.end()) {
        int ii = aiMapper->mapInActive (id);
        int jj = aiMapper->map (id);
        adjacencyPattern[ii].insert (jj);
      }
    }
  }
  
  for (auto row = B.createbegin(); row != B.createend(); ++row) {
    if (row.index() < total) {
      for (auto col = A[row.index()].begin(); col != A[row.index()].end(); ++col)
        row.insert (col.index());
    } else if (row.index() < total+n_I) {
      row.insert (row.index());
    } else if (row.index() < total+n_I+n_A) {
      row.insert (row.index());
        //row.insert (row.index()+n_A);  // HACK: zero columns to have a square matrix
    } else if (row.index() < B.N()) {
      for (const auto& it : adjacencyPattern[row.index() - total - n_I - n_A])
        row.insert (it);
    }
  }
  
  B = 0.0;
  
  for (auto row = B.begin(); row != B.end(); ++row) {
    if (row.index() < total) {
      for (auto col = A[row.index()].begin(); col != A[row.index()].end(); ++col)
        B[row.index()][col.index()] = A[row.index()][col.index()];
    } else if (row.index() < total+n_I) {
      B[row.index()][row.index()] = I;
    } else if (row.index() < total+n_I+n_A) {
      B[row.index()][row.index()] = I;
        //B[row.index()][row.index()+n_A] = 0;// HACK: zero columns to have a square matrix
    }
  }
  
    // Fill the lines:
    //             |    0      0      0      0   T_A |
    //             |    0      0    N_A      0     0 |
  
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
            int kk = aiMapper->mapInBoundary (id);
            coord_t nr = FMatrixHelp::mult (D[kk][kk], is->centerUnitOuterNormal());
              //cout << "kk= " << kk << ", nr= " << nr << "\n";
            n_d[kk] += nr*(1.0/vnum); // FIXME: we shouldn't compute this here
            
            if (active.find (id) != active.end()) {
                //cout << "Found active: " << id << " -> " << kk << "\n";
              auto ipos = is->inside()->template subEntity<dim>(subi)->geometry().center();

                // HACK: Build temporary matrices with repeated elements.
                // We could flatten the whole matrix B instead to scalar blocks.
                // The problem is that some rows must be matrices, some vectors.
              
              coord_t tg;
              tg[0] = -nr[1]; tg[1] = nr[0];
              for (int r=2; r<dim; ++r)
                tg[r] = 0;
              
              block_t T, N;
              for (int k=0; k<dim; ++k) {
                for (int l=0; l<dim; ++l) {
                  T[k][l] = tg[l];
                  N[k][l] = nr[l];
                }
              }
                // first the tangential stuff for the multipliers.
              int ii = total + kk;
              B[ii][ii] += T*(1.0/vnum);
              c[ii] = 0.0;

                // now the last rows
              ii = total + n_A + n_I + aiMapper->mapInActive (id);
              int jj = aiMapper->map (id);
                //cout << "ii= " << ii << ", jj= " << jj << "\n";
              B[ii][jj] += N*(1.0/vnum);
              c[ii] = gap(ipos);  // copies the scalar dim times (part of HACK)
            }
          }
        }
      }
    }
  }

  bench().report ("Stepping", " done.");

    //printmatrix(cout, B, "B", "");
    //printvector(cout, c, "c", "");
  writeMatrixToMatlab (B, "/tmp/stiff");
  writeVectorToFile (c, "/tmp/rhs");

  bench().report ("Stepping", "Solving", false);
  
  InverseOperatorResult stats;
  MatrixAdapter<BlockMatrix, CoordVector, CoordVector> op (B);

  SeqSSOR<BlockMatrix, CoordVector, CoordVector> ssor (B, 1, 0.95);
    //BiCGSTABSolver<CoordVector> bcgs (op, ssor, 1e-8, 100, 0);
  RestartedGMResSolver<CoordVector> rgmres (op, ssor, 1e-8, gv.size(dim)/2, 200, 0);
  rgmres.apply (uu, c, stats);
  
  /*
    // This one seems to work
  SeqSSOR<BlockMatrix, CoordVector, CoordVector> ssor (B, 1, 0.96);
  LoopSolver<CoordVector> lsol (op, ssor, 1e-8, 200, 1);
  lsol.apply (uu, c, stats);
  */

  /*
   //The SeqILUn preconditioner stalls as soon as any node is active
  SeqILUn<BlockMatrix, CoordVector, CoordVector> ilu1 (B, 1, 0.96);
  BiCGSTABSolver<CoordVector> bcgs (op, ilu1, 1e-8, 200, 0);
  bcgs.apply (uu, c, stats);
   */

  /*
    //This solver produces lots of NaNs as soon as any node is active
  SeqSSOR<BlockMatrix, CoordVector, CoordVector> ssor (B, 1, 0.96);
  CGSolver<CoordVector> cgs (op, ssor, 1e-8, 200, 2);
  cgs.apply (uu, c, stats);
  */
  
  bench().report ("Stepping", " done.");

    //// FIXME: the following stuff doesn't belong here.
  
    // Gap at boundary for the computation of the active index set
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it)
    if (gap.isSupported (it->geometry()))
      g[aiMapper->mapInBoundary (gids.id (*it))] = gap (it->geometry().center());
  
  u = uu;
  for (const auto& x : inactive) {
    int i = aiMapper->mapInBoundary(x);
    int j = aiMapper->map(x);
    n_u[i] = n_d[i] * uu[j];
    n_m[i] = n_d[i] * uu[total+i];
  }
  for (const auto& x : active) {
    int i = aiMapper->mapInBoundary(x);
    int j = aiMapper->map(x);
    n_u[i] = n_d[i] * uu[j];
    n_m[i] = n_d[i] * uu[total+i];
  }
  
  /* This must be wrong:
  for (int i = 0; i < n_I+n_A; ++i) {
    n_u[i] = n_d[i] * uu[n_N+i];
    n_m[i] = n_d[i] * uu[total+i];
  }
   */
    //printvector (cout, n_d, "D * Normal", "");
  printvector (cout, n_u, "Solution Normal", "");
  printvector (cout, n_m, "Multiplier Normal", "");
}

template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::solve ()
{
  PostProcessor<TGV, TET, AIMapper, TSS> post (gv, *aiMapper, a);

  for (int cnt = 0; cnt < 5; ++cnt) {
    bench().start ("Active set initialization", false);
    determineActive();
    bench().stop ("Active set initialization");
    bench().start ("Adjacency computation", false);
    setupMatrices();
    bench().stop ("Adjacency computation");
    bench().start ("Assembly");
    assemble();
    bench().stop ("Assembly");
    bench().start ("Stepping");
    step();
    bench().stop ("Stepping");
    bench().start ("Postprocessing", false);
    /*
    double error = post.computeError (u);
    post.computeVonMisesSquared ();
    (void) post.writeVTKFile ("/tmp/SignoriniFEM", cnt);
     */
    bench().stop ("Postprocessing");
  }
}
#endif /* defined (SIGNORINI_ACTIVESET_HPP) */
