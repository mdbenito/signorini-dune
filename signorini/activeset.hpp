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
 
 This implements the algorithm described in [HW05].
 
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
  typedef FieldMatrix<ctype, 1, 1>             scalarmat_t;
  typedef BCRSMatrix<scalarmat_t>             ScalarMatrix;
    //typedef BCRSMatrix<coord_t>                 VectorMatrix; // WRONG
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
  CoordVector  n_d; //!< D * normal at the gap nodes
  ScalarVector g;   //!< Gap functor evaluated at the gap nodes
  ScalarVector n_u; //!< Normal component of the solution at gap nodes
  ScalarVector n_m; //!< Normal component of the lagrange multiplier at gap nodes
  
  int quadratureOrder;
  IdSet        active;
  IdSet      inactive;
  IdSet         other;
  AIMapper*  aiMapper;

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

    //// Initialize the AIMapper setting all nodes in the gap to inactive:
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    auto id = gids.id (*it);
    if (gap.isSupported (it->geometry()))  inactive << id;
    else                                      other << id;
  }
  aiMapper = new AIMapper (gv, active, inactive, other);
  
    //// Other initializations (inactive.size() at init =  size of gap)
  g.resize (inactive.size());
  n_d.resize (inactive.size());
  n_u.resize (inactive.size());
  n_m.resize (inactive.size());
  n_d = 0.0;
  n_u = 0.0;
  n_m = 0.0;
  g   = 0.0;
}

/*! Initializes the stiffness matrix.
 
 At step k=0, the system matrix (without boundary conditions) is
 
  /             \
  | A_NN   A_NS |
  | A_SN   A_SS |
  \             /

 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::setupMatrices ()
{
  const int  n_T = gv.size (dim);
  const auto ingap = active.size() + inactive.size();
    //cout << "total= " << n_T << ", ingap= " << ingap << "\n";

  std::vector<std::set<int> > adjacencyPattern (n_T);
  
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
    }
  }

    //// Initialize default values

  A.setSize (n_T, n_T);
  A.setBuildMode (BlockMatrix::random);
  D.setSize (ingap, ingap);
  D.setBuildMode (BlockMatrix::random);
  
  b.resize (n_T + inactive.size() + active.size(), false);
  u.resize (n_T + inactive.size() + active.size(), false);

  for (int i = 0; i < n_T; ++i)
    A.setrowsize (i, adjacencyPattern[i].size ());
  A.endrowsizes ();
  
  for (int i = 0; i < n_T; ++i)
    for (const auto& it : adjacencyPattern[i])
      A.addindex (i, it);
  A.endindices ();
  
  for (int i = 0; i < ingap; ++i)
    D.setrowsize (i, 1);
  D.endrowsizes();
  for (int i = 0; i < ingap; ++i)
    D.addindex (i, i);
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
  u = 0.0;
}


/*!
 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::determineActive ()
{
  AIFunctor contact (g, n_u, n_m, 1.0); // all nodes inactive for initial values == 0
  
  active.clear();
  inactive.clear();

  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    auto id = gids.id (*it);
    if (gap.isSupported (it->geometry())) {
      if (contact.isSupported (aiMapper->mapInBoundary (id)))
        active << id;
      else
        inactive << id;
    }
  }
  
  aiMapper->update (active, inactive, other);
  
  cout << "\nInactive: " << inactive.size();
    //for (auto x : inactive) cout << x << "->" << aiMapper->map(x) << " ";
  cout << "\nActive: " << active.size();
    //for (auto x : active)   cout << x << "->" << " " << aiMapper->map(x) << " ";
  cout << "\n";
}


/*! Assemble the system matrix from the matrices A, D.
 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::assemble ()
{
  const auto& multBasis = TLM::instance();
  const auto&     basis = TSS::instance();
  
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
   Replace the associated line of A and b with a trivial one.
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
          auto v = it->template subEntity<dim> (subi)->geometry().center();
          int ii = aiMapper->map (*it, subi , dim);
          if (v[1]==1 && v[0]>0 && v[0]<1) {  // HACK! replace with functor!!
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
  
  g = 0.0;
    //// Calculate submatrix D and  the gap at boundary for the computation of
    //// the active index set
    
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
            int kk = aiMapper->mapInBoundary (id);
              //auto   v = it->template subEntity<dim> (subi)->geometry().center();
              //cout << "Setting index for D: " << kk << " at " << v << "\n";
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
              D[kk][kk] += I * basis[subi].evaluateFunction (local) *
                           multBasis[subi].evaluateFunction (local) *
                           x.weight () *
                           igeo.integrationElement (x.position ());
              g[kk] += gap (global) *
                       multBasis[subi].evaluateFunction (local) *
                       x.weight () *
                       igeo.integrationElement (x.position ());
              boundaryVisited.insert (aiMapper->map (id));
            }
          }
        }
      }
    }
  }
 
  bench().report ("Assembly", "\tdone.");
    //printmatrix (cout, D, "D", "");
}

/*! HACK of the HACKS...
 */
template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::step ()
{

  bench().report ("Stepping", "Gluing", false);
  const int n_T = gv.size (dim);
  const int n_N = (int)other.size();
  const int n_A = (int)active.size();
  const int n_I = (int)inactive.size();
  const int total = n_T + n_I + n_A;
  
  assert (n_T == n_N + n_A + n_I);

  ScalarVector uu, c;
  c.resize (total*dim, false);
  uu.resize (total*dim, false);
  uu  = 0.0;
  n_d = 0.0;
  n_m = 0.0;
  n_u = 0.0;

    // Copy RHS vector
  for (int i = 0; i < n_T*dim; ++i)
    for (int j = 0; j < dim; ++j)
      c[i*dim+j] = b[i][j];
  for (int i = n_T*dim; i < total*dim; ++i)
    c[i] = 0.0; // we initialize the last n_A elements to gap() later.

  /*
   We copy matrix A and append some columns and lines. B should be
   
         /                                 \
         | A_NN   A_NI   A_NS      0       |
         | A_IN   A_II   A_IA    D_I     0 |
         | A_AN   A_AI   A_AA      0   D_A |
     B = |    0      0      0   Id_I     0 |
         |    0      0      0      0   T_A |
         |    0      0    N_A      0     0 |
         \                                /
   
   In TWO DIMENSIONS THIS IS SQUARE as a scalar matrix!!!!
   */

  ScalarMatrix B;
  B.setBuildMode (ScalarMatrix::row_wise);
  B.setSize (total*dim, total*dim);
    //cout << " B is " << B.N() << " x " << B.M() << "\n";

    // Flatten the adjacency pattern of A to scalar entries
  std::vector<std::set<int> > adjacencyPattern (total*dim);
  for (auto row = A.begin(); row != A.end(); ++row) {
    for (auto col = (*row).begin(); col != (*row).end(); ++col) {
      auto r = row.index();
      auto c = col.index();
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          if (A[r][c][i][j] != 0.0)
            adjacencyPattern[r*dim+i].insert((int)c*dim+j);
    }
  }
  
    // Build the adjacency pattern for the entries in B to the right of and below A
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    IdType id = gids.id (*it);
    if (inactive.find (id) != inactive.end()) {
      int ii = aiMapper->mapInBoundary (id);
      for (int i = 0; i < dim; ++i) {
        adjacencyPattern[(n_T+ii)*dim+i].insert ((n_T+ii)*dim+i);  // Id_I
        adjacencyPattern[(n_N+ii)*dim+i].insert ((n_T+ii)*dim+i);  // D_I
      }
    } else if (active.find (id) != active.end()) {
      int ia = aiMapper->mapInActive (id);
      for (int i = 0; i < dim; ++i)
        adjacencyPattern[(n_N+n_I+ia)*dim+i].insert((n_T+n_I+ia)*dim+i); // D_A

      int ii = (n_T+n_I)*dim + ia;
       // T_A
      adjacencyPattern[ii].insert ((n_T+aiMapper->mapInBoundary (id))*dim);
      adjacencyPattern[ii].insert ((n_T+aiMapper->mapInBoundary (id))*dim+1);
        // N_A
      adjacencyPattern[ii+n_A].insert (aiMapper->map (id)*dim);
      adjacencyPattern[ii+n_A].insert (aiMapper->map (id)*dim+1);
    }
  }
    //cout << "Adjacency computed.\n";
  
  for (auto row = B.createbegin(); row != B.createend(); ++row)
    for (const auto& col : adjacencyPattern[row.index()])
        row.insert (col);
  
  B = 0.0;

    // Copy A
  for (auto row = A.begin(); row != A.end(); ++row) {
    auto r = row.index();
    for (auto col = A[r].begin(); col != A[r].end(); ++col) {
      auto c = col.index();
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          if (A[r][c][i][j] != 0.0)
            B[r*dim+i][c*dim+j] = A[r][c][i][j];
    }
  }
    // Copy D
  for (int r = 0; r < n_I+n_A; ++r)
    for (int i = 0; i < dim; ++i)
      B[(r+n_N)*dim+i][(r+n_T)*dim+i] = D[r][r][i][i];

    // Copy Id_I
  for (int r = 0; r < n_I; ++r)
    for (int i = 0; i < dim; ++i)
      B[(r+n_T)*dim+i][(r+n_T)*dim+i] = 1.0;
  
    //cout << "B initialized (1/2).\n";

    // Fill the lines:
    //             |    0      0      0      0   T_A |
    //             |    0      0    N_A      0     0 |
  
  const auto& multBasis = TLM::instance();
  g = 0.0;

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
            int ib = aiMapper->mapInBoundary (id);
            coord_t   nr = is->centerUnitOuterNormal();
            coord_t D_nr = FMatrixHelp::mult (D[ib][ib], nr);
              //cout << "ib= " << ib << ", nr= " << nr << "\n";
            
              // FIXME: we shouldn't compute this here
            n_d[ib] += D_nr;//*(1.0/vnum);
            
              // nor this:
              // Gap at boundary for the computation of the active index set
              // FIXME: yet another stupid recomputation. Should reorder instead.
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (is->type(), quadratureOrder)) {
                // Transform relative (dim-1)-dimensional coord. in local coord.
              const auto& global = igeo.global (x.position ());
              const auto&  local = it->geometry().local (global);
              g[ib] += gap (global) *
              multBasis[subi].evaluateFunction (local) *
              x.weight () *
              igeo.integrationElement (x.position ());
            }
            
            if (active.find (id) != active.end()) {
              int ia = aiMapper->mapInActive (id);
                //cout << "Found active: " << id << " -> " << ia << "\n";
              auto ipos = is->inside()->template subEntity<dim>(subi)->geometry().center();
              
              coord_t tg;
              tg[0] = -nr[1]; tg[1] = nr[0]; for (int r=2; r<dim; ++r) tg[r] = 0;

                // first the tangential stuff for the multipliers.
              int ii = (n_T + n_I)*dim + ia;
                //cout << "ii= " << ii << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][(n_T+ib)*dim+j] += tg[j];//*(1.0/vnum);
              c[ii] = 0.0;
              
                // now the last rows
              ii = (n_T + n_I)*dim + n_A + ia;
              int jj = aiMapper->map (id)*dim;
                //cout << "ii= " << ii << ", jj= " << jj << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][jj+j] += D_nr[j];//*(1.0/vnum);
              c[ii] = gap(ipos);
            }
          }
        }
      }
    }
  }
    //cout << "B initialized (2/2).\n";
  bench().report ("Stepping", " done.");

    //printmatrix(cout, B, "B", "");
    //printvector(cout, c, "c", "");
  writeMatrixToMatlab (B, "/tmp/B");
  writeVectorToFile (c, "/tmp/c");

  bench().report ("Stepping", "Solving", false);

    //HACK: multiply the system by transpose(B) by the left to make it symmetric.
    // and remove zeroes from the diagonal or ILUn will hang
  ScalarMatrix BB;
  BB.setSize (total*dim, total*dim);
  transposeMatMultMat(BB, B, B);
  ScalarVector cc;
  cc.resize (total*dim, false);
  B.mtv (c, cc);
 
  writeMatrixToMatlab (BB, "/tmp/BB");
  writeVectorToFile (cc, "/tmp/cc");

  try {
    InverseOperatorResult stats;
    MatrixAdapter<ScalarMatrix, ScalarVector, ScalarVector> op (BB);
    
      //DummyPreconditioner<ScalarMatrix, ScalarVector, ScalarVector> pre;
      //SeqSSOR<BlockMatrix, ScalarVector, ScalarVector> ssor (BB, 1, 0.95);
      //BiCGSTABSolver<ScalarVector> bcgs (op, ilu, 1e-9, 400, 1);
    SeqILUn<ScalarMatrix, ScalarVector, ScalarVector> ilu (BB, 1, 0.90);
    CGSolver<ScalarVector> cgs (op, ilu, 1e-10, 50000, 1);
    cgs.apply (uu, cc, stats);
    /*
     InverseOperatorResult stats;
     MatrixAdapter<ScalarMatrix, ScalarVector, ScalarVector> op (B);
     SeqSSOR<ScalarMatrix, ScalarVector, ScalarVector> ssor (B, 1, 0.95);
     BiCGSTABSolver<ScalarVector> bcgs (op, ssor, 1e-17, 300, 2);
     bcgs.apply (uu, c, stats);
     */
    /*
     SeqSSOR<ScalarMatrix, ScalarVector, ScalarVector> ssor (B, 1, 0.95);
     RestartedGMResSolver<ScalarVector> rgmres (op, ssor, 1e-8, n_T/2, 200, 0);
     rgmres.apply (uu, c, stats);
     */
    /*
     // This one seemed to work for a while
     SeqSSOR<ScalarMatrix, ScalarVector, ScalarVector> ssor (B, 1, 0.96);
     LoopSolver<ScalarVector> lsol (op, ssor, 1e-8, 200, 1);
     lsol.apply (uu, c, stats);
     */
    
    /*
     SeqSSOR<ScalarMatrix, ScalarVector, ScalarVector> ssor (B, 1, 0.96);
     CGSolver<ScalarVector> cgs (op, ssor, 1e-8, 200, 2);
     cgs.apply (uu, c, stats);
     */
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
    exit(1);
  }

  bench().report ("Stepping", " done.");

    //// FIXME: the following stuff doesn't belong here.
  
  for (int i = 0; i < total; ++i)
    for (int j = 0; j < dim; ++j)
      u[i][j] = uu[i*dim+j];

  for (const auto& x : inactive) {
    int i = aiMapper->mapInBoundary(x);
    int j = aiMapper->map(x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
  for (const auto& x : active) {
    int i = aiMapper->mapInBoundary(x);
    int j = aiMapper->map(x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
    //printvector (cout, n_d, "D * Normal", "");
  printvector (cout, n_u, "Solution Normal", "");
  printvector (cout, n_m, "Multiplier Normal", "");
}

template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::solve ()
{
  PostProcessor<TGV, TET, AIMapper, TSS> post (gv, *aiMapper, a);
  const int maxiter = 20;
  int cnt=0;
  while (true && ++cnt < maxiter) {
    bench().start ("Active set initialization", false);
    determineActive();  // needs g, n_u, n_m computed from last iteration or =0
    bench().stop ("Active set initialization");
    bench().start ("Adjacency computation", false);
    setupMatrices();  // resets sparse matrix info according to new active/inactive sets
    bench().stop ("Adjacency computation");
    bench().start ("Assembly");
    assemble();  // reassembles matrices FIXME: should just reorder.
    bench().stop ("Assembly");
    bench().start ("Stepping");
    step();
    bench().stop ("Stepping");
    bench().start ("Postprocessing", false);
    (void) post.computeError (u);
    post.computeVonMisesSquared ();
    (void) post.writeVTKFile ("/tmp/SignoriniFEM", cnt);
    bench().stop ("Postprocessing");
  }
}
#endif /* defined (SIGNORINI_ACTIVESET_HPP) */
