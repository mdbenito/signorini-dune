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
#include <dune/istl/superlu.hh>

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
  typedef LeafMultipleCodimMultipleGeomTypeMapper
            <typename TGV::Grid, MCMGVertexLayout> VertexMapper;
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
  const int n_T = gv.size (dim);
  const int n_A = static_cast<int> (active.size());
  const int n_I = static_cast<int> (inactive.size());
  const auto ingap = n_A + n_I;
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
  
  b.resize (n_T + n_I + n_A, false);
  u.resize (n_T + n_I + n_A, false);

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
   // all nodes inactive for initial values == 0
  AIFunctor contact (g, n_u, n_m);
  
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

    //bench().report ("Assembly", "Feeding gremlins...", false);
  
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
            A[ii] = 0.0;
            A[ii][ii] = I;
            b[ii] = dirichlet;
          }
        }
      }
    }
  }
  
    //// Calculate submatrix D and  the gap at boundary for the computation of
    //// the active index set
  
  g = 0.0;
  
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
            int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids.subId (*it, subi, dim);
            int    ib = aiMapper->mapInBoundary (id);
              //auto   v = it->template subEntity<dim> (subi)->geometry().center();
              //cout << "Setting index for D: " << ib << " at " << v << "\n";
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
              D[ib][ib] += I * basis[subi].evaluateFunction (local) *
                           multBasis[subi].evaluateFunction (local) *
                           x.weight () *
                           igeo.integrationElement (x.position ());
              g[ib] += gap (global) *
                       multBasis[subi].evaluateFunction (local) *
                       x.weight () *
                       igeo.integrationElement (x.position ());
            }
          }
        }
      }
    }
  }
 
  bench().report ("Assembly", "\tdone.");
  writeMatrixToMatlab(A, "/tmp/A");
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
  const int total = n_T + n_I + n_A;  // number of dim*dim blocks
  
  assert (n_T == n_N + n_A + n_I);
  assert (A.M() == A.N());
  assert (A.N() == n_T);
  assert (D.N() == D.M());
  assert (D.N() == n_A + n_I);

  ScalarVector uu, c;
  c.resize (total*dim, false);
  uu.resize (total*dim, false);
  uu  = 0.0;
  n_d = 0.0;
  n_m = 0.0;
  n_u = 0.0;

  std::vector<int> n_d_count (n_d.size(), 0);  // count of vertices contributing to the computation of each n_d[i]
  int i;
    // Copy RHS vector
  for (i = 0; i < n_T*dim; ++i)
    for (int j = 0; j < dim; ++j)
      c[i*dim+j] = b[i][j];
  cout << "*i= " << i << LF;
    // Rest of entries:
  for (i = n_T*dim; i < (n_T+n_I)*dim+n_A*(dim-1); ++i)  // Rows corresponding to Id_I and T_A
    c[i] = 0.0;
  cout << "**i= " << i << LF;
  for (i = 0; i < n_A; ++i)                  // Rows corresponding to N_A
    // recall that g[] uses aiMapper.inBoundary() ordering.
    c[i+(n_T+n_I)*dim+n_A*(dim-1)] = g[i+n_I];
  cout << "***i= " << i << LF;
  /*
   We copy matrix A and append some columns and lines. B should be
   
         /                                 \
         | A_NN   A_NI   A_NS      0     0 |
         | A_IN   A_II   A_IA    D_I     0 |
         | A_AN   A_AI   A_AA      0   D_A |
     B = |    0      0      0   Id_I     0 |
         |    0      0      0      0   T_A |
         |    0      0    N_A      0     0 |
         \                                /
   
   Recall that T_A and N_A are *not* block matrices but scalar, and that B is
   a square scalar matrix! In three dimensions this is achieved by inserting
   two rows per node in T_A, one per vector in an orthonormal basis
   of the plane whose normal is the unit outer normal vector.
   */

  ScalarMatrix B;
  B.setBuildMode (ScalarMatrix::row_wise);
  B.setSize (total*dim, total*dim);
  cout << " B is " << B.N() << " x " << B.M() << "\n";
  cout << " A is " << A.N() << " x " << A.M() << "\n";
  
    // Flatten the adjacency pattern of A to scalar entries
  std::vector<std::set<int> > adjacencyPattern (total*dim);
  for (auto row = A.begin(); row != A.end(); ++row) {
    auto ri = row.index();
    for (auto col = (*row).begin(); col != (*row).end(); ++col) {
      auto ci = col.index();
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          if (A[ri][ci][i][j] != 0.0)
            adjacencyPattern[ri*dim+i].insert((int)ci*dim+j);
    }
  }

    // Build the adjacency pattern for the entries in B to the right of and below A
  for (auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {
    IdType id = gids.id (*it);
    if (inactive.find (id) != inactive.end()) {
      int im = aiMapper->mapInBoundary (id);
      for (int i = 0; i < dim; ++i) {
        adjacencyPattern[(n_T+im)*dim+i].insert ((n_T+im)*dim+i);  // Id_I
        adjacencyPattern[(n_N+im)*dim+i].insert ((n_T+im)*dim+i);  // D_I
      }
    } else if (active.find (id) != active.end()) {
      int im = aiMapper->map (id);
      int ia = aiMapper->mapInActive (id);
      for (int i = 0; i < dim; ++i)
        adjacencyPattern[(n_N+n_I+ia)*dim+i].insert((n_T+n_I+ia)*dim+i); // D_A

      int ii = (n_T+n_I)*dim;
       // T_A: dim-1 tangent vectors to determine the tangent hyperplane at each node
      for (int j = 0; j < dim - 1; ++j)
        for (int k = 0; k < dim; ++k)
          adjacencyPattern[ii + (dim-1)*ia + j].insert (ii + ia*dim + k);   // FIXME: offsets ok??? <==========================
      assert (ii + ia*(dim-1) + dim-2 < total*dim);
      assert (ii + ia*dim + dim-1 < total*dim);
        // N_A
      for (int k = 0; k < dim; ++k)
        adjacencyPattern[ii+(dim-1)*n_A+ia].insert (im*dim + k);
      assert (ii + (dim-1)*n_A + ia < total*dim);
      assert ((im+1)*dim < total*dim);
    }
  }
  cout << "Adjacency computed.\n";
  
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
  
  cout << "B initialized (1/2).\n";

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
            int    ib = aiMapper->mapInBoundary (id);
            coord_t   nr = is->centerUnitOuterNormal();
            coord_t D_nr = FMatrixHelp::mult (D[ib][ib], nr);
            
            n_d[ib] += D_nr; // FIXME: we shouldn't compute this here
            n_d_count[ib] += 1;
            
              //cout << "ib= " << ib << ", nr= " << nr << "\n";
            
            if (active.find (id) != active.end()) {
              int ia = aiMapper->mapInActive (id);
//              cout << "Found active: " << id << " -> " << ia << "\n";
              
                // first the tangential stuff for the multipliers.
              std::vector<coord_t> tg = basisOfPlaneNormalTo (nr);  // dim-1 items
//              cout << "nr= " << nr << LF;
//              for (int i=0; i<dim-1; ++i)
//                cout << "\t\t tg[" << i << "]= " << tg[i]
//                     << ", norm= " << tg[i].two_norm() << LF;

              int ii = (n_T + n_I)*dim + ia*(dim-1);
                //cout << "ii= " << ii << "\n";
              for (int j = 0; j < dim - 1; ++j)
                for (int k = 0; k < dim; ++k)
                  B[ii+j][(n_T+ib)*dim+k] += tg[j][k]; //FIXME: why a sum?!
              
                // now the last rows
              ii = (n_T + n_I)*dim + n_A*(dim-1) + ia;
              int jj = aiMapper->map (id)*dim;
                //cout << "ii= " << ii << ", jj= " << jj << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][jj+j] += D_nr[j]; //FIXME: why a sum?!
            }
          }
        }
      }
    }
  }

  cout << "Scaling n_d.\n";
  
    // Fix the computation of n_d averaging through the number of vertices
    // contributing to each node.
  for (int i = 0; i < n_d.size(); ++i)
    if (n_d_count[i] > 1)
      n_d[i] /= static_cast<double> (n_d_count[i]);

  cout << "B initialized (2/2).\n";
  bench().report ("Stepping", " done.");

  writeMatrixToMatlab (B, "/tmp/B");
  writeVectorToFile (c, "/tmp/c");

  bench().report ("Stepping", "Solving", false);

  try {
    dinfo.attach (cout);
    InverseOperatorResult stats;
    SuperLU<ScalarMatrix> slu (B, true);
    slu.apply (uu, c, stats);
    dinfo.flush();
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
    //printvector (cout, n_u, "Solution Normal", "");
    //printvector (cout, n_m, "Multiplier Normal", "");
}

template<class TGV, class TET, class TFT, class TTT, class TGF, class TSS, class TLM>
void SignoriniIASet<TGV, TET, TFT, TTT, TGF, TSS, TLM>::solve ()
{
  PostProcessor<TGV, TET, AIMapper, TSS> post (gv, *aiMapper, a);
  const int maxiter = 10;
  int cnt=0;
  while (true && ++cnt <= maxiter) {
    bench().start ("Active set initialization", false);
    determineActive();  // needs g, n_u, n_m computed from last iteration or =0
    bench().stop ("Active set initialization");
    bench().start ("Adjacency computation", false);
    setupMatrices();  // resets sparse matrix info according to new active/inactive sets
    bench().stop ("Adjacency computation");
    bench().start ("Assembly", false);
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
