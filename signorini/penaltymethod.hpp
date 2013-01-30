/******************************************************************************
 * penaltymethod.hpp                                                          *
 ******************************************************************************/

#include "config.h"
#include <string>
#include <iostream>
#include <vector>
#include <set>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fassign.hh>     // operator <<= for FieldVectors
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include "shapefunctions.hpp"
#include "traverse.hpp"
#include "utils.hpp"
#include "benchmark.hpp"
#include "postprocessor.hpp"

using namespace Dune;
using std::string;


/*! Solution of the one Signorini problem using a penalty method.

 Implements an iterative scheme for the approximate computation of the solutions
 to the penalized problem.
 
 TODO:
  - Better <TFunctor>::isSupported(): specify directly on the mesh where the
    different boundaries are (and use a Mapper)
  - Use some sort of mask to reduce the number of calculations in the penalty
  - CHECK: "Big penalty parameters => ill-conditioned systems"
  - Generalize to higher order elements.
 
 Template type names:
    TGV: TGridView
    TET: TElasticityTensor
    TFF: TForcesFunctor
    TTF: TTractionsFunctor
    TGF: TGapFunctor
    TSS: TShapeFunctionSet
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
class SignoriniFEPenalty
{
public:
  static const int dim = TGV::dimension;
  
  typedef            typename TGV::ctype ctype;
  typedef      FieldVector<ctype, dim> coord_t;
  typedef FieldMatrix<ctype, dim, dim> block_t;
  typedef      BCRSMatrix<block_t> BlockMatrix;
  typedef     BlockVector<coord_t> CoordVector;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                  MCMGVertexLayout> VertexMapper;
  
private:
  const TGV& gv;  //!< Grid view
  
  const TET& a;   //!< Elasticity tensor
  const TFT& f;   //!< Volume forces
  const TTT& p;   //!< Boundary forces
  const TGT& g;   //!< Normal gap function (scalar)
  
  block_t     I;  //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  BlockMatrix A;  //!< Stiffness matrix
  BlockMatrix P;  //!< Penalty matrix
  CoordVector b;  //!< RHS: volume forces and tractions
  CoordVector r;  //!< RHS: penalty contributions
  CoordVector u;  //!< Solution
  
  int quadratureOrder;
  double          eps; //!< Penalty parameter
public:
  SignoriniFEPenalty (const TGV& _gv,  const TET& _a, const TFT& _f,
                      const TTT& _p, const TGT& _g, double _eps,
                      int _quadratureOrder = 4);

  void initialize ();
  void solve (int maxsteps, double tolerance);

  const CoordVector& solution() const { return u; }
  
protected:
  void assembleMain ();
  void assemblePenalties ();
  void solveSystem ();
};


/******************************************************************************
 * Implementation                                                             *
 ******************************************************************************/

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::SignoriniFEPenalty
(const TGV& _gv,  const TET& _a, const TFT& _f, const TTT& _p, const TGT& _g,
 double _eps, int _quadratureOrder)
: gv (_gv), a (_a), f (_f), p (_p), g(_g), quadratureOrder(_quadratureOrder), eps(_eps)
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
}


/*! One-time initialization: Sets the sizes and sparsity patterns for the 
 stiffness and penalty matrices, then assembles the stiffness matrix.
 
 FIXME: How do I set the number of nonzeros for the matrices?
 The value N + 2*gv.size (dim-1) causes an exception to be thrown:
   "Specified number of nonzeros ... not sufficient for calculated nonzeros ..."
 
 NOTE: for the computation of adjacency information, the previous version in
 dunetests was traversing faces of each element, using the intersection iterator
 then storing V⨯V in the adjacency list, where V = {vertices of the face}.
 For triangular elements this is ok, but for cuadrilateral ones this method
 results clearly in a subset of
 
 U = {(x,y): x,y ∈ V={vertices of the element}}
 
 so here we traverse all vertices of each leaf of codim 0 instead.
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::initialize ()
{
  bench().start ("Initialization", false);
  
  
    //// 1. Compute adjacency information.
  
  bench().report ("Initialization", "Imbuing with adjacency intuitions...", false);
  
  const auto& set = gv.indexSet ();
  const int N = gv.size (dim);
  std::vector<std::set<int> > adjacencyPattern (N);
  
    //For each element we traverse all its vertices and set them as adjacent
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    
    int vnum = ref.size (dim);
    for (int i = 0; i < vnum; ++i) {
      for (int j = 0; j < vnum; ++j) {
        auto& cur = adjacencyPattern[set.subIndex (*it, i, dim)];
        cur.insert (set.subIndex (*it, j, dim));
      }
    }
  }
  
  bench().report ("Initialization", "\tdone.");
  
    //// 2. Initialize matrices and vectors
  
  bench().report ("Initialization", "Randomizing inverse flow generator...", false);

  A.setSize (N, N);
  A.setBuildMode (BlockMatrix::random);
  
  P.setSize (N, N);
  P.setBuildMode (BlockMatrix::random);
  
  b.resize (N, false);
  r.resize (N, false);
  u.resize (N, false);

  for (int i = 0; i < N; ++i) {
    A.setrowsize (i, adjacencyPattern[i].size ());
    P.setrowsize (i, adjacencyPattern[i].size ());
  }
  A.endrowsizes ();
  P.endrowsizes ();
  
  for (int i = 0; i < N; ++i) {
    for (const auto& it : adjacencyPattern[i]) {
      A.addindex (i, it);
      P.addindex (i, it);
    }
  }
  
  A.endindices ();
  P.endindices ();
  
  A = 0.0;
  P = 0.0;
  b = 0.0;
  r = 0.0;
  u = 0.0;
  
  bench().report ("Initialization", "\tdone.");
  
    //// 3. Assemble stiffness matrix

  assembleMain ();
  
  bench().stop ("Initialization");
}


/*! Assemble stiffness matrix A and right side b.
 
 BOUNDARY CONDITIONS:
 
 We traverse all faces of all elements with the intersection iterator gv.ibegin()...
 Boundary conditions must be set so that the Neumann and Signorini boundaries have
 no vertices in common. For every Neumann face we compute the contributions
 to the system.
 
 Last we apply Dirichlet conditions in a separate loop through all boundary faces, 
 because of the destructive way in which we impose the Dirichlet condition: even 
 if we did allow overriding of Dirichlet nodes by the other kinds, we'd have
 empty rows in the stiffness matrix.
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::assembleMain ()
{
  bench().report ("Initialization", "Feeding gremlins...", false);
  
  std::set<int> boundaryVisited;
  const auto& basis = TSS::instance ();
  VertexMapper mapper (gv.grid());
  
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
          
          int ii = mapper.map (*it, i, dim);
          int jj = mapper.map (*it, j, dim);
            //cout << "ii= " << ii << ", jj= " << jj << "\n";
          try {
            A[ii][jj] +=
              a (grad2, grad1) * x.weight () * geo.integrationElement (x.position ());
          } catch (ISTLError& e) {       // The adjacencyPattern does not match.
            cout << "FAILED setting data for A[" << ii << ", " << jj << "] = "
                 << "(" << geo.corner(i) << ") x (" << geo.corner (j) << ")\n";
          }
        }
      }
        //// Integrand of the RHS
      
      for (int i = 0 ; i < vnum; ++i)
        b[mapper.map (*it, i, dim)] +=
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
            int   ii = mapper.map (*it, subi, dim);
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
        } else if (g.isSupported (igeo)) {    // Signorini conditions.
          //cout << "Signorini'ing: "; printCorners (ig);
          for (int i = 0 ; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            int   ii = mapper.map (*it, subi, dim);
              //cout << "Signorini'ing node: " << ii << "\n";
            boundaryVisited.insert (ii);
            
            for (int j = 0; j < ivnum; ++j) {
              int subj = ref.subEntity (is->indexInInside (), 1, j, dim);
              int   jj = mapper.map (*it, subj, dim);
                //cout << "Signorini'ing node: " << jj << "\n";
              boundaryVisited.insert (jj);
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
    //dirichlet <<= zero;
  dirichlet[0] = 0; dirichlet[1] = -0.07;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
          //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
        
        for (int i = 0; i < ivnum; ++i) {
          auto subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            //auto v = it->template subEntity<dim> (rsub)->geometry().center();
          int ii = mapper.map (*it, subi , dim);
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
  bench().report ("Initialization", "\tdone.");
}


/*! Assemble the penalty matrix and rhs vector.
 
 At each quadrature point, check whether for the previous solution
 
    u_{t-1}*n-g 
 
 was positive at the node we are considering and add
 
    u_t*n-g
 
 only in that case
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::assemblePenalties ()
{
  bench().report ("Penalty matrix assembly", "Coercing bystanders...", false);
  
  VertexMapper mapper (gv.grid());
  const auto& basis = TSS::instance ();
  
  eps = 1.0 / eps;
  P   = 0.0;
  r   = 0.0;
  
    // Stuff just for display:
  unsigned int pen = 0;
    //unsigned int cnt = 0; const int div = gv.size (dim) / 10;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType typ = it->type ();
    const auto&  ref = GenericReferenceElements<ctype, dim>::general (typ);
    const int   vnum = ref.size (0, 1, dim);
    
      // Iterate through all intersections
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
        //cout << ((++cnt%div == 0) ? "." : "");
      if (is->boundary ()) {
        const auto& igeo  = is->geometry ();
                
        if (g.isSupported (igeo)) {  // Possible contact zone.
          const auto& ityp = is->type ();
          const auto& rule = QuadratureRules<ctype, dim-1>::rule (ityp, 2*(dim-1));
            //const auto& iref = GenericReferenceElements<ctype, dim-1>::general (ityp);
          
          block_t penalty (0.0);

          for (int i = 0 ; i < vnum; ++i) {
            int   isub = ref.subEntity (is->indexInInside (), 1, i, dim);
            int     ii = mapper.map (*it, isub, dim);
            auto iipos = is->inside()->template subEntity<dim>(isub)->geometry().center();

            for (auto& x : rule) {
              auto global = igeo.global (x.position ());
              auto local  = it->geometry().local (global);
              coord_t   n = is->unitOuterNormal (x.position());
              
              if (n * u[ii] - g (iipos) > 0) {
                ++pen;  //cout << " " << ii;
                r[ii] += n * basis[isub].evaluateFunction (local) * g (global) *
                         x.weight() * igeo.integrationElement (x.position ()) * eps;
                         
                for (int j = 0; j < vnum; ++j) {
                  int jsub = ref.subEntity (is->indexInInside (), 1, j, dim);
                  int   jj = mapper.map (*it, jsub, dim);
                  for (int k = 0; k < dim; ++k) {
                    for (int l = 0; l < dim; ++l) {
                      penalty[k][l] = n[k] * basis[isub].evaluateFunction (local) *
                                      n[l] * basis[jsub].evaluateFunction (local) *
                                      x.weight() * eps *
                                      igeo.integrationElement (x.position ());
                    }
                    P[ii][jj] += penalty;
                  }
                }
              } else {
                  //cout << "\nNot penalizing. n= " << n << ", u[" << ii << "]= "
                  // << u[ii] << ", g= " << g(global) << ", result= " << n * u[ii] - g (global);
              }
            }
          }
        }
      }
    }
  }
  
  bench().report ("Penalty matrix assembly",
                  string (string (" (") + pen).append(" nodes constrained)"));
}


/*! Solve (A+P)u = b+r for u.

 Arguments for the BiCGSTABSolver:
 
 op        The operator we solve.
 prec      The preconditioner to apply in each iteration of the loop.
 reduction The relative defect reduction to achieve when applying the operator.
 maxit     The maximum number of iteration steps allowed when applying the operator.
 verbose   The verbosity level. (0,1,2)
 
 */
template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::solveSystem ()
{
  
  BlockMatrix B = A;  B += P;      // Add penalties
  CoordVector c = b;  c += r;      // Add penalties

  /*
  cout << "Writing sparse matrix to /tmp/stiff\n";
  writeMatrixToMatlab (B, "/tmp/stiff");

  cout << "Writing rhs vector to /tmp/rhs\n";
  writeVectorToFile (c, "/tmp/rhs");
   
  printmatrix (cout, P, "Penalty matrix","");
  */
  
  InverseOperatorResult stats;                                      // statistics of the solver
  MatrixAdapter<BlockMatrix, CoordVector, CoordVector> op (B);      // make linear operator with A+P
  
  SeqILUn<BlockMatrix, CoordVector, CoordVector> ilu1 (B, 1, 0.96); // initialize preconditioner
  BiCGSTABSolver<CoordVector> bcgs (op, ilu1, 1e-15, 500, 0);       // Bi-conjugate gradient solver
  bcgs.apply (u, c, stats);
  
  /* Using another solver...
  
  SeqSSOR<BlockMatrix, CoordVector, CoordVector> ssor (B, 1, 1.0);  // SSOR preconditioner
  CGSolver<CoordVector> cgs (op, ssor, 1E-10, 5000, 2);             // CG solver
  cgs.apply (u, c, stats);
   
   */

}

template<class TGV, class TET, class TFT, class TTT, class TGT, class TSS>
void SignoriniFEPenalty<TGV, TET, TFT, TTT, TGT, TSS>::solve (int maxsteps,
                                                              double tolerance)
{
  PostProcessor<TGV, TET, TSS> post (gv, a);
  
  int     step = 0;
  double error = 1.0;
    // HACK: don't stop in the first iterations
  while (step++ < maxsteps &&
         (error > tolerance || (error <= tolerance && step < 5))) {
    bench().start (string("Iteration ") + step);
    
    bench().start ("Penalty matrix assembly", false);
    assemblePenalties ();
    bench().stop ("Penalty matrix assembly");
    
    bench().start ("System resolution", false);
    solveSystem ();
    bench().stop ("System resolution");
    
    bench().start ("Postprocessing", false);
    error = post.computeError (u);
    post.computeVonMisesSquared ();
    (void) post.writeVTKFile ("/tmp/SignoriniFEM", step);
    bench().stop ("Postprocessing");
    
    bench().stop (string("Iteration ") + step);
  }
}


