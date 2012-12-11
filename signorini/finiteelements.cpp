/******************************************************************************
 * finiteelements.cpp                                                         *
 ******************************************************************************/

#include "config.h"
#include <iostream>
#include <vector>
#include <set>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fassign.hh>     // operator <<= for vectors
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/albertagrid.hh>

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

using namespace Dune;
using std::cout;

/*! SignoriniFEPenalty for the Signorini problem using a penalty method.

 Implements an iterative scheme for the approximate computation of the solutions
 to the penalized problem.
 
 TODO:
  - Drop <TFunctor>::isSupported(): specify directly on the mesh where the
    different boundaries are.
  - Use some sort of mask to reduce the number of calculations in the penalty
  - CHECK: "Big penalty parameters => ill-conditioned systems"
  - Make the elements (P1, Q1) parametrizable.
  - Generalize to higher order elements.
 
 Template type names:
    TGV: TGridView
    THT: THookeTensor
    TFF: TForcesFunctor
    TTF: TTractionsFunctor
    TGF: TGapFunctor
 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
class SignoriniFEPenalty
{
public:
  static const int dim = TGV::dimension;
  
  typedef            typename TGV::ctype ctype;
  typedef      FieldVector<ctype, dim> coord_t;
  typedef FieldMatrix<ctype, dim, dim> block_t;
  typedef      BCRSMatrix<block_t> BlockMatrix;
  typedef     BlockVector<coord_t> CoordVector;
  
private:
  const TGV& gv;
  
  const THT& a;   //!< Hooke Tensor
  const TFT& f;   //!< Volume forces
  const TTT& p;   //!< Boundary forces
  const TGT& g;   //!< Normal gap function (scalar)
  
  block_t I;      //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  
  BlockMatrix A;  //!< Stiffness matrix
  BlockMatrix P;  //!< Penalty matrix
  //BlockMatrix M;  //!< Mask matrix for the computation of the positive part of the penalty
  
  CoordVector b;  //!< RHS: volume forces and tractions
  CoordVector r;  //!< RHS: penalty contributions
  //CoordVector m;  //!< RHS: mask vector for the computation of the positive part of the penalty
  
  CoordVector u;  //!< Solution
  
  std::vector<std::set<int> > adjacencyPattern;
  
public:
  SignoriniFEPenalty (const TGV& _gv,  const THT& _a, const TFT& _f, const TTT& _p, const TGT& _g);
  
  void determineAdjacencyPattern ();
  void initialize ();
  
  void assembleMain ();
  void assemblePenalties (double epsilon);
  
  void solve ();
  
  const CoordVector& solution() const { return u; }
  std::vector<ctype> solutionAsVector() const {
    std::vector<ctype> ret;
    for (auto& p : u)
      for (auto& c : p)
        ret.push_back (c);  // HACK, TEST
    return ret;
  }
};

template<class TGV, class THT, class TFT, class TTT, class TGT>
SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::SignoriniFEPenalty
(const TGV& _gv,  const THT& _a, const TFT& _f, const TTT& _p, const TGT& _g)
: gv (_gv), a (_a), f (_f), p (_p), g(_g)
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
}


/*! Store adjacency information in a vector of sets.
 
 For each element traverse all its vertices and set them as adjacent to one
 another.
 
 NOTE: the previous version in dunetests was traversing faces of each element,
 using the intersection iterator, then storing V⨯V in the adjacency list, where
 V = {vertices of the face}. This is clearly a subset of
 
   U = {(x,y): x,y ∈ V={vertices of the element}}
 
 but assemble() inserts values in the stiffnes matrix A for every (x,y) ∈ U,
 which here (Signorini, vector valued functions), but not there (Poisson, scalar
 valued) causes exceptions to be thrown...
 
 WTF??
 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
void SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::determineAdjacencyPattern ()
{
  cout << "Following adjacency intuitions... ";
  const auto& set = gv.indexSet ();
  const auto    N = gv.size (dim);
  
  adjacencyPattern.resize (N);
  
  unsigned int cnt = 0;
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
    
    int vnum = ref.size(dim);
    for (int i = 0; i < vnum; ++i)
      for (int j = 0; j < vnum; ++j, ++cnt)
        adjacencyPattern[set.subIndex (*it, i, dim)].insert (set.subIndex (*it, j, dim));

    /* This worked with 1-dim block matrices in the stiffness matrix and the laplace operator...
       Why??? And what have I overlooked?
     
     Instead we traverse all vertices above.
     
      // traverse all intersections (codim-1 entities) of the current element and
      // store all pairs of vertices in adjacencyPattern
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      int vertexsize = ref.size (is->indexInInside (), 1, dim);
      for (int i = 0; i < vertexsize; ++i) {
        int ii = set.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
        for (int j = 0; j < vertexsize; ++j) {
          int jj = set.subIndex (*it, ref.subEntity (is->indexInInside (), 1, j, dim), dim);
          adjacencyPattern[ii].insert (jj);
        }
      }
    }*/
  }

  cout << " Performed " << cnt << " insertions.\n";
}


/*! One-time initialization.

 Sets the sizes and sparsity patterns for the stiffness and penalty matrices.
 Initializes all matrices and vectors to zero.
 
 FIXME! How do I set the number of nonzeros for the matrices?
 The value N + 2*gv.size (dim-1) causes an exception to be thrown:
   "Specified number of nonzeros ... not sufficient for calculated nonzeros ..."
 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
void SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::initialize ()
{
  determineAdjacencyPattern();
  
  cout << "Randomizing inverse flow generator... ";
  
  const int N = gv.size (dim);
  
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
  b = 0.0;
  u = 0.0;
  
  cout << "ok.\n";
}


/*! Assemble stiffness matrix A and right side b.
 
 BOUNDARY CONDITIONS:
 
 We traverse all faces of all elements with the intersection iterator gv.ibegin()...
 Boundary conditions must be so set, that the Neumann and Signorini boundaries have
 no vertices in common. For every Neumann face we compute the contributions
 to the system and tag all the vertices as visited. Finally we set all untagged
 vertices as being of Dirichlet kind.
 
 We must do the latter in a separate loop through all boundary faces, because of the
 destructive way in which we impose the Dirichlet condition: even if we did allow
 overriding of Dirichlet nodes by the other kinds, we'd have empty rows in the stiffness
 matrix.
 
 FIXME: the method with the std::set (which find()s a key in log time) is probably
 awfully inefficient.
 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
void SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::assembleMain ()
{
  cout << "Gathering gremlins... ";
  std::set<int> boundaryVisited;
  const auto&  iset = gv.indexSet ();
  const auto& basis = Q1ShapeFunctionSet<ctype, dim>::instance ();
  
    //cout << "*** Traversing codim 0 leaves:" << endl;
  
  /*
   We compute the stiffness matrix in an intuitive way: for each entry A_ij one
   must integrate (a function of) the gradients on all the domain, but instead
   of traversing the whole grid for every entry, we traverse once and compute the
   contribution of each of the basis functions on each element to their corresponding
   matrix entry. This in turn uses the fact that we use P1 or Q1 elements with
   exactly as many nodes as vertices has an element (that's why we may use vnum)
   */
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType gt = it->type ();
    const auto& ref = GenericReferenceElements<ctype, dim>::general (gt);
    const auto&  ig = it->geometry();
      //printCorners (ig);
    
    int vnum = ref.size (dim);   // Number of vertices
    
      // get a quadrature rule of order two for the geometry type, then
      // compute transformed gradients, then
      // gain global indices of vertices i and j and update associated matrix entry
    
    for (auto& x : QuadratureRules<ctype, dim>::rule (gt, 2)) {
      block_t jacInvTra = ig.jacobianInverseTransposed (x.position ());
      coord_t grad1, grad2;

      for (int i = 0; i < vnum; ++i) {
        jacInvTra.mv (basis[i].evaluateGradient (x.position ()), grad1);

        for (int j = 0; j < vnum; ++j) {
          jacInvTra.mv (basis[j].evaluateGradient (x.position ()), grad2);
          
          auto ii = iset.subIndex (*it, i, dim);
          auto jj = iset.subIndex (*it, j, dim);
          try {
            A[ii][jj] +=
              a (grad1, grad2) * x.weight () * ig.integrationElement (x.position ());
          } catch (ISTLError& e) {
              // If we are here, the adjacencyPattern does not match.
            cout << "FAILED setting data for A[" << ii << ", " << jj << "] = ";
            cout << "(" << ig.corner(i) << ") x (" << ig.corner (j) << ")\n";
          }
        }
      }
    }

      //// Integrand of the RHS (using a quadrature rule of order two)
    
    for (auto& x : QuadratureRules<ctype, dim>::rule (gt, 2))
      for (int i = 0 ; i < vnum; ++i)
        b[iset.subIndex (*it, i, dim)] += f (it->geometry ().global (x.position ())) *
                                          basis[i].evaluateFunction (x.position ()) *
                                          x.weight () *
                                          it->geometry ().integrationElement (x.position ());

      //// Boundary conditions. See the discussion in the comments to the method
      
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int  vnum = ref.size (is->indexInInside (), 1, dim);  // Number of vertices
        const auto& ig  = is->geometry ();
        const auto& igt = is->type ();
        
        /* Wouldn't it be nice to have lambdas?
         
        for (int i = 0 ; i < vertexnum; ++i) {
          int ii = set.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
          auto integrand = [&] (const coord_t& local) -> coord_t {  // assuming integrate() handles vectors...
            return basis[i].evaluateFunction (it->geometry().local (ig.global (local))) *
                   g (ig.global (local)) * is->unitOuterNormal (local); };

          b[ii] = integrate (integrand, ig);
        }
         */
        
        if (p.isSupported (ig)) {            // Neumann conditions.
          //cout << "Neumann'ing: "; printCorners (ig);
          for (int i = 0 ; i < vnum; ++i) {
            int ii = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
              //cout << "Neumann'ing node: " << ii << "\n";
            boundaryVisited.insert (ii);
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (igt, 2)) {
              b[ii] += p (ig.global (x.position ())) *
                       basis[i].evaluateFunction (it->geometry().local (ig.global (x.position ()))) *
                       x.weight () *
                       ig.integrationElement (x.position ());
            }
          }
        }/*  else if (g.isSupported (ig)) {    // Signorini conditions.
          //cout << "Signorini'ing: "; printCorners (ig);
          for (int i = 0 ; i < vnum; ++i) {
            int ii = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
              //cout << "Signorini'ing node: " << ii << "\n";
            boundaryVisited.insert (ii);
            
            for (int j = 0; j < vnum; ++j) {
              int jj = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, j, dim), dim);
                //cout << "Signorini'ing node: " << jj << "\n";
              boundaryVisited.insert (jj);
            }
          }
        }*/
      }
    }
  }
  
    /* Dirichlet boundary conditions.
     
     For unvisited vertices on the boundary, replace the associated line of A 
     and b with a trivial one.
     
     As stated before, this must be done in a separate step.
     */

  coord_t dirichlet;
  dirichlet <<= zero; // 0.0, -0.07;
  
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType gt = it->type ();
    const auto& ref = GenericReferenceElements<ctype, dim>::general (gt);
    
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const int vertexnum = ref.size (is->indexInInside (), 1, dim);
          //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
        
        for (int i = 0; i < vertexnum; ++i) {
          int ii = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
          if (boundaryVisited.count (ii) > 0)
            continue;
            //cout << "Dirichlet'ing node: " << ii << "\n";
          boundaryVisited.insert(ii);
          A[ii] = 0.0;
          A[ii][ii] = I;
          b[ii] = dirichlet;
        }
      }
    }
  }

  cout << " Done. Total nodes constrained: " << boundaryVisited.size() << ".\n";
}


/*! Assemble the penalty matrix and rhs vector.
 

 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
void SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::assemblePenalties (double eps)
{
  const auto&  iset = gv.indexSet ();
  const auto& basis = Q1ShapeFunctionSet<ctype, dim>::instance ();
  
  eps = 1.0 / eps;
  P   = 0.0;
  r   = 0.0;
  
    //cout << "Penalizing nodes: ";
  for (auto it = gv.template begin<0>(); it != gv.template end<0>(); ++it) {
    GeometryType gt = it->type ();
    const auto& ref = GenericReferenceElements<ctype, dim>::general (gt);
    const int vnum = ref.size (0, 1, dim);
    
      // Iterate through all intersections
    for (auto is = gv.ibegin (*it) ; is != gv.iend (*it) ; ++is) {
      if (is->boundary ()) {
        const auto& ig  = is->geometry ();
                
        if (g.isSupported (ig)) {  // Possible contact zone.
          const auto&  igt = is->type ();
          const auto& rule = QuadratureRules<ctype, dim-1>::rule (igt, 2);
          block_t  penalty (0.0);

          for (int i = 0 ; i < vnum; ++i) {
            int ii = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, i, dim), dim);
            
              // At each quadrature point, check whether for the previous solution
              // u_{t-1}*n-g was positive and add u_t*n-g only in that case
            for (auto& x : rule) {
              coord_t n = is->unitOuterNormal (x.position ());
                // FIXME! Shouldn't I try to interpolate the solution to evaluate it in the quadrature point?
              if (n * u[ii] - g (ig.global (x.position())) > 0) {
                  //cout << " " << ii;
                r[ii] += n *
                         basis[i].evaluateFunction (it->geometry().local (ig.global (x.position ()))) *
                         g (ig.global (x.position ())) *
                         x.weight() *
                         ig.integrationElement (x.position ()) *
                         eps;
                for (int j = 0; j < vnum; ++j) {
                  int jj = iset.subIndex (*it, ref.subEntity (is->indexInInside (), 1, j, dim), dim);
                  for (int k = 0; k < dim; ++k) {
                    for (int l = 0; l < dim; ++l) {
                      
                      penalty[k][l] = n[k] *
                      basis[i].evaluateFunction (it->geometry().local (ig.global (x.position ()))) *
                      n[l] *
                      basis[j].evaluateFunction (it->geometry().local (ig.global (x.position ()))) *
                      x.weight() *
                      ig.integrationElement (x.position ()) *
                      eps;
                    }
                    P[ii][jj] += penalty;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
    //cout << ".\n";
}


/*! Solve (A+P)u = b+r for u.

 Arguments for the BiCGSTABSolver:
 
 op        The operator we solve.
 prec      The preconditioner to apply in each iteration of the loop.
 reduction The relative defect reduction to achieve when applying the operator.
 maxit     The maximum number of iteration steps allowed when applying the operator.
 verbose   The verbosity level. (0,1,2)
 
 */
template<class TGV, class THT, class TFT, class TTT, class TGT>
void SignoriniFEPenalty<TGV, THT, TFT, TTT, TGT>::solve ()
{
  BlockMatrix B = A;//  B += P;      // Add penalties
  CoordVector c = b;//  c += r;      // Add penalties

  /*
  cout << "Writing sparse matrix to /tmp/stiff\n";
  writeMatrixToMatlab (B, "/tmp/stiff");

  cout << "Writing rhs vector to /tmp/rhs\n";
  writeVectorToFile (c, "/tmp/rhs");
  */
  
  InverseOperatorResult stats;                                      // statistics of the solver
  MatrixAdapter<BlockMatrix, CoordVector, CoordVector> op (B);      // make linear operator with A+P
  
  SeqILUn<BlockMatrix, CoordVector, CoordVector> ilu1 (B, 1, 0.96); // initialize preconditioner
  BiCGSTABSolver<CoordVector> bcgs (op, ilu1, 1e-15, 500, 0);      // Bi-conjugate gradient solver
  bcgs.apply (u, c, stats);
  
  /* Using another solver...
  
  SeqSSOR<BlockMatrix, CoordVector, CoordVector> ssor (B, 1, 1.0);  // SSOR preconditioner
  CGSolver<CoordVector> cgs (op, ssor, 1E-10, 5000, 2);             // CG solver
  cgs.apply (u, c, stats);
   
   */
}
