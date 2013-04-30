/******************************************************************************
 * twobodies.hpp                                                              *
 * Solve the two body contact problem with an active-inactive set strategy.   *
 *                                                                            *
 * This implements the algorithm described in [HW05].                         *
 ******************************************************************************/

#ifndef TWOBODIES_ACTIVESET_HPP
#define TWOBODIES_ACTIVESET_HPP

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
#include "twobodymapper.hpp"

using namespace Dune;
using std::cout;

#define MASTER 0
#define SLAVE  1
#define DOBOTH(x) for (int x=0; x < 2; ++x)


/*!
 
 Template type names:
 TGV: TGridView
 TET: TElasticityTensor
 TFF: TForcesFunctor
 TDF: TDirichletFunctor
 TTF: TTractionsFunctor
 TGF: TGapFunctor
 TSS: TShapeSet
 TLM: TLagrangeMultipliersShapeSet
 */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
class TwoBodiesIASet
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
  typedef BCRSMatrix<block_t>                  BlockMatrix;
  typedef BlockVector<scalar_t>               ScalarVector;
  typedef BlockVector<coord_t>                 CoordVector;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector>  AIFunctor;
  typedef TwoBodyMapper<dim, TGV>                     TwoMapper;
  typedef LeafMultipleCodimMultipleGeomTypeMapper
          <typename TGV::Grid, MCMGVertexLayout>   VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGF>   GapVertexMapper;
  
private:
  TwoRefs<TGV> gv;           //!< Grid view
  TwoRefs<GlobalIdSet> gids; //!< Element ids
  
  TwoRefs<TET> a;    //!< Elasticity tensor
  TwoRefs<TFT> f;    //!< Volume forces
  TwoRefs<TDF> dir;  //!< Dirichlet conditions
  TwoRefs<TTF> p;    //!< Boundary forces
  TwoRefs<TGF> gap;  //!< Normal gap function
  
  block_t      I;    //!< Identity matrix block (Dune::DiagonalMatrix not working?)
  BlockMatrix  A;    //!< Stiffness matrix
  BlockMatrix  D;    //!< See [HW05] (Should be Dune::BDMatrix, but cannot build it)
  BlockMatrix  MM;   //!< Coupling matrix of the mortar method
  BlockMatrix  Q;    //!< Coordinate transformation matrix
  CoordVector  b;    //!< RHS: volume forces and tractions
  CoordVector  u;    //!< Solution
  CoordVector  n_d;  //!< D * normal at the gap nodes
  ScalarVector g;    //!< Gap functor evaluated at the gap nodes
  ScalarVector n_u;  //!< Normal component of the solution at gap nodes
  ScalarVector n_m;  //!< Normal component of the lagrange multiplier at gap nodes
  
  int    quadratureOrder;
  IdSet      inactive[2];
  IdSet        active[2];
  IdSet         other[2];  //!< All nodes (boundary or not) which are not part of the gap
  TwoMapper*   twoMapper;

public:
  TwoBodiesIASet (const TGV& _mgv, const TGV& _sgv, const TET& _ma, const TET& _sa,
                  const TFT& _mf, const TFT& _sf, const TDF& _md, const TDF& _sd,
                  const TTF& _mp, const TTF& _sp, const TGF& _mgap, const TGF& _sgap,
                  int _quadratureOrder = 4);
  
  void setupMatrices ();
  void assemble ();
  void determineActive ();
  void step (int cnt);
  void solve ();
  
  const CoordVector& solution() const { return u; }
};

/******************************************************************************
 * Implementation                                                             *
 ******************************************************************************/

template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::TwoBodiesIASet (const TGV& _mgv,
                                                                        const TGV& _sgv,
                                                                        const TET& _ma,
                                                                        const TET& _sa,
                                                                        const TFT& _mf,
                                                                        const TFT& _sf,
                                                                        const TDF& _md,
                                                                        const TDF& _sd,
                                                                        const TTF& _mp,
                                                                        const TTF& _sp,
                                                                        const TGF& _mgap,
                                                                        const TGF& _sgap,
                                                                        int _quadratureOrder)
:
  gv (TwoRefs<TGV> (_mgv, _sgv)),
  gids (TwoRefs<GlobalIdSet> (_mgv.grid().globalIdSet(), _sgv.grid().globalIdSet())),
  a (TwoRefs<TET> (_ma, _sa)),
  f (TwoRefs<TFT> (_mf, _sf)),
  dir (TwoRefs<TDF> (_md, _sd)),
  p (TwoRefs<TTF> (_mp, _sp)),
  gap (TwoRefs<TGF> (_mgap, _sgap)),
  quadratureOrder (_quadratureOrder),
  inactive(), active(), other()
{
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
  
    //// Initialize the mapper setting all nodes in the gap to inactive:
    /* 
     Careful: vertices may be shared by faces in and out of the gap.
     We mustn't count them twice, and the gap has "priority", meaning that a
     vertex of a face in the gap must always be added to the inactive set,
     hence the use of erase().
     */
  DOBOTH (body) {
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
        const int ivnum = ref.size (is->indexInInside (), 1, dim);
        if (is->boundary() && gap[body].isSupported (*is)) {
          for (int i = 0; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids[body].subId (*it, subi, dim);
            if (other[body].find (id) != other[body].end())
              other[body].erase (id);
            inactive[body] << id;
          }
        } else {
          for (int i = 0; i < ivnum; ++i) {
            int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids[body].subId (*it, subi, dim);
            if (inactive[body].find (id) == inactive[body].end())
              other[body] << id;
          }
        }
      }
    }
  }
  
  DOBOTH (body) {
//    cout << "inactive[" << body << "] = " << inactive[body].size() << LF;
//    cout << "active[" << body << "] = " << active[body].size() << LF;
//    cout << "other[" << body << "] = " << other[body].size() << LF;
    assert (inactive[body].size() + active[body].size() + other[body].size() == gv[body].size (dim));
  }
  
  twoMapper = new TwoMapper (gv, active, inactive, other);

    //// Other initializations
  g.resize (inactive[SLAVE].size());
  n_d.resize (inactive[SLAVE].size());
  n_u.resize (inactive[SLAVE].size());
  n_m.resize (inactive[SLAVE].size());
  n_d = 0.0;
  n_u = 0.0;
  n_m = 0.0;
  g   = 0.0;
}


/*! Initializes the stiffness matrix.
  */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::setupMatrices ()
{
  const auto& basis = TSS::instance();
  
  const int n_Tm = gv[MASTER].size (dim);
  const int n_Ts = gv[SLAVE].size (dim);
  const int  n_T = n_Tm + n_Ts;
  const int n_Am = static_cast<int> (active[MASTER].size());
  const int n_As = static_cast<int> (active[SLAVE].size());
  const int n_Im = static_cast<int> (inactive[MASTER].size());
  const int n_Is = static_cast<int> (inactive[SLAVE].size());
  const auto ingap_m = n_Am + n_Im;
  const auto ingap_s = n_As + n_Is;
//  cout << "total= " << n_T << ", ingap_s= " << ingap_s << "\n";
  
  bench().start ("Adjacency for stiffness matrix", false);
  std::vector<std::set<int> > adjacencyPattern (n_T);
  DOBOTH (body) {
      //For each element we traverse all its vertices and set them as adjacent
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
      const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type ());
      int vnum = ref.size (dim);
      for (int i = 0; i < vnum; ++i) {
        int ii = twoMapper->map (body, *it, i, dim);
        for (int j = 0; j < vnum; ++j) {
          int jj = twoMapper->map (body, *it, j, dim);
          adjacencyPattern[ii].insert (jj);
        }
      }
    }
  }
  bench().stop ("Adjacency for stiffness matrix");

    /* Adjacency pattern for M: the mindless approach.

     Traverse the boundary of the slave, and at each boundary element in the gap
     traverse the boundary of the master looking for basis functions whose support
     intersects.
     
     ****
     
     FIXME: This is just a hack! We should actually project from the master onto
     the slave along the normal and... [complete this!]
     
     FIXME: we check whether the support of a basis function in the master side
     contains the coordinates of a node in the slave side, but we should rather
     check at the *quadrature points* in the slave side to avoid marking pairs
     slave node / master node which actually won't count
     
     FIXME: what happens when the gap is bigger than the mesh size? In the
     following situation, everything is fine: (slave to the left, master to the
     right)
     
             ^  |\
            /|\ | \ φ
        Ψ  / | \|  \
          /  |  \   \
         /   |  |\   \
        /    |  | \   \
     --*-----*  *------*------*---
         h_s             h_m
             |__|
              |
              gap
     
     But if the gap were larger than h_s, then the support of Ψ wouldn't
     intersect the boundary of the master and the matrix M would be empty.
     
     ****
     
     NOTE: for regular SGrids we could simply traverse all nodes of the slave
     then all of the master and if h is the grid width parameter, then
     if |x_0-x_1|+|y_0-y_1| <= h, (the supports of the basis functions intersect)
     add the entry. This should be the same, right?
     
     NOTE2: I should be able to exploit this computation later when actually
     building M in order not to have to traverse all the grid (squared) again, 
     but I just don't feel like it...
     */
  bench().start ("Adjacency for coupling matrix", false);
  std::vector<std::set<int> > couplingPattern (ingap_s);
  for (auto it_s = gv[SLAVE].template begin<0>(); it_s != gv[SLAVE].template end<0>(); ++it_s) {
    for (auto is_s = gv[SLAVE].ibegin (*it_s) ; is_s != gv[SLAVE].iend (*it_s) ; ++is_s) {
      if (is_s->boundary () && gap[SLAVE].isSupported (*is_s)) {
        const auto& ref_s = GenericReferenceElements<ctype, dim>::general (it_s->type());
//        cout << " ********\n";
        for (auto it_m = gv[MASTER].template begin<0>(); it_m != gv[MASTER].template end<0>(); ++it_m) {
          for (auto is_m = gv[MASTER].ibegin (*it_m) ; is_m != gv[MASTER].iend (*it_m) ; ++is_m) {
            if (is_m->boundary () && gap[MASTER].isSupported (*is_m)) {
//              cout << " ------\n";
              const auto& ref_m = GenericReferenceElements<ctype, dim>::general (it_m->type());
              
              const int ivnum_s = ref_s.size (is_s->indexInInside (), 1, dim);
              const int ivnum_m = ref_m.size (is_m->indexInInside (), 1, dim);
              
              for (int i_s = 0 ; i_s < ivnum_s; ++i_s) {
                int subi_s = ref_s.subEntity (is_s->indexInInside (), 1, i_s, dim);
                auto   v_s = it_s->template subEntity<dim> (subi_s)->geometry().center();
                
                for (int i_m = 0 ; i_m < ivnum_m; ++i_m) {
                  int subi_m = ref_m.subEntity (is_m->indexInInside (), 1, i_m, dim);
                  auto   v_m = it_m->template subEntity<dim> (subi_m)->geometry().center();
                  const auto& local_m = it_m->geometry().local (v_s);
                  
                  if (basis[subi_m].isSupported (local_m)) {
//                    cout << "*** " << v_s << " *** Supported by *** " << subi_m
//                         << " *** at " << v_m << "\n";
                    const int ii_s = twoMapper->mapInBoundary (SLAVE, *it_s, subi_s, dim);
                    const int ii_m = twoMapper->mapInBoundary (MASTER, *it_m, subi_m, dim);
                    couplingPattern[ii_s].insert (ii_m);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  bench().stop ("Adjacency for coupling matrix");
    //// Initialize default values
  
  A.setSize (n_T, n_T);
  A.setBuildMode (BlockMatrix::random);
  D.setSize (ingap_s, ingap_s);
  D.setBuildMode (BlockMatrix::random);
  MM.setSize (ingap_s, ingap_m);
  MM.setBuildMode (BlockMatrix::random);
  
  b.resize (n_T, false);
  u.resize (n_T + n_Is + n_As, false); // we use the extra room (+n_Is+n_As) in step()
  
  for (int i = 0; i < n_T; ++i)
    A.setrowsize (i, adjacencyPattern[i].size ());
  A.endrowsizes ();
  
  for (int i = 0; i < n_T; ++i)
    for (const auto& it : adjacencyPattern[i])
      A.addindex (i, it);
  A.endindices ();
  
  for (int i = 0; i < ingap_s; ++i)
    MM.setrowsize (i, couplingPattern[i].size ());
  MM.endrowsizes ();
  
  for (int i = 0; i < ingap_s; ++i)
    for (const auto& it : couplingPattern[i])
      MM.addindex (i, it);
  MM.endindices ();
  
  for (int i = 0; i < ingap_s; ++i)
    D.setrowsize (i, 1);
  D.endrowsizes();
  for (int i = 0; i < ingap_s; ++i)
    D.addindex (i, i);
  D.endindices();

  A = 0.0;
  D = 0.0;
  MM = 0.0;
  b = 0.0;
  u = 0.0;
}


/*!
 */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::determineActive ()
{
   // all nodes inactive for default values == 0
  AIFunctor contact (g, n_u, n_m);
  
  active[SLAVE].clear();
  inactive[SLAVE].clear();
  /*
  for (auto it_s = gv[SLAVE].template begin<dim>(); it_s != gv[SLAVE].template end<dim>(); ++it_s) {
    auto id_s = gids[SLAVE].id (*it_s);
    if (gap[SLAVE].isSupported (it_s->geometry())) {
      if (contact.isSupported (twoMapper->mapInBoundary (SLAVE, id_s)))
        active[SLAVE] << id_s;
      else
        inactive[SLAVE] << id_s;
    }
  }
*/
  for (auto it = gv[SLAVE].template begin<0>(); it != gv[SLAVE].template end<0>(); ++it) {
    for (auto is = gv[SLAVE].ibegin (*it) ; is != gv[SLAVE].iend (*it) ; ++is) {
      const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
      const int ivnum = ref.size (is->indexInInside (), 1, dim);
      if (is->boundary() && gap[SLAVE].isSupported (*is)) {
        for (int i = 0; i < ivnum; ++i) {
          int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          IdType id = gids[SLAVE].subId (*it, subi, dim);
          if (contact.isSupported (twoMapper->mapInBoundary (SLAVE, id)))
            active[SLAVE] << id;
          else
            inactive[SLAVE] << id;
        }
      }
    }
  }

  
  twoMapper->update (active, inactive, other);
  
  cout << "\nInactive: " << inactive[SLAVE].size() << ": ";
    //for (auto x : inactive[SLAVE]) cout << twoMapper->map (SLAVE, x) << " ";
  cout << "\nActive: " << active[SLAVE].size() << ": ";
    //for (auto x : active[SLAVE])   cout << twoMapper->map (SLAVE, x) << " ";
  cout << "\n";
}


/*! Assemble the system matrix.
 */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::assemble ()
{
  const auto& multBasis = TLM::instance();
  const auto&     basis = TSS::instance();
  
    //// Stiffness matrix:
  
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
  DOBOTH (body) {
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
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
            
            int ii = twoMapper->map (body, *it, i, dim);
            int jj = twoMapper->map (body, *it, j, dim);
              //cout << "ii= " << ii << ", jj= " << jj << "\n";
            try {
              A[ii][jj] += a[body] (grad2, grad1) * x.weight () * geo.integrationElement (x.position ());
            } catch (ISTLError& e) {       // The adjacencyPattern does not match.
              cout << "FAILED setting data for A[" << ii << ", " << jj << "] = "
              << "(" << geo.corner (i) << ") x (" << geo.corner (j) << ")\n";
            }
          }
        }
          //// Integrand of the RHS
        
        for (int i = 0 ; i < vnum; ++i)
          b[twoMapper->map (body, *it, i, dim)] +=
                    f[body] (it->geometry ().global (x.position ())) *
                    basis[i].evaluateFunction (x.position ()) *
                    x.weight () *
                    it->geometry ().integrationElement (x.position ());
      }
      
        //// Neumann Boundary conditions.
      
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        if (is->boundary () && p[body].isSupported (*is)) {
          const int  ivnum = ref.size (is->indexInInside (), 1, dim);
          const auto& igeo = is->geometry ();
          const auto& ityp = is->type ();
            //cout << "Neumann'ing: "; printCorners (igeo);
          for (int i = 0 ; i < ivnum; ++i) {
            int subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            int   ii = twoMapper->map (body, *it, subi, dim);
              //auto   v = it->template subEntity<dim> (subi)->geometry().center();
              //cout << "Neumann'ing node: " << ii << " at " << v << "\n";
            
            for (auto& x : QuadratureRules<ctype, dim-1>::rule (ityp, quadratureOrder)) {
                // Transform relative (dim-1)-dimensional coord. in local coord.
              const auto& global = igeo.global (x.position ());
              const auto& local  = it->geometry().local (global);
              b[ii] += p[body] (global) *
              basis[subi].evaluateFunction (local) *
              x.weight () *
              igeo.integrationElement (x.position ());
            }
          }
        }
      }
    }
  
    /* Dirichlet boundary conditions.
     For unvisited vertices on the boundary, replace the associated line of A
     and b with a trivial one.
     */
    
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
      const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
      
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        if (is->boundary () && dir[body].isSupported (*is)) {
          const int ivnum = ref.size (is->indexInInside (), 1, dim);
            //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
          
          for (int i = 0; i < ivnum; ++i) {
            auto subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            auto global = it->geometry().global (ref.position (subi, dim));
            int ii = twoMapper->map (body, *it, subi , dim);
              //cout << "Dirichlet'ing node: " << ii << " at " << v << "\n";
            A[ii] = 0.0;
            A[ii][ii] = I;
            b[ii] = dir[body] (global);            
          }
        }
      }
    }
  }  // end of DOBOTH
  
  g = 0.0;
  
    //// Compute submatrix D and  the gap at boundary for the computation of
    //// the active index set (slave nodes!)
  
    // Recall that the integral is over the gap boundary, so we need not integrate
    // the basis functions of nodes outside it.
  
  std::set<int> check;
  for (auto it = gv[SLAVE].template begin<0>(); it != gv[SLAVE].template end<0>(); ++it) {
    for (auto is = gv[SLAVE].ibegin (*it) ; is != gv[SLAVE].iend (*it) ; ++is) {
      if (is->boundary () && gap[SLAVE].isSupported (*is)) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        
        for (int i = 0 ; i < vnum; ++i) {
          int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          IdType id = gids[SLAVE].subId (*it, subi, dim);
          int    ib = twoMapper->mapInBoundary (SLAVE, id);
          check << ib;
            //auto   v = it->template subEntity<dim> (subi)->geometry().center();
            //cout << "Setting index for D: " << kk << " at " << v << "\n";
          for (auto& x : QuadratureRules<ctype, dim-1>::rule (is->type(), quadratureOrder)) {
              // Transform relative (dim-1)-dimensional coord. in local coord.
            const auto& global = igeo.global (x.position ());
            const auto&  local = it->geometry().local (global);
            
              //               cout << "        quadrature point= " << global << " (" << local
              //               << ")\n" << "           basis eval[" << subi
              //               << "]= " << basis[subi].evaluateFunction (local)
              //               << "\n" << "       multbasis eval[" << subi
              //               << "]= " << multBasis[subi].evaluateFunction (local) << "\n";
            D[ib][ib] += I * basis[subi].evaluateFunction (local) *
                         multBasis[subi].evaluateFunction (local) *
                         x.weight () *
                         igeo.integrationElement (x.position ());
            
            g[ib] += gap[SLAVE] (global) *
                     multBasis[subi].evaluateFunction (local) *
                     x.weight () *
                     igeo.integrationElement (x.position ());
          }
        }
      }
    }
  }
//  cout << "check= " << check.size() << ", D.N()= " << D.N() << LF;
  assert (check.size() >= D.N()); // we should've gone through all vertices in the slave gap
  
    //// Compute mortar coupling matrix MM
  
  /**********   FIXME FIXME FIXME FIXME FIXME *********

   because we traverse intersections, nodes in the master contribute
   four times (one if the node is at the boundary of the gap) to each node in
   the slave. Using the hack with std::set<int> visited to fix this seems not to
   work very well.
   
   */
  for (auto it_s = gv[SLAVE].template begin<0>(); it_s != gv[SLAVE].template end<0>(); ++it_s) {
    for (auto is_s = gv[SLAVE].ibegin (*it_s) ; is_s != gv[SLAVE].iend (*it_s) ; ++is_s) {
      if (is_s->boundary () && gap[SLAVE].isSupported (*is_s)) {
        const auto& igeo_s = is_s->geometry ();
        const auto& ref_s = GenericReferenceElements<ctype, dim>::general (it_s->type());
        
        for (auto it_m = gv[MASTER].template begin<0>(); it_m != gv[MASTER].template end<0>(); ++it_m) {
          for (auto is_m = gv[MASTER].ibegin (*it_m) ; is_m != gv[MASTER].iend (*it_m) ; ++is_m) {
            if (is_m->boundary () && gap[MASTER].isSupported (*is_m)) {
              const auto&  ref_m = GenericReferenceElements<ctype, dim>::general (it_m->type());
              
              const int ivnum_s = ref_s.size (is_s->indexInInside (), 1, dim);
              const int ivnum_m = ref_m.size (is_m->indexInInside (), 1, dim);
              
              for (int i_s = 0 ; i_s < ivnum_s; ++i_s) {
                int subi_s = ref_s.subEntity (is_s->indexInInside (), 1, i_s, dim);
                auto   v_s = it_s->template subEntity<dim> (subi_s)->geometry().center();
                const auto& local_v_m = it_m->geometry().local (v_s);
                const int ii_s = twoMapper->mapInBoundary (SLAVE, *it_s, subi_s, dim);
                for (int i_m = 0 ; i_m < ivnum_m; ++i_m) {
                  int subi_m = ref_m.subEntity (is_m->indexInInside (), 1, i_m, dim);
                  const int ii_m = twoMapper->mapInBoundary (MASTER, *it_m, subi_m, dim);
                  
                  for (auto& x : QuadratureRules<ctype, dim-1>::rule (is_s->type(), quadratureOrder)) {
                    const auto&   global = igeo_s.global (x.position ());
                    const auto&  local_s = it_s->geometry().local (global);
                    const auto&  local_m = it_m->geometry().local (global);
                    
                      //                          cout << " global= " << global
                      //                          << "    local_m= " << local_m
                      //                          << "         local_s= " << local_s << LF
                      //                          << "    basis[" << subi_m
                      //                          << "]= " << basis[subi_m].evaluateFunction (local_m)
                      //                          << "    multbasis[" << subi_s
                      //                          << "]= " << multBasis[subi_s].evaluateFunction (local_s) << LF;
                    if (basis[subi_m].isSupported (local_v_m))
                      MM[ii_s][ii_m] += I * basis[subi_m].evaluateFunction (local_m) *
                                        multBasis[subi_s].evaluateFunction (local_s) *
                                        x.weight () *
                                        igeo_s.integrationElement (x.position ());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
    //printmatrix (cout, D, "D", "");
}


/*! HACK of the HACKS...
 
 This is terrible: we assembled a block matrix and now we flatten it. We should've
 started with a scalar matrix!!! Also, instead of using Q we should carry out the
 multiplication by blocks. Currently, with a mesh of ~20000 nodes each of the
 matMultMats involving A and Q require ~30 sec!!!
 
 At step k=0, the system matrix (without boundary conditions) is
 
     /                           \
     | A_NN   A_NM   A_NS      0 |
     | A_MN   A_MM   A_MS  -MM^T |
     | A_SN   A_SM   A_SS      D |
     \                           /
 
 where N stands for those vertices not in the possible contact zones of both
 bodies; M for the vertices on the gap zone of the master body and S on the
 slave; and MM is the mortar mass matrix coupling the basis functions on both
 bodies at the interface.
 
 Instead of using this matrix, we transform the stiffness matrix A using Q

     /                 \
     | Id    0     0   |
 Q = |  0   Id   ~MM^T |
     |  0    0    Id   |
     \                 /
 
 with ~MM = D^-1*MM. We set Â=Q*A*Q^T and append some columns and lines.
 The final system matrix B is
 
     /                                        \
     | Â_NN   Â_NM   Â_NI   Â_NS     0      0 |
     | Â_MN   Â_MM   Â_MI   Â_MA     0      0 |
     | Â_IN   Â_IM   Â_II   Â_IA    D_I     0 |
     | Â_AN   Â_AM   Â_AI   Â_AA     0    D_A |
 B = |    0      0      0      0   Id_I     0 |
     |    0      0      0      0      0   T_A |
     |    0      0      0    N_A      0     0 |
     \                                        /
 
 Note that this is a SQUARE scalar matrix, because T_A has n_As*(dim-1) rows and
 n_As*dim columns, and N_A has n_As rows and n_As*dim columns.
 
 */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::step (int cnt)
{
  
  bench().report ("Stepping", "Gluing");
  const int   n_T = static_cast<int> (gv[MASTER].size (dim) + gv[SLAVE].size (dim));
  const int   n_N = static_cast<int> (other[MASTER].size() + other[SLAVE].size());
  const int   n_M = static_cast<int> (inactive[MASTER].size() + active[MASTER].size());
  const int   n_S = static_cast<int> (inactive[SLAVE].size() + active[SLAVE].size());
  const int   n_I = static_cast<int> (inactive[SLAVE].size());
  const int   n_A = static_cast<int> (active[SLAVE].size());
  const int total = n_T + n_I + n_A;  // matrix size
  
  assert (n_T == n_N + n_M + n_S);
  assert (A.M() == A.N());
  assert (A.N() == n_T);
  assert (D.N() == D.M());
  assert (D.N() == n_S);
  assert (n_S > 1);  // The slave boundary must have at least two [?] nodes!
  assert (n_M > 0);

  cout << "n_T= " << n_T << ", n_N= " << n_N << ", n_M= " << n_M
       << ", n_S= " << n_S << ", n_I= " << n_I << ", n_A= " << n_A << LF;
  
    // M = D^-1*MM
  BlockMatrix M;
  BlockMatrix DD (D);
  for (int i=0; i < DD.N(); ++i)
    DD[i][i].invert();

  matMultMat (M, DD, MM);
  
  BlockMatrix MT;
  transpose (MT, M);
  cout << "M was transposed.\n";
  BlockMatrix Q;
  Q.setSize (n_T, n_T);
  Q.setBuildMode (BlockMatrix::row_wise);
  
  for (auto row = Q.createbegin(); row != Q.createend(); ++row) {
    auto r = row.index();
    row.insert (r);
    if (n_N <= r && r < n_N+n_M)
      for (auto col = MT[r-n_N].begin(); col != MT[r-n_N].end(); ++col)
        row.insert (col.index()+n_N+n_M);
  }

  Q = 0.0;  // necessary?

  for (int r = 0; r < n_T; ++r) {
    Q[r][r] = I;
    if (n_N <= r && r < n_N+n_M)
      for (auto col = MT[r-n_N].begin(); col != MT[r-n_N].end(); ++col)
        Q[r][col.index()+n_N+n_M] = MT[r-n_N][col.index()];
  }

  writeMatrixToMatlab (A, "/tmp/A");
  writeMatrixToMatlab (D, "/tmp/D");
  writeMatrixToMatlab (DD, "/tmp/DD");
  writeMatrixToMatlab (Q, "/tmp/Q");
  writeMatrixToMatlab (M, "/tmp/M");
  writeMatrixToMatlab (MM, "/tmp/MM");
  writeMatrixToMatlab (MT, "/tmp/MT");

    // The brute force approach:
  BlockMatrix AA;
  matMultTransposeMat (AA, A, Q);
  matMultMat (A, Q, AA);  // A is now the new Â.
  
  writeMatrixToMatlab (A, "/tmp/AA");
  cout << "We have Â\n";
  
  ScalarVector uu, c;
  c.resize (total*dim, false);
  uu.resize (total*dim, false);
  c   = 0.0;
  uu  = 0.0;
  n_d = 0.0;
  n_m = 0.0;
  n_u = 0.0;
  
  CoordVector bb (b);
  Q.mv (bb, b);  // b is now the new ^b=Q*b
  
    // Copy RHS vector
  for (int i = 0; i < n_T; ++i)
    for (int j = 0; j < dim; ++j)
      c[i*dim+j] = b[i][j];
  
    // Rest of entries.
    // No need for this loop since we already set c = 0.0.
//  for (int i = n_T*dim; i < (n_T+n_I)*dim+n_A*(dim-1); ++i)
//    c[i] = 0.0;
  
    // recall that g[] uses twoMapper.inBoundary (SLAVE,...) ordering.
  for (int i = 0; i < n_A; ++i) {
//    cout << " c[" << (n_T+n_I)*dim + n_A + i << "]= g[" << i+n_I << "]= " << g[i+n_I] << LF;
    c[(n_T+n_I)*dim + n_A*(dim-1) + i] = g[i+n_I];
  }
  ScalarMatrix B;
  B.setBuildMode (ScalarMatrix::row_wise);
  B.setSize (total*dim, total*dim);
  cout << " B is " << B.N() << " x " << B.M() << "\n";
  
    // Flatten the adjacency pattern of A to scalar entries
  std::vector<std::set<int> > adjacencyPattern (total*dim);
  for (auto row = A.begin(); row != A.end(); ++row) {
    auto r = row.index();
    for (auto col = (*row).begin(); col != (*row).end(); ++col) {
      auto c = col.index();
      for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
          if (A[r][c][i][j] != 0.0)
            adjacencyPattern[r*dim+i].insert ((int)c*dim+j);
    }
  }
  
    // Build the adjacency pattern for the entries in B to the right of and below A
  for (auto it = gv[SLAVE].template begin<dim>(); it != gv[SLAVE].template end<dim>(); ++it) {
    IdType id = gids[SLAVE].id (*it);
    if (inactive[SLAVE].find (id) != inactive[SLAVE].end()) {
      int ib = twoMapper->mapInBoundary (SLAVE, id);
      for (int i = 0; i < dim; ++i) {
        adjacencyPattern[(n_T+ib)*dim+i].insert ((n_T+ib)*dim+i);  // Id_I
        adjacencyPattern[(n_N+n_M+ib)*dim+i].insert ((n_T+ib)*dim+i);  // D_I
      }
    } else if (active[SLAVE].find (id) != active[SLAVE].end()) {
      int ia = twoMapper->mapInActive (SLAVE, id);
      int ib = twoMapper->mapInBoundary (SLAVE, id);
      int im = twoMapper->map (SLAVE, id);
      int ii = (n_T+n_I)*dim;
      
      for (int i = 0; i < dim; ++i)
        adjacencyPattern[(n_N+n_M+n_I+ia)*dim+i].insert ((n_T+n_I+ia)*dim+i); // D_A
      
      // T_A: dim-1 tangent vectors to determine the tangent hyperplane at each node
      for (int j = 0; j < dim - 1; ++j)
        for (int k = 0; k < dim; ++k)
          adjacencyPattern[ii + (dim-1)*ia + j].insert ((n_T + ib)*dim + k);
      
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
  
  B = 0.0; // Careful! we will be using operator+= later, so this is important.
  
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
  
  cout << "A was copied.\n";
  
    // Copy D
  for (int r = 0; r < n_I+n_A; ++r)
    for (int i = 0; i < dim; ++i)
      B[(r+n_N+n_M)*dim+i][(r+n_T)*dim+i] = D[r][r][i][i];
  
  
  cout << "D was copied.\n";
  
    // Copy Id_I
  for (int r = 0; r < n_I; ++r)
    for (int i = 0; i < dim; ++i)
      B[(r+n_T)*dim+i][(r+n_T)*dim+i] = 1.0;
  
  cout << "B initialized (1/2).\n";
  
    // Fill the lines:
    //             |    0      0      0      0   T_A |
    //             |    0      0    N_A      0     0 |
  
  std::vector<int> n_d_count (n_d.size(), 0);  // count of vertices contributing to the computation of each n_d[i]

  for (auto it = gv[SLAVE].template begin<0>(); it != gv[SLAVE].template end<0>(); ++it) {
    for (auto is = gv[SLAVE].ibegin (*it) ; is != gv[SLAVE].iend (*it) ; ++is) {
      if (is->boundary () && gap[SLAVE].isSupported (*is)) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        for (int i = 0 ; i < vnum; ++i) {
          int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
          IdType id = gids[SLAVE].subId (*it, subi, dim);
          int    ib = twoMapper->mapInBoundary (SLAVE, id);
          coord_t   nr = is->centerUnitOuterNormal();
          coord_t D_nr = FMatrixHelp::mult (D[ib][ib], nr);
          
          n_d[ib] += D_nr; // FIXME: we shouldn't compute this here
          n_d_count[ib] += 1;
          
            //cout << "ib= " << ib << ", nr= " << nr << "\n";
          if (active[SLAVE].find (id) != active[SLAVE].end()) {
            int ia = twoMapper->mapInActive (SLAVE, id);
              //cout << "Found active: " << id << " -> " << ia << "\n";
              //auto ipos = is->inside()->template subEntity<dim>(subi)->geometry().center();
            
              // first the tangential stuff for the multipliers.
            std::vector<coord_t> tg = basisOfPlaneNormalTo (nr);  // dim-1 items
            
            int ii = (n_T + n_I)*dim + ia*(dim-1);
              //cout << "ii= " << ii << "\n";
            for (int j = 0; j < dim-1; ++j)
              for (int k = 0; k < dim; ++k)
                B[ii+j][(n_T+ib)*dim+k] += tg[j][k]; //FIXME: why a sum?
            
            c[ii] = 0.0;
            
              // now the last rows
            ii = (n_T + n_I)*dim + n_A*(dim-1) + ia;
            int jj = twoMapper->map (SLAVE, id)*dim;
              //cout << "ii= " << ii << ", jj= " << jj << "\n";
            for (int j = 0; j < dim; ++j)
              B[ii][jj+j] += D_nr[j]; //FIXME: why a sum?
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
  
  /*
  for (int i = 0; i < total*dim; ++i) {
    if (adjacencyPattern[i].size() > 0) {
      cout << "WTF?!?!? Columns not filled at row " << i << ": ";
      for (const auto& x : adjacencyPattern[i])
        cout << " " << x;
      cout << LF;
    }
  }
  */
  writeMatrixToMatlab (B, string ("/tmp/B") + cnt);
  writeVectorToFile (c, string ("/tmp/c") + cnt);
  
  bench().report ("Stepping", "Solving", false);
  
  try {
    InverseOperatorResult stats;
    SuperLU<ScalarMatrix> slu (B, true);
    slu.apply (uu, c, stats);
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
    exit (1);
  }
  
  bench().report ("Stepping", " done.");
  
    //// FIXME: the following stuff doesn't belong here.
  
  for (int i = 0; i < total; ++i)
    for (int j = 0; j < dim; ++j)
      u[i][j] = uu[i*dim+j];

    //// Copy results
  for (const auto& x : inactive[SLAVE]) {
    int i = twoMapper->mapInBoundary (SLAVE, x);
    int j = twoMapper->map (SLAVE, x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
  for (const auto& x : active[SLAVE]) {
    int i = twoMapper->mapInBoundary (SLAVE, x);
    int j = twoMapper->map (SLAVE, x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
    //printvector (cout, n_d, "D * Normal", "");
    //printvector (cout, n_u, "Solution Normal", "");
    //printvector (cout, n_m, "Multiplier Normal", "");
  
    //// Transform back to original coordinates

    // Is this the same as the code below with the loops?
//    CoordVector tu (u);
//    Q.mtv (tu, u);

  CoordVector tu (n_T), tu2 (n_T);
  for (int i = 0; i < n_T; ++i) tu[i] = u[i];
  Q.mtv (tu, tu2);
  for (int i = 0; i < n_T; ++i) u[i] = tu2[i];
  
  writeVectorToFile (u, string ("/tmp/u") + cnt);
}


template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::solve ()
{
  typedef TwoToOneBodyMapper<dim, TGV> MapperAdapter;
  MapperAdapter mapper_m (MASTER, *twoMapper);
  MapperAdapter mapper_s (SLAVE, *twoMapper);

  PostProcessor<TGV, TET, MapperAdapter, TSS> post_m (gv[MASTER], mapper_m, a[MASTER]);
  PostProcessor<TGV, TET, MapperAdapter, TSS> post_s (gv[SLAVE], mapper_s, a[SLAVE]);
    //TwoRefs<PostProcessor<TGV, TET, MapperAdapter, TSS> > post (post_m, post_s);

  std::string filename[2];
  filename[MASTER] = "/tmp/TwoBodiesIA_MASTER";
  filename[SLAVE] = "/tmp/TwoBodiesIA_SLAVE";
  
  CoordVector cu[2];
  DOBOTH (body) {
    cu[body].resize (gv[body].size (dim));
  }
  
  const int maxiter = 10;
  int cnt=0;
  bench().start ("Solving");
  while (true && ++cnt <= maxiter) {
    bench().start ("Active set initialization", false);
    determineActive();  // needs g, n_u, n_m computed from last iteration or =0
    bench().stop ("Active set initialization");
    bench().start ("Adjacency computation");
    setupMatrices();  // resets sparse matrix info according to new active/inactive sets
    bench().stop ("Adjacency computation");
    bench().start ("Assembly", false);
    assemble();  // reassembles matrices FIXME: should just reorder.
    bench().stop ("Assembly");
    bench().start ("Stepping");
    step (cnt);
    bench().stop ("Stepping");

    bench().start ("Postprocessing", false);
    DOBOTH (body) {
      cu[body] = 0.0;
      for (auto it = gv[body].template begin<dim>(); it != gv[body].template end<dim>(); ++it) {
        int from = twoMapper->map (body, *it);
        int to = twoMapper->mapInBody (body, *it);
//        cout << "from: " << from << ", to: " << to << LF;
        cu[body][to] = u[from];
      }
      writeVectorToFile (cu[body], string ("/tmp/cu-")+body);
    }
    (void) post_m.computeError (cu[MASTER]);
    post_m.computeVonMisesSquared ();
    post_m.writeVTKFile (filename[MASTER], cnt);
    (void) post_s.computeError (cu[SLAVE]);
    post_s.computeVonMisesSquared ();
    post_s.writeVTKFile (filename[SLAVE], cnt);

    bench().stop ("Postprocessing");
  }
  bench().stop ("Solving");
}
#endif /* defined (TWOBODIES_ACTIVESET_HPP) */

/*
void oldStuff()
{
  
 // Build transformed matrix Â:
   
 //Mindless copying all around... It would've been better to define Q
 
  if (false) {
    BlockMatrix tmp, tmp2;
    BlockMatrix AA_NN, AA_NM, AA_NS;
    BlockMatrix AA_MN, AA_MM, AA_MS;
    BlockMatrix AA_SN, AA_SM, AA_SS;
    
    subMatrix (A, AA_NN, 0,       0,       n_N, n_N);
    subMatrix (A, AA_NM, 0,       n_N,     n_N, n_M);
    subMatrix (A, AA_NS, 0,       n_N+n_M, n_N, n_S);
    subMatrix (A, AA_MN, n_N,     0,       n_M, n_N);
    subMatrix (A, AA_MM, n_N,     n_N,     n_M, n_M);
    subMatrix (A, AA_MS, n_N,     n_N+n_M, n_M, n_S);
    subMatrix (A, AA_SN, n_N+n_M, 0,       n_S, n_N);
    subMatrix (A, AA_SM, n_N+n_M, n_N,     n_S, n_M);
    subMatrix (A, AA_SS, n_N+n_M, n_N+n_M, n_S, n_S);
    
    cout << "one" << LF;
    matMultMat (tmp, AA_NS, M);
    AA_NM += tmp;  // WRONG! sparsity pattern of tmp must be a subset of that of AA_NM!!
    cout << "two" << LF;
    transposeMatMultMat (tmp, M, AA_SN);
    AA_MN += tmp;
    cout << "three" << LF;
    matMultMat (tmp, AA_MS, M);
    AA_MM += tmp;
    transposeMatMultMat (tmp, M, AA_SM);
    AA_MM += tmp;
    matMultMat (tmp, A, M);
    transposeMatMultMat (tmp2, M, tmp);
    AA_MM += tmp2;
    cout << "four" << LF;
    transposeMatMultMat (tmp, M, AA_SS);
    AA_MS += tmp;
    cout << "five" << LF;
    matMultMat (tmp, AA_SS, M);
    AA_SM += tmp;
    cout << "six" << LF;
    BCRSMatrix<BlockMatrix> big;
    
    big.setSize (3,3);
    big.setBuildMode (BCRSMatrix<BlockMatrix>::row_wise);
    for (auto row = big.createbegin(); row != big.createend(); ++row)
      for (int i=0;i<3;++i)
        row.insert (i);
    cout << "seven" << LF;
    big[0][0] = AA_NN;
    big[0][1] = AA_NM;
    big[0][2] = AA_NS;
    big[1][0] = AA_MN;
    big[1][1] = AA_MM;
    big[1][2] = AA_MS;
    big[2][0] = AA_SN;
    big[2][1] = AA_SM;
    big[2][2] = AA_SS;
    cout << "eight" << LF;
    flattenMatrix (big, A);
    cout << "nine" << LF;
    
    
    CoordVector b_M, b_S;
    b_M.resize (n_M);
    b_S.resize (n_S);
    for (int i=0; i<n_M; ++i) b_M[i] = b[i+n_N];
    for (int i=0; i<n_S; ++i) b_S[i] = b[i+n_N+n_M];
    
    M.umtv (b_S, b_M); // b_M = b_M + M^T*b_S
    
    ScalarVector uu, c;
    c.resize (total*dim, false);
    uu.resize (total*dim, false);
    uu  = 0.0;
    n_d = 0.0;
    n_m = 0.0;
    n_u = 0.0;
    
      // Copy RHS vector
    for (int i = 0; i < n_N*dim; ++i)
      for (int j = 0; j < dim; ++j)
        c[i*dim+j] = b[i][j];
    for (int i = 0; i < n_M*dim; ++i)
      for (int j = 0; j < dim; ++j)
        c[(n_N+i)*dim+j] = b_M[i][j];
    for (int i = 0; i < n_S*dim; ++i)
      for (int j = 0; j < dim; ++j)
        c[(n_N+n_M+i)*dim+j] = b_S[i][j];
    
      // Rest of entries.
    for (int i = n_T*dim; i < (n_T+n_I)*dim; ++i)
      c[i] = 0.0;
      // recall that g[] uses aiMapper.inBoundary() ordering.
    for (int i = 0; i < n_A; ++i)
      c[i+n_A+(n_T+n_I)*dim] = g[i+n_I];
    
  }
  
  
  
  
  bench().report ("Stepping", "Gluing");
  const int   n_T = static_cast<int> (gv[MASTER].size (dim) + gv[SLAVE].size (dim));
  const int   n_N = static_cast<int> (other[MASTER].size() + other[SLAVE].size());
  const int   n_M = static_cast<int> (inactive[MASTER].size() + active[MASTER].size());
  const int   n_S = static_cast<int> (inactive[SLAVE].size() + active[SLAVE].size());
  const int   n_I = static_cast<int> (inactive[SLAVE].size());
  const int   n_A = static_cast<int> (active[SLAVE].size());
  const int total = n_T + n_I + n_A;  // matrix size
  
  assert (n_T == n_N + n_M + n_S);
  
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
    // Rest of entries.
  for (int i = n_T*dim; i < (n_T+n_I)*dim; ++i)
    c[i] = 0.0;
    // recall that g[] uses aiMapper.inBoundary() ordering.
  for (int i = 0; i < n_A; ++i)
    c[i+n_A+(n_T+n_I)*dim] = g[i+n_I];
  
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
            adjacencyPattern[r*dim+i].insert ((int)c*dim+j);
    }
  }
  
    // Build the adjacency pattern for the entries in B to the right of and below A
  for (auto it = gv[SLAVE].template begin<dim>(); it != gv[SLAVE].template end<dim>(); ++it) {
    IdType id = gids[SLAVE].id (*it);
    if (inactive[SLAVE].find (id) != inactive[SLAVE].end()) {
      int ii = twoMapper->mapInBoundary (SLAVE, id);
      for (int i = 0; i < dim; ++i) {
        adjacencyPattern[(n_T+ii)*dim+i].insert ((n_T+ii)*dim+i);  // Id_I
        adjacencyPattern[(n_N+n_M+ii)*dim+i].insert ((n_T+ii)*dim+i);  // D_I
      }
    } else if (active[SLAVE].find (id) != active[SLAVE].end()) {
      int ia = twoMapper->mapInActive (SLAVE, id);
      for (int i = 0; i < dim; ++i)
        adjacencyPattern[(n_N+n_M+n_I+ia)*dim+i].insert ((n_T+n_M+n_I+ia)*dim+i); // D_A
      
      int ii = (n_T+n_I)*dim + ia;
        // T_A
      adjacencyPattern[ii].insert ((n_T+twoMapper->mapInBoundary (SLAVE, id))*dim);
      adjacencyPattern[ii].insert ((n_T+twoMapper->mapInBoundary (SLAVE, id))*dim+1);
        // N_A
      adjacencyPattern[ii+n_A].insert (twoMapper->map (SLAVE, id)*dim);
      adjacencyPattern[ii+n_A].insert (twoMapper->map (SLAVE, id)*dim+1);
    }
  }
    //cout << "Adjacency computed.\n";
  
  for (auto row = B.createbegin(); row != B.createend(); ++row)
    for (const auto& col : adjacencyPattern[row.index()])
      row.insert (col);
  
    //B = 0.0;  // necessary?
  
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
      B[(r+n_N+n_M)*dim+i][(r+n_T)*dim+i] = D[r][r][i][i];
  
    // Copy Id_I
  for (int r = 0; r < n_I; ++r)
    for (int i = 0; i < dim; ++i)
      B[(r+n_T)*dim+i][(r+n_T)*dim+i] = 1.0;
  
    //cout << "B initialized (1/2).\n";
  
    // Fill the lines:
    //             |    0      0      0      0   T_A |
    //             |    0      0    N_A      0     0 |
  
  for (auto it = gv[SLAVE].template begin<0>(); it != gv[SLAVE].template end<0>(); ++it) {
    for (auto is = gv[SLAVE].ibegin (*it) ; is != gv[SLAVE].iend (*it) ; ++is) {
      if (is->boundary ()) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        
        if (gap[SLAVE].isSupported (igeo)) {
          for (int i = 0 ; i < vnum; ++i) {
            int  subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids[SLAVE].subId (*it, subi, dim);
            int    ib = twoMapper->mapInBoundary (SLAVE, id);
            coord_t   nr = is->centerUnitOuterNormal();
            coord_t D_nr = FMatrixHelp::mult (D[ib][ib], nr);
            
              // FIXME: for nodes at the boundary of the gap the division by vnum
              // will be wrong since there won't be as many sums
            n_d[ib] += D_nr*(1.0/vnum); // FIXME: we shouldn't compute this here
            
              //cout << "ib= " << ib << ", nr= " << nr << "\n";
            if (active[SLAVE].find (id) != active[SLAVE].end()) {
              int ia = twoMapper->mapInActive (SLAVE, id);
                //cout << "Found active: " << id << " -> " << ia << "\n";
                //auto ipos = is->inside()->template subEntity<dim>(subi)->geometry().center();
              
              coord_t tg;
              tg[0] = -nr[1]; tg[1] = nr[0];
              
                // first the tangential stuff for the multipliers.
              int ii = (n_T + n_I)*dim + ia;
                //cout << "ii= " << ii << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][(n_T+ib)*dim+j] += tg[j];// *(1.0/vnum);
              c[ii] = 0.0;
              
                // now the last rows
              ii = (n_T + n_I)*dim + n_A + ia;
              int jj = twoMapper->map (SLAVE, id)*dim;
                //cout << "ii= " << ii << ", jj= " << jj << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][jj+j] += D_nr[j];// *(1.0/vnum);
            }
          }
        }
      }
    }
  }
    //cout << "B initialized (2/2).\n";
  bench().report ("Stepping", " done.");
  
  writeMatrixToMatlab (B, "/tmp/B");
  writeVectorToFile (c, "/tmp/c");
  
  bench().report ("Stepping", "Solving", false);
  
  try {
    InverseOperatorResult stats;
    SuperLU<ScalarMatrix> slu (B, true);
    slu.apply (uu, c, stats);
  } catch (Exception& e) {
    cout << "DEAD! " << e.what() << "\n";
    exit (1);
  }
  
  bench().report ("Stepping", " done.");
  
    //// FIXME: the following stuff doesn't belong here.
  
  for (int i = 0; i < total; ++i)
    for (int j = 0; j < dim; ++j)
      u[i][j] = uu[i*dim+j];
  
  CoordVector tu (u);
  Q.mtv (tu, u);
  
  for (const auto& x : inactive[SLAVE]) {
    int i = twoMapper->mapInBoundary (SLAVE, x);
    int j = twoMapper->map (SLAVE, x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
  for (const auto& x : active[SLAVE]) {
    int i = twoMapper->mapInBoundary (SLAVE, x);
    int j = twoMapper->map (SLAVE, x);
    n_u[i] = n_d[i] * u[j];
    n_m[i] = n_d[i] * u[n_T+i];
  }
    //printvector (cout, n_d, "D * Normal", "");
  printvector (cout, n_u, "Solution Normal", "");
  printvector (cout, n_m, "Multiplier Normal", "");
}
*/