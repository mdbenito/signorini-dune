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
  // HACK: Iterate thorugh arrays for both bodies
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
    //typedef BCRSMatrix<coord_t>                 VectorMatrix; // WRONG
  typedef BCRSMatrix<block_t>                  BlockMatrix;
  typedef BlockVector<scalar_t>               ScalarVector;
  typedef BlockVector<coord_t>                 CoordVector;
  typedef ActiveSetFunctor<ctype, dim, ScalarVector>  AIFunctor;
  typedef TwoBodyMapper<dim, TGV>                TwoMapper;
  typedef LeafMultipleCodimMultipleGeomTypeMapper<typename TGV::Grid,
                                                MCMGVertexLayout> VertexMapper;
  typedef FunctorSupportMapper<dim, TGV, TGF> GapVertexMapper;
  
private:
  TwoRefs<TGV> gv;    //!< Grid view
  TwoRefs<GlobalIdSet> gids; //!<
  
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
  IdSet         other[2];
  TwoMapper*   twoMapper;

public:
  TwoBodiesIASet (const TGV& _mgv, const TGV& _sgv, const TET& _ma, const TET& _sa,
                  const TFT& _mf, const TFT& _sf, const TDF& _md, const TDF& _sd,
                  const TTF& _mp, const TTF& _sp, const TGF& _mgap, const TGF& _sgap,
                  int _quadratureOrder = 4);
  
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
  assert (dim == 2);
  
  I = 0.0;
  for (int i=0; i < dim; ++i)
    I[i][i] = 1.0;
  
    //// Initialize the mapper setting all nodes in the gap to inactive:

  DOBOTH (body) {
    for (auto it = gv[body].template begin<dim>(); it != gv[body].template end<dim>(); ++it) {
      auto id = gids[body].id (*it);
      if (gap[body].isSupported (it->geometry())) inactive[body] << id;
      else                                           other[body] << id;
    }
  }
  /*  Really complicated way of doing the very same thing as above:
   Traversing the boundary brings little benefit when dim=2 and we have so few nodes.
  DOBOTH(body) {
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        if (is->boundary ()) {
          const auto& igeo = is->geometry();
          const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
          const int ivnum = ref.size (is->indexInInside (), 1, dim);
          if (gap[body].isSupported (igeo)) {
            for (int i_m = 0 ; i_m < ivnum; ++i_m) {
              int subi = ref.subEntity (is->indexInInside (), 1, i_m, dim);
              IdType id = gids[body].subId (*it, subi, dim);
              inactive[body].insert(id);
              other[body].erase(id);
            }
          } else {
            for (int i_m = 0 ; i_m < ivnum; ++i_m) {
              int subi = ref.subEntity (is->indexInInside (), 1, i_m, dim);
              IdType id = gids[body].subId (*it, subi, dim);
              if (inactive[body].find(id) == inactive[body].end())
                other[body].insert(id);
            }
          }
        }
      }
    }
  }
  */
  
  /*
  cout << "SizeI = " << inactive[MASTER].size() << LF;
  cout << "SizeI = " << inactive[SLAVE].size() << LF;
  cout << "SizeA = " << active[MASTER].size() << LF;
  cout << "SizeA = " << active[SLAVE].size() << LF;
  cout << "SizeO = " << other[MASTER].size() << LF;
  cout << "SizeO = " << other[SLAVE].size() << LF;
*/
  
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
  const auto& multBasis = TLM::instance();
  
  const int n_Tm = gv[MASTER].size (dim);
  const int n_Ts = gv[SLAVE].size (dim);
  const int  n_T = n_Tm + n_Ts;
  const int n_Am = static_cast<int>(active[MASTER].size());
  const int n_As = static_cast<int>(active[SLAVE].size());
  const int n_Im = static_cast<int>(inactive[MASTER].size());
  const int n_Is = static_cast<int>(inactive[SLAVE].size());
  const auto ingap_m = n_Am + n_Im;
  const auto ingap_s = n_As + n_Is;

    //cout << "total= " << n_T << ", ingap= " << ingap << "\n";
  
  std::vector<std::set<int> > adjacencyPattern (n_T);
  
    //For each element we traverse all its vertices and set them as adjacent
  DOBOTH (body) {
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
  
    /* Adjacency pattern for M: the mindless approach.

     Traverse the boundary of the slave, and at each boundary element in the gap
     traverse the boundary of the master looking for basis functions whose support
     ntersects.
     FIXME: we check whether the support of a basis function in the master side
     contains the coordinates of a node in the slave side, but we should rather
     check at the *quadrature points* in the slave side to avoid marking pairs
     slave node / master node which actually won't count
     
     NOTE: for regular SGrids we could simply traverse all nodes of the slave
     then all of the master and if h is the grid width parameter, then
     if |x_0-x_1|+|y_0-y_1| <= h, (the supports of the basis functions intersect)
     add the entry.
     
     Should be the same.
     
     NOTE2: We should be able to exploit this computation later when actually
     building M in order not to have to traverse all the grid (squared) again, 
     but I just don't feel like it...
     */
  bench().start ("Adjacency for coupling matrix", false);
  std::vector<std::set<int> > couplingPattern (ingap_s);
  for (auto it_s = gv[SLAVE].template begin<0>(); it_s != gv[SLAVE].template end<0>(); ++it_s) {
    for (auto is_s = gv[SLAVE].ibegin (*it_s) ; is_s != gv[SLAVE].iend (*it_s) ; ++is_s) {
      if (is_s->boundary ()) {
        const auto& igeo_s = is_s->geometry ();
        if (gap[SLAVE].isSupported (igeo_s)) {
          const auto& ref_s = GenericReferenceElements<ctype, dim>::general (it_s->type());
          
          for (auto it_m = gv[MASTER].template begin<0>(); it_m != gv[MASTER].template end<0>(); ++it_m) {
            for (auto is_m = gv[MASTER].ibegin (*it_m) ; is_m != gv[MASTER].iend (*it_m) ; ++is_m) {
              if (is_m->boundary ()) {
                const auto& igeo_m = is_m->geometry();
                if (gap[MASTER].isSupported (igeo_m)) {
                  const auto& ref_m = GenericReferenceElements<ctype, dim>::general (it_m->type());
                  
                  const int ivnum_s = ref_s.size (is_s->indexInInside (), 1, dim);
                  const int ivnum_m = ref_m.size (is_m->indexInInside (), 1, dim);
                  
                  for (int i_s = 0 ; i_s < ivnum_s; ++i_s) {
                    int subi_s = ref_s.subEntity (is_s->indexInInside (), 1, i_s, dim);
                    
                    for (int i_m = 0 ; i_m < ivnum_m; ++i_m) {
                      int subi_m = ref_m.subEntity (is_m->indexInInside (), 1, i_m, dim);
                      auto   v_s = it_s->template subEntity<dim> (subi_s)->geometry().center();
                        //const auto& global = igeo_s.global (v_s); // center() is already global?
                      const auto& local_s = it_s->geometry().local (v_s);
                      if (multBasis[subi_m].isSupported (local_s)) {
                        const int ii_s = twoMapper->mapInBoundary (SLAVE, *it_s, subi_s, dim);
                        const int ii_m = twoMapper->mapInBoundary (MASTER, *it_m, subi_m, dim);
                        couplingPattern[ii_s].insert(ii_m);
                      }
                    }
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
  
  b.resize (n_T + n_Is + n_As, false);
  u.resize (n_T + n_Is + n_As, false);
  
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
  
  for (auto it_s = gv[SLAVE].template begin<dim>(); it_s != gv[SLAVE].template end<dim>(); ++it_s) {
    auto id_s = gids[SLAVE].id (*it_s);
    if (gap[SLAVE].isSupported (it_s->geometry())) {
      if (contact.isSupported (twoMapper->mapInBoundary (SLAVE, id_s)))
        active[SLAVE] << id_s;
      else
        inactive[SLAVE] << id_s;
    }
  }
  
  twoMapper->update (active, inactive, other);
  
  cout << "\nInactive: " << inactive[SLAVE].size();
    //for (auto x : inactive) cout << x << "->" << aiMapper->map(x) << " ";
  cout << "\nActive: " << active[SLAVE].size();
    //for (auto x : active)   cout << x << "->" << " " << aiMapper->map(x) << " ";
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
              << "(" << geo.corner(i) << ") x (" << geo.corner (j) << ")\n";
            }
          }
        }
          //// Integrand of the RHS
        
        for (int i = 0 ; i < vnum; ++i)
          b[twoMapper->map (body, *it, i, dim)] += f[body] (it->geometry ().global (x.position ())) *
                                                   basis[i].evaluateFunction (x.position ()) *
                                                   x.weight () *
                                                   it->geometry ().integrationElement (x.position ());
      }
      
        //// Neumann Boundary conditions.
      
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        if (is->boundary ()) {
          const int  ivnum = ref.size (is->indexInInside (), 1, dim);
          const auto& igeo = is->geometry ();
          const auto& ityp = is->type ();
          
          if (p[body].isSupported (igeo)) {
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
    }
    
    /* Dirichlet boundary conditions.
     For unvisited vertices on the boundary, replace the associated line of A
     and b with a trivial one.
     */
    
    for (auto it = gv[body].template begin<0>(); it != gv[body].template end<0>(); ++it) {
      const auto& ref = GenericReferenceElements<ctype, dim>::general (it->type());
      
      for (auto is = gv[body].ibegin (*it) ; is != gv[body].iend (*it) ; ++is) {
        if (is->boundary ()) {
          const int ivnum = ref.size (is->indexInInside (), 1, dim);
            //cout << "Dirichlet'ing: "; printCorners (is->geometry ());
          
          for (int i = 0; i < ivnum; ++i) {
            auto subi = ref.subEntity (is->indexInInside (), 1, i, dim);
            auto global = it->geometry().global (ref.position (subi, dim));
            int ii = twoMapper->map (body, *it, subi , dim);
            if (dir[body].isSupported (global)) {
                //cout << "Dirichlet'ing node: " << ii << " at " << v << "\n";
              A[ii] = 0.0;
              A[ii][ii] = I;
              b[ii] = dir[body] (global);
            }
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
  
  for (auto it = gv[SLAVE].template begin<0>(); it != gv[SLAVE].template end<0>(); ++it) {
    for (auto is = gv[SLAVE].ibegin (*it) ; is != gv[SLAVE].iend (*it) ; ++is) {
      if (is->boundary ()) {
        const auto&   in = is->inside();
        const auto&  ref = GenericReferenceElements<ctype, dim>::general (in->type());
        const int   vnum = ref.size (is->indexInInside (), 1, dim);
        const auto& igeo = is->geometry ();
        
        if (gap[SLAVE].isSupported (igeo)) {
          for (int i = 0 ; i < vnum; ++i) {
            int subi  = ref.subEntity (is->indexInInside (), 1, i, dim);
            IdType id = gids[SLAVE].subId (*it, subi, dim);
            int ib = twoMapper->mapInBoundary (SLAVE, id);
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
  }
  
    //// Compute mortar coupling matrix MM
  
  for (auto it_s = gv[SLAVE].template begin<0>(); it_s != gv[SLAVE].template end<0>(); ++it_s) {
    for (auto is_s = gv[SLAVE].ibegin (*it_s) ; is_s != gv[SLAVE].iend (*it_s) ; ++is_s) {
      if (is_s->boundary ()) {
        const auto& igeo_s = is_s->geometry ();
        if (gap[SLAVE].isSupported (igeo_s)) {
          const auto& ref_s = GenericReferenceElements<ctype, dim>::general (it_s->type());
          
          for (auto it_m = gv[MASTER].template begin<0>(); it_m != gv[MASTER].template end<0>(); ++it_m) {
            for (auto is_m = gv[MASTER].ibegin (*it_m) ; is_m != gv[MASTER].iend (*it_m) ; ++is_m) {
              if (is_m->boundary ()) {
                const auto& igeo_m = is_m->geometry();
                if (gap[MASTER].isSupported (igeo_m)) {
                  const auto& ref_m = GenericReferenceElements<ctype, dim>::general (it_m->type());
                  
                  const int ivnum_s = ref_s.size (is_s->indexInInside (), 1, dim);
                  const int ivnum_m = ref_m.size (is_m->indexInInside (), 1, dim);
                  
                  for (int i_s = 0 ; i_s < ivnum_s; ++i_s) {
                    int subi_s = ref_s.subEntity (is_s->indexInInside (), 1, i_s, dim);
                    
                    for (int i_m = 0 ; i_m < ivnum_m; ++i_m) {
                      int subi_m = ref_m.subEntity (is_m->indexInInside (), 1, i_m, dim);
                      auto   v_s = it_s->template subEntity<dim> (subi_s)->geometry().center();
                        //const auto& global = igeo_s.global (v_s); // center() is already global?
                      const auto& local_s = it_s->geometry().local (v_s);
                      if (multBasis[subi_m].isSupported (local_s)) {
                        const int ii_s = twoMapper->mapInBoundary (SLAVE, *it_s, subi_s, dim);
                        const int ii_m = twoMapper->mapInBoundary (MASTER, *it_m, subi_m, dim);
                        for (auto& x : QuadratureRules<ctype, dim-1>::rule (is_s->type(), quadratureOrder)) {
                          const auto&   global = igeo_s.global (x.position ());
                          const auto&  local_s = it_s->geometry().local (global);
                          const auto&  local_m = it_m->geometry().local (global);
                            //cout << "    local_s= " << local_s
                            // << "    local_m= " << local_m << LF;
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
      }
    }
  }
  
    //printmatrix (cout, D, "D", "");
}


/*! HACK of the HACKS...
 
 This is terrible: we assembled a block matrix and now we flatten it. Should've
 started with a scalar matrix!!!
 
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
 B should be
 
     /                                        \
     | Â_NN   Â_NM   Â_NI   Â_NS     0      0 |
     | Â_MN   Â_MM   Â_MI   Â_MA     0      0 |
     | Â_IN   Â_IM   Â_II   Â_IA    D_I     0 |
     | Â_AN   Â_AM   Â_AI   Â_AA     0    D_A |
 B = |    0      0      0      0   Id_I     0 |
     |    0      0      0      0      0   T_A |
     |    0      0      0    N_A      0     0 |
     \                                        /
 
 In TWO DIMENSIONS this is a SQUARE scalar matrix!
 
 */
template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::step ()
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

  cout << "n_T= " << n_T << ", n_N= " << n_N << ", n_M= " << n_M
       << ", n_S= " << n_S << ", n_I= " << n_I << ", n_A= " << n_A << LF;
  
    // M = D^-1*MM
  BlockMatrix M;
  BlockMatrix DD(D);
  for (int i=0; i < DD.N(); ++i)
    DD[i][i].invert();
  matMultMat(M, DD, MM);
  
  BlockMatrix MT;
  transpose (MT, M);
  cout << "Tranposed MT is " << MT.N() << " x " << MT.M() << LF;
  BlockMatrix Q;
  Q.setSize (n_T, n_T);
  Q.setBuildMode (BlockMatrix::row_wise);
  
  for (auto row = Q.createbegin(); row != Q.createend(); ++row) {
    auto r = row.index();
    row.insert(r);
    if (n_N <= r && r < n_N+n_M)
      for (auto col = MT[r-n_N].begin(); col != MT[r-n_N].end(); ++col)
        row.insert (col.index()+n_N+n_M);
  }

  for (int r = 0; r < n_T; ++r) {
    Q[r][r] = I;
    if (n_N <= r && r < n_N+n_M)
      for (auto col = MT[r-n_N].begin(); col != MT[r-n_N].end(); ++col)
        Q[r][col.index()+n_N+n_M] = MT[r-n_N][col.index()];
  }

  writeMatrixToMatlab(A, "/tmp/A");
  writeMatrixToMatlab(D, "/tmp/D");
  writeMatrixToMatlab(DD, "/tmp/DD");
  writeMatrixToMatlab(Q, "/tmp/Q");
  writeMatrixToMatlab(M, "/tmp/M");
  writeMatrixToMatlab(MM, "/tmp/MM");
  writeMatrixToMatlab(MT, "/tmp/MT");

  BlockMatrix AA;
  matMultTransposeMat(AA, A, Q);
  matMultMat(A, Q, AA);  // A is now the new Â.
  
  cout << "We have Â\n";
  ScalarVector uu, c;
  c.resize (total*dim, false);
  uu.resize (total*dim, false);
  uu  = 0.0;
  n_d = 0.0;
  n_m = 0.0;
  n_u = 0.0;
  
  CoordVector bb(b);
  Q.mv(bb, b);  // b is now the new ^b=Q*b
  
    // Copy RHS vector
  for (int i = 0; i < n_T; ++i)
    for (int j = 0; j < dim; ++j)
      c[i*dim+j] = b[i][j];
  
    // Rest of entries.
  for (int i = n_T*dim; i < (n_T+n_I)*dim; ++i)
    c[i] = 0.0;
    // recall that g[] uses twoMapper.inBoundary(SLAVE,...) ordering.
  for (int i = 0; i < n_A; ++i) {
    cout << " c[" << (n_T+n_I)*dim + n_A + i << "]= g[" << i+n_I << "]= " << g[i+n_I] << LF;
    c[(n_T+n_I)*dim + n_A + i] = g[i+n_I];
  }
  ScalarMatrix B;
  B.setBuildMode (ScalarMatrix::row_wise);
  B.setSize (total*dim, total*dim);
  cout << " B is " << B.N() << " x " << B.M() << "\n";
  
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
        adjacencyPattern[(n_N+n_M+n_I+ia)*dim+i].insert((n_T+n_I+ia)*dim+i); // D_A
      
      int ii = (n_T+n_I)*dim + ia;
        // T_A
      adjacencyPattern[ii].insert ((n_T+twoMapper->mapInBoundary (SLAVE, id))*dim);
      adjacencyPattern[ii].insert ((n_T+twoMapper->mapInBoundary (SLAVE, id))*dim+1);
        // N_A
      adjacencyPattern[ii+n_A].insert (twoMapper->map (SLAVE, id)*dim);
      adjacencyPattern[ii+n_A].insert (twoMapper->map (SLAVE, id)*dim+1);
    }
  }
  cout << "Adjacency computed.\n";
  
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
                B[ii][(n_T+ib)*dim+j] += tg[j];//*(1.0/vnum);
              c[ii] = 0.0;
              
                // now the last rows
              ii = (n_T + n_I)*dim + n_A + ia;
              int jj = twoMapper->map (SLAVE, id)*dim;
                //cout << "ii= " << ii << ", jj= " << jj << "\n";
              for (int j = 0; j < dim; ++j)
                B[ii][jj+j] += D_nr[j];//*(1.0/vnum);
            }
          }
        }
      }
    }
  }
  cout << "B initialized (2/2).\n";
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
    exit(1);
  }
  
  bench().report ("Stepping", " done.");
  
    //// FIXME: the following stuff doesn't belong here.
  
  for (int i = 0; i < total; ++i)
    for (int j = 0; j < dim; ++j)
      u[i][j] = uu[i*dim+j];

  CoordVector tu(u);
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


template<class TGV, class TET, class TFT, class TDF, class TTF, class TGF, class TSS, class TLM>
void TwoBodiesIASet<TGV, TET, TFT, TDF, TTF, TGF, TSS, TLM>::solve ()
{
  typedef TwoToOneBodyMapper<dim, TGV> MapperAdater;
  MapperAdater mapper_m (MASTER, *twoMapper);
  MapperAdater mapper_s (SLAVE, *twoMapper);
  PostProcessor<TGV, TET, MapperAdater, TSS> post_m (gv[MASTER], mapper_m, a[MASTER]);
  PostProcessor<TGV, TET, MapperAdater, TSS> post_s (gv[SLAVE], mapper_s, a[SLAVE]);
    //TwoRefs<PostProcessor<TGV, TET, MapperAdater, TSS> > post(post_m, post_s);
  
  VertexMapper defaultMapper_m (gv[MASTER].grid());
  VertexMapper defaultMapper_s (gv[SLAVE].grid());
  TwoRefs<VertexMapper> defaultMapper (defaultMapper_m, defaultMapper_s);
  
  std::string filename[2];
  filename[MASTER] = "/tmp/TwoBodiesIA_MASTER";
  filename[SLAVE] = "/tmp/TwoBodiesIA_SLAVE";
  
  CoordVector cu[2];
  DOBOTH (body) {
    cu[body].resize(gv[body].size(dim));
  }
  
  const int maxiter = 10;
  int cnt=0;
  while (true && ++cnt < maxiter) {
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
    step();
    bench().stop ("Stepping");

    bench().start ("Postprocessing", false);
    DOBOTH (body) {
      cu[body] = 0.0;
      for (auto it = gv[body].template begin<dim>(); it != gv[body].template end<dim>(); ++it) {
        int from = twoMapper->map (body, *it);
        int to = defaultMapper[body].map (*it);
          //cout << "Mapping vertex " << from << " to " << to << LF;
        cu[body][to] = u[from];  // FIXMMEEEEEE: is this ok?
      }
    }
    (void) post_m.computeError (cu[MASTER]);
    post_m.computeVonMisesSquared ();
    (void) post_m.writeVTKFile (filename[MASTER], cnt);
    (void) post_s.computeError (cu[SLAVE]);
    post_s.computeVonMisesSquared ();
    (void) post_s.writeVTKFile (filename[SLAVE], cnt);

    bench().stop ("Postprocessing");

  }
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
    
    subMatrix(A, AA_NN, 0,       0,       n_N, n_N);
    subMatrix(A, AA_NM, 0,       n_N,     n_N, n_M);
    subMatrix(A, AA_NS, 0,       n_N+n_M, n_N, n_S);
    subMatrix(A, AA_MN, n_N,     0,       n_M, n_N);
    subMatrix(A, AA_MM, n_N,     n_N,     n_M, n_M);
    subMatrix(A, AA_MS, n_N,     n_N+n_M, n_M, n_S);
    subMatrix(A, AA_SN, n_N+n_M, 0,       n_S, n_N);
    subMatrix(A, AA_SM, n_N+n_M, n_N,     n_S, n_M);
    subMatrix(A, AA_SS, n_N+n_M, n_N+n_M, n_S, n_S);
    
    cout << "one" << LF;
    matMultMat (tmp, AA_NS, M);
    AA_NM += tmp;  // WRONG! sparsity pattern of tmp must be a subset of that of AA_NM!!
    cout << "two" << LF;
    transposeMatMultMat (tmp, M, AA_SN);
    AA_MN += tmp;
    cout << "three" << LF;
    matMultMat(tmp, AA_MS, M);
    AA_MM += tmp;
    transposeMatMultMat(tmp, M, AA_SM);
    AA_MM += tmp;
    matMultMat(tmp, A, M);
    transposeMatMultMat(tmp2, M, tmp);
    AA_MM += tmp2;
    cout << "four" << LF;
    transposeMatMultMat(tmp, M, AA_SS);
    AA_MS += tmp;
    cout << "five" << LF;
    matMultMat(tmp, AA_SS, M);
    AA_SM += tmp;
    cout << "six" << LF;
    BCRSMatrix<BlockMatrix> big;
    
    big.setSize(3,3);
    big.setBuildMode(BCRSMatrix<BlockMatrix>::row_wise);
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
    b_M.resize(n_M);
    b_S.resize(n_S);
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
            adjacencyPattern[r*dim+i].insert((int)c*dim+j);
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
        adjacencyPattern[(n_N+n_M+n_I+ia)*dim+i].insert((n_T+n_M+n_I+ia)*dim+i); // D_A
      
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
    exit(1);
  }
  
  bench().report ("Stepping", " done.");
  
    //// FIXME: the following stuff doesn't belong here.
  
  for (int i = 0; i < total; ++i)
    for (int j = 0; j < dim; ++j)
      u[i][j] = uu[i*dim+j];
  
  CoordVector tu(u);
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