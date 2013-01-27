/******************************************************************************
 * integration.hpp                                                            *
 *                                                                            *
 * Old stuff that I probably don't need any more...                           *
 ******************************************************************************/

#ifndef SIGNORINI_INTEGRATION_HPP
#define SIGNORINI_INTEGRATION_HPP

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace Dune;

#include <iostream>
using std::cout; using std::cerr; using std::endl; using std::string;

/*! Compute the approximate integral of a function over entity with given order
 
 We first extract a reference to a Dune::QuadratureRule from the Dune::QuadratureRules 
 singleton which is a container containing quadrature rules for all the different
 reference element types and different orders of approximation. Both classes are
 parametrized by dimension and the basic type used for the coordinate positions.
 Dune::QuadratureRule in turn is a container of Dune::QuadraturePoint supplying
 positions ξ_i and weights w_i.
 For each quadrature point i the function value at the transformed position,
 the weight and the integration element are computed and summed.
 */
template<class TEntity, class TFunctor>
double integrateEntity (const TEntity &entity, const TFunctor &f, int p)
{
  typedef typename TEntity::ctype ctype;
  const int                             dim = TEntity::dimension;
  const typename TEntity::Geometry geometry = entity.geometry();
  const GeometryType                     gt = geometry.type();
  
    // get quadrature rule of order p
  const auto& rule = QuadratureRules<ctype, dim>::rule (gt, p);
  
    // ensure that rule has at least the requested order
  if (rule.order() < p)
    DUNE_THROW(Exception,"order not available");
  
    // compute approximate integral looping over all quadrature points in the quadrature rule
  double result=0;
  for (auto i = rule.begin(); i != rule.end(); ++i) {
    double fval   = f (geometry.global (i->position()));
    double weight = i->weight();
    double detjac = geometry.integrationElement (i->position());
    result       += fval * weight * detjac;
  }

  return result;
}

/*! Parameter predicate for mapper class.
 
 Includes entities with codimension zero.
 
 This class is only here to show what such a class looks like -- it does
 exactly the same as Dune::MCMGElementLayout.
 */
template <int dimgrid>
struct CubeLayout
{
  bool contains (Dune::GeometryType gt) { return (gt.dim() == dimgrid); }
};

/*! Parameter predicate for mapper class.
 
 Includes entities with dimension zero.
 
 This class is only here to show what such a class looks like -- it does
 exactly the same as Dune::MCMGVertexLayout.
 */
template <int dimgrid>
struct VertexLayout
{
  bool contains (Dune::GeometryType gt) { return (gt.dim() == 0); }
};


template <template <int> class TPolicy, class TGrid>
LeafMultipleCodimMultipleGeomTypeMapper<TGrid, TPolicy>* makeMapper(TGrid& grid)
{
  return new LeafMultipleCodimMultipleGeomTypeMapper<TGrid, TPolicy> (grid);
}

template <typename TFunctor>
const string zeropad (int num, int max)
{
  std::stringstream tmp;
  const static int width = ceil (log10 (++max));  // max number of digits
  tmp << "integrate" << TFunctor::name() << std::setfill('0') << std::setw(width) << num;
  shared_ptr<string> ss(new string(tmp.str()));
  return *ss;
}


/*! Uniform refinement test.
 
 Integrate the given functor in the grid and refine the given number of steps.
 For each step produce an output file.
 
 NEWBIE WARNING: if using the CubeLayout predicate for the mapper, then to 
 inspect the array "integral" in ParaView, use something like
 
   pdi = self.GetInput()
   sol = pdi.GetCellData().GetArray('solution')
 
 NOT GetPointData(), because the mapper attaches data to cells, not points. ok?
 
 FIXME: use an id-based mapper to transfer data from one grid to the next one.
 */
template<class TFunctor, class TGrid>
double integrateGridTestError (TGrid& grid, int steps, int order)
{
  const TFunctor f;
  const int dim = TGrid::dimension;

  /* //no c++11 => no lambdas :(
  static std::stringstream basename ("integrateGridTestError");
  auto zeropad = [&] (int num) -> const std::ostream& {
    const static int width = ceil (log10 (steps+1));
    return basename << std::setfill('0') << std::setw(width) << num;
  };
  */
  auto gridView = grid.leafView();
  double value = 0, oldvalue = 1E100;
  for (int k = 0; k < steps; k++, value = 0.0) {
    gridView = grid.leafView();
      // create mapper for codim0 entities and allocate vector for data
    LeafMultipleCodimMultipleGeomTypeMapper<TGrid, VertexLayout>* mapper;
    mapper = makeMapper<VertexLayout> (grid);
    std::vector<double> c;
    c.resize (mapper->size());
      //cout << "(Re)Created std::vector of size: " << c.size() << endl;
    
       // compute integral with some order
    for (auto it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it) {
      double tmp = integrateEntity (*it, f, order);
      
        //auto g = it->geometry();
        //cout << "Traversing " << it->template count<dim> () << " vertices for " << g.type() << " at " << g.corner(0) << endl;
      
      for (int i=0; i < it->template count<dim> (); ++i) {
          //auto ep = it->template subEntity<dim>(i);  // EntityPointer
          //cout << "Storing value for " <<  ep->geometry().type() << " at iteration " << i << " with pos " << ep->geometry().corner(0) << " at vector offset: "  << mapper->map(*ep) << endl;
        c[mapper->map(*(it->template subEntity<dim>(i)))] = tmp;
      }
      
      value += tmp;
        //cout << "\tcomputed= " << tmp << "\n\t\ttotal= " << value << endl;
    }
      // print result and error estimate
    cout << "elements= " << std::setw(8) << std::right << grid.size(0)
         << " integral= " << std::scientific << std::setprecision(12) << value
         << " error= " << std::abs(value-oldvalue) << endl;

      // generate a VTK file
    VTKWriter<typename TGrid::LeafGridView> vtkwriter(gridView);
    vtkwriter.addVertexData (c, "integral");
    vtkwriter.write (zeropad<TFunctor> (k, steps), VTKOptions::binaryappended);
    
      // save result for the next estimate and refine the grid
    oldvalue = value;
    grid.globalRefine(1);
    
    delete mapper;
    c.empty();
  }
  
  return value;
}

#endif
