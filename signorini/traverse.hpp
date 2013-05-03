/******************************************************************************
 * traverse.hpp                                                               *
 *                                                                            *
 * More old stuff that I probably don't need any more...                      *
 ******************************************************************************/

#ifndef SIGNORINI_TRAVERSE_HPP
#define SIGNORINI_TRAVERSE_HPP

using namespace Dune;

#include <iostream>
using std::cout; using std::cerr; using std::endl;

#include "integration.hpp"

/* Helper function */
template <class TGeometry>
void printCorners (const TGeometry& g)
{
  cout << " visiting " << g.type() << " with vertices at (" << g.corner(0) << ")";
  
  for (int i = 1; i < g.corners(); ++i)
    cout << ", (" << g.corner(i) << ")";
  
  cout << endl;
}


#define THERE_ARE(x) cout << "There " << (((x) != 1) ? "are " : "is ") << (x)

/*! Traverse the entities of a given mesh in various ways */
template <class TGrid>
void printTraversalInfo (TGrid& grid)
{
  int cnt;
  const int dim = TGrid::dimension;
  
    ///// Leaf Traversal of co-dimension 0
  
  cout << "*** Traversing codim 0 leaves:" << endl;
  
  auto leafView = grid.leafView(); // get the instance of the LeafGridView
  cnt = 0;
  for (auto it = leafView.template begin<0>(); it != leafView.template end<0>(); ++it, ++cnt)
    printCorners (it->geometry());
  
  THERE_ARE(cnt) << " leaf element(s)" << endl;
  
    ///// Leafwise traversal of co-dimension dim (vertices)
  
  cout << endl << "*** Traversing codim " << dim << " leaves:" << endl;
  cnt = 0;
  for (auto it = leafView.template begin<dim>(); it != leafView.template end<dim>(); ++it, ++cnt)
    printCorners (it->geometry());
  
  THERE_ARE(cnt) << " leaf vertices(s)" << endl;
  
    ///// Levelwise traversal of codim 0
  
  cout << endl << "*** Traversing codim 0 for each level:" << endl;
  
  for (int level=0; level <= grid.maxLevel(); ++level) {
    auto levelView = grid.levelView (level);
       // iterate through all entities of codim 0 on the given level
    cnt = 0;
    for (auto it = levelView.template begin<0>(); it != levelView.template end<0>(); ++it, ++cnt)
      printCorners (it->geometry());
    
    THERE_ARE(cnt) << " element(s) on level " << level << endl;
  }
}
#undef THERE_ARE

template <class TFunctor, class TGrid>
void printTraversalQuadratures (TGrid& grid)
{
  auto leafView = grid.leafView();
  const TFunctor f;

  for (auto it = leafView.template begin<0>(); it != leafView.template end<0>(); ++it) {
    cout << "Quadrature for " << it->geometry().type()
         << " with first vertex at ("
         << std::fixed  << std::setprecision(2) << it->geometry().corner(0) << ") = "
         << std::scientific << std::setprecision(12) << integrateEntity(*it, f, 1) << endl;
  }
}

#endif  //SIGNORINI_TRAVERSE_HPP
