/******************************************************************************
 * physicalgroupdescriptors.hpp                                               *
 *                                                                            *
 ******************************************************************************/

#ifndef SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP
#define SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP

#include <iostream>
#include <vector>
#include <set>

#include <dune/grid-glue/extractors/extractorpredicate.hh>

using namespace Dune;
using std::cout;
using std::vector;
using std::set;


/*! FIXME: why level views? 
 */
template <class TGV, class TGF>
class PhysicalFaceDescriptor : public ExtractorPredicate<TGV, 1>
{
  static const int dim = TGV::dimension;
  
  typedef typename TGV::Traits::template Codim<0>::EntityPointer EntityPointer;
  
  const TGF&     gf;   //!< Grid factory
  vector<int> bi2pe;
  set<int>   groups;
  
public:
  PhysicalFaceDescriptor (const TGF& _gf,
                          const vector<int>& boundary_id_to_physical_entity,
                          const set<int>& _groups)
  : gf (_gf), bi2pe (boundary_id_to_physical_entity), groups (_groups)
  { }
  
  virtual bool contains (const EntityPointer& ep, unsigned int face) const
  {
//      cout << "     ep has boundary intersections. Looking for face " << face << LF;
    for (auto is = ep->ileafbegin(); is != ep->ileafend(); ++is) {
        //        cout << "          intersection #" << is->indexInInside();
      if (is->indexInInside() == face && is->boundary()) {
        unsigned int idx = gf.insertionIndex (*is);
          //          cout << "               was inserted as " << idx << LF;
        return (idx < bi2pe.size() && groups.find (bi2pe[idx]) != groups.end());
      } //else cout << LF;
    }
    return false;
  }
};

/*! FIXME: why level views?
 */
template <class TGV, class TGF>
class PhysicalVolumeDescriptor : public ExtractorPredicate<TGV, 0>
{
  static const int dim = TGV::dimension;
  
  typedef typename TGV::Traits::template Codim<0>::EntityPointer EntityPointer;
  
  const TGF&     gf;   //!< Grid factory
  vector<int> ei2pe;
  set<int>   groups;
  
public:
  PhysicalVolumeDescriptor (const TGF& _gf,
                            const vector<int>& element_index_to_physical_entity,
                            const set<int>& _groups)
  : gf (_gf), ei2pe (element_index_to_physical_entity), groups (_groups)
  { }

  virtual bool contains (const EntityPointer& ep) const
  {
    const auto idx = gf.insertionIndex (*ep);
    return (idx >=0 && idx < ei2pe.size() && groups.find (ei2pe[idx]) != groups.end());  // FIXME: ok? test!
  }
};

#endif /* defined(SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP) */
