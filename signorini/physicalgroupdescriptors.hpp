/******************************************************************************
 * physicalgroupdescriptors.hpp                                               *
 *                                                                            *
 ******************************************************************************/

#ifndef SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP
#define SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP

#include <vector>
#include <set>

#include <dune/grid-glue/extractors/extractorpredicate.hh>

using namespace Dune;

/*! FIXME: why level views? 
 */
template <class TGV>
class PhysicalFaceDescriptor : public ExtractorPredicate<TGV, 1>
{
  static const int dim = TGV::dimension;
  
  typedef typename TGV::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef LevelMultipleCodimMultipleGeomTypeMapper
            <typename TGV::Grid, MCMGElementLayout> ElementMapper;
  
  const TGV& gv;   //!< Grid view
  ElementMapper* mapper;
  std::vector<int> bi2pe;
  std::set<int> groups;
  
public:
  PhysicalFaceDescriptor (const TGV& _gv,
                          const std::vector<int>& boundary_id_to_physical_entity,
                          std::set<int> _groups)
  : gv (_gv), bi2pe (boundary_id_to_physical_entity), groups (_groups)
  {
    mapper = new ElementMapper (gv.grid(), 0);  // FIXME: use actual level
  }
  
  virtual ~PhysicalFaceDescriptor ()
  {
    delete mapper;
  }
  
  virtual bool contains (const EntityPointer& ep, unsigned int face) const
  {
    const int bid = mapper->map (*ep, face, 1);
    return (groups.find (bi2pe[bid]) != groups.end());
  }
};

/*! FIXME: why level views?
 */
template <class TGV>
class PhysicalVolumeDescriptor : public ExtractorPredicate<TGV, 0>
{
  static const int dim = TGV::dimension;
  
  typedef typename TGV::Traits::template Codim<0>::EntityPointer EntityPointer;
  typedef LevelMultipleCodimMultipleGeomTypeMapper
          <typename TGV::Grid, MCMGElementLayout> ElementMapper;
  
  const TGV& gv;   //!< Grid view
  ElementMapper* mapper;
  std::vector<int> ei2pe;
  std::set<int> groups;
  
public:
  PhysicalVolumeDescriptor (const TGV& _gv,
                            const std::vector<int>& element_index_to_physical_entity,
                            std::set<int> _groups)
  : gv (_gv), ei2pe (element_index_to_physical_entity), groups (_groups)
  {
    mapper = new ElementMapper (gv.grid(), 0);  // FIXME: use actual level
  }
  
  virtual ~PhysicalVolumeDescriptor ()
  {
    delete mapper;
  }
  
  virtual bool contains (const EntityPointer& ep) const
  {
    const int eid = mapper->map (*ep);
    return (groups.find (ei2pe[eid]) != groups.end());
  }
};




#endif /* defined(SIGNORINI_PHYSICALGROUPDESCRIPTORS_HPP) */
