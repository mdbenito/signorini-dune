/******************************************************************************
 * activeinactivemapper.hpp                                                   *
 ******************************************************************************/

#ifndef SIGNORINI_ACTIVEINACTIVEMAPPER
#define SIGNORINI_ACTIVEINACTIVEMAPPER

#include <map>
#include <dune/grid/common/mapper.hh>

/*! A mapper to order entities according to their status as active or inactive.
 
 The constructor is given three sets of entities representing those which are
 active "A", inactive "I" and the rest "N". The final ordering returned by the
 mapper will be:
 
    N_1 ... N_r I_1 ... I_m A_1 ... A_n 
 
 Template parameters:
 TGV: TGridView
 */
template <int codim, class TGV>
class ActiveInactiveMapper
: public Mapper<typename TGV::Grid, ActiveInactiveMapper<codim, TGV> >
{
  typedef typename TGV::template Codim<codim>::Entity Entity;
  typedef typename TGV::EntityPointer EntityPointer;
  typedef std::vector<EntityPointer> EntityPointerVector;
  
  const TGV&      gv;
  const typename TGV::IndexSet& iset;
  int*     indices;

public:
  ActiveInactiveMapper (const TGV& _gv,
                        const EntityPointerVector& active,
                        const EntityPointerVector& inactive,
                        const EntityPointerVector& others);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (const Entity& e, int i, unsigned int cc) const;
  
  bool contains (const Entity& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;
  
  int size () const;
  void update (const EntityPointerVector& active,
               const EntityPointerVector& inactive,
               const EntityPointerVector& others);
};


template <int codim, class TGV>
ActiveInactiveMapper<codim, TGV>::ActiveInactiveMapper (const TGV& _gv,
                                                        const EntityPointerVector& active,
                                                        const EntityPointerVector& inactive,
                                                        const EntityPointerVector& others)
: gv (_gv), iset (_gv.indexSet()), indices (NULL)
{
  update (active, inactive, others);
}


template <int codim, class TGV>
template<class EntityType>
int ActiveInactiveMapper<codim, TGV>::map (const EntityType& e) const
{
    //int idx = indices[iset.index (e)];
    //return (idx < 0) ? (-1*idx + gapEntities.size()-1) : idx;
  return indices[iset.index (e)];
}


template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (const Entity& e, int i,
                                           unsigned int cc) const
{
    //int idx = indices[iset.index (e, i, cc)];
    //return (idx < 0) ? (-1*idx + gapEntities.size()-1) : idx;
  indices[iset.index (e, i, cc)];
}


template <int codim, class TGV>
bool ActiveInactiveMapper<codim, TGV>::contains (const Entity& e, int i,
                                                 int cc, int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV>
template <class EntityType>
bool ActiveInactiveMapper<codim, TGV>::contains (const EntityType& e,
                                                 int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::size () const
{
  return gv.size (codim);
}


template <int codim, class TGV>
void ActiveInactiveMapper<codim, TGV>::update (const EntityPointerVector& active,
                                               const EntityPointerVector& inactive,
                                               const EntityPointerVector& others)
{
  delete indices; indices = new int[gv.size (codim)];
  int cnt = 0;
  for (auto& x : others)   indices[iset.index(x)] = cnt++;
  for (auto& x : inactive) indices[iset.index(x)] = cnt++;
  for (auto& x : active)   indices[iset.index(x)] = cnt++;
  
  /*
  gapEntities.clear();
  total = gv.size (codim);
  delete indices;
  indices = new int[total];
  int cnt = 0;

  for (auto it = gv.template begin<codim>(); it != gv.template end<codim>(); ++it) {
    if (gap.isSupported (it->geometry())) {
      gapEntities << *it;
    } else {
      indices[iset.index(*it)] = cnt++;
    }
  }
  
  if (gapSupported == 0)
    DUNE_THROW (Exception, "functor had empty support");

  activeEntities.clear();
  int active = 0;
  int cnt2 = -1;
  for (auto& x : gapEntities) {
    if (contact.isSupported (x.geometry())) {
      activeEntities << *it;
      indices[iset.index(x)] = cnt + gapEntities.size() + cnt2--;
    } else {
      indices[iset.index(x)] = cnt + active++;
    }
  }
   */
}


#endif /* defined (SIGNORINI_ACTIVEINACTIVEMAPPER) */
