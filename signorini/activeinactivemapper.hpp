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
  typedef typename TGV::template Codim<0>::Entity Element;
  typedef typename TGV::template Codim<codim>::Entity Entity;
  typedef std::vector<const Entity*> EntityVector;
  
  typedef typename TGV::Grid::GlobalIdSet      GlobalIdSet;
  typedef typename TGV::Grid::GlobalIdSet::IdType   IdType;
  typedef std::set<IdType>                           IdSet;
 
  const TGV& gv;
  const GlobalIdSet& gids;
  const typename TGV::IndexSet& iset;
  
  std::map<IdType, int> indices;
  int offsetInner;
  int offsetActive;
public:
  ActiveInactiveMapper (const TGV& _gv,
                        const IdSet& active,
                        const IdSet& inactive,
                        const IdSet& other);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (IdType id) const;
  int map (const Element& e, int i, unsigned int cc) const;
  int map (Element& e, int i, unsigned int cc) const;
  int mapInBoundary (IdType id) const;
  int mapInBoundary (const Element& e, int i, unsigned int cc) const;
  int mapInBoundary (Element& e, int i, unsigned int cc) const;
  int mapInActive (IdType id) const;
  int mapInActive (const Element& e, int i, unsigned int cc) const;
  int mapInActive (Element& e, int i, unsigned int cc) const;
  
  bool contains (const Element& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;
  
  int size () const;
  void update (const IdSet& active,
               const IdSet& inactive,
               const IdSet& other);
};


template <int codim, class TGV>
ActiveInactiveMapper<codim, TGV>::ActiveInactiveMapper (const TGV& _gv,
                                                        const IdSet& active,
                                                        const IdSet& inactive,
                                                        const IdSet& other)
: gv (_gv), gids(_gv.grid().globalIdSet()), iset (_gv.indexSet()),
  offsetInner(0), offsetActive(0)
{
  update (active, inactive, other);
}


template <int codim, class TGV>
template<class EntityType>
int ActiveInactiveMapper<codim, TGV>::map (const EntityType& e) const
{
  return indices.at (gids.id(e));
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (const IdType id) const
{
  return indices.at (id);
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (const Element& e, int i,
                                           unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc));
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (Element& e, int i,
                                           unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc));
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInBoundary (const IdType id) const
{
  return indices.at (id)-offsetInner;
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInBoundary (const Element& e, int i,
                                                     unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc))-offsetInner;
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInBoundary (Element& e, int i,
                                                     unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc))-offsetInner;
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInActive (const IdType id) const
{
  return indices.at (id)-offsetActive;
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInActive (const Element& e, int i,
                                                   unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc))-offsetActive;
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::mapInActive (Element& e, int i,
                                                     unsigned int cc) const
{
  return indices.at (gids.subId (e, i, cc))-offsetActive;
}

template <int codim, class TGV>
bool ActiveInactiveMapper<codim, TGV>::contains (const Element& e, int i,
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
void ActiveInactiveMapper<codim, TGV>::update (const IdSet& active,
                                               const IdSet& inactive,
                                               const IdSet& other)
{
  indices.clear();
  int cnt = 0;
  for (auto x : other)   indices[x] = cnt++;
  offsetInner = cnt;
  for (auto x : inactive) indices[x] = cnt++;
  offsetActive = cnt;
  for (auto x : active)   indices[x] = cnt++;
}


#endif /* defined (SIGNORINI_ACTIVEINACTIVEMAPPER) */
