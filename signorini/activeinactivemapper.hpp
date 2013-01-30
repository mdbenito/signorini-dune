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

public:
  ActiveInactiveMapper (const TGV& _gv,
                        const IdSet& active,
                        const IdSet& inactive,
                        const IdSet& others);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (IdType id) const;
  int map (const Element& e, int i, unsigned int cc) const;
  int map (Element& e, int i, unsigned int cc) const;
  
  bool contains (const Element& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;
  
  int size () const;
  void update (const IdSet& active,
               const IdSet& inactive,
               const IdSet& others);
};


template <int codim, class TGV>
ActiveInactiveMapper<codim, TGV>::ActiveInactiveMapper (const TGV& _gv,
                                                        const IdSet& active,
                                                        const IdSet& inactive,
                                                        const IdSet& others)
: gv (_gv), gids(_gv.grid().globalIdSet()), iset (_gv.indexSet())
{
  update (active, inactive, others);
}


template <int codim, class TGV>
template<class EntityType>
int ActiveInactiveMapper<codim, TGV>::map (const EntityType& e) const
{
  return indices.at(gids.id(e));
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (const IdType id) const
{
  return indices.at(id);
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (const Element& e, int i,
                                           unsigned int cc) const
{
  return indices.at(gids.subId (e, i, cc));
}

template <int codim, class TGV>
int ActiveInactiveMapper<codim, TGV>::map (Element& e, int i,
                                           unsigned int cc) const
{
  return indices.at(gids.subId (e, i, cc));
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
                                               const IdSet& others)
{
  indices.clear();
  int cnt = 0;
  for (auto x : others)   indices[x] = cnt++;
  for (auto x : inactive) indices[x] = cnt++;
  for (auto x : active)   indices[x] = cnt++;
}


#endif /* defined (SIGNORINI_ACTIVEINACTIVEMAPPER) */
