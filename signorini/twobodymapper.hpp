/******************************************************************************
 * twobodymapper.hpp                                                          *
 ******************************************************************************/

#ifndef TWOBODY_ACTIVEINACTIVEMAPPER
#define TWOBODY_ACTIVEINACTIVEMAPPER

#include <map>
#include <dune/grid/common/mapper.hh>

  //template <int codim, class TGV> class TwoToOneBodyMapper;

/*! A mapper to order entities according to their status as active or inactive.
 
 For each body the constructor is given three sets of entities representing those
 which are active "A", inactive "I" and the rest "N". The final ordering returned
 by the mapper will be:
 
 N1_1 ... N1_r1 N2_1 ... N2_r2 I1_1 ... I1_m1 I2_1 ... I2_m2 A1_1 ... A1_n1 A2_1 ... A2_n2
 
 Template parameters:
 TGV: TGridView
 */
template <int codim, class TGV>
class TwoBodyMapper
: public Mapper<typename TGV::Grid, TwoBodyMapper<codim, TGV> >
{
  typedef typename TGV::template Codim<0>::Entity Element;
  typedef typename TGV::template Codim<codim>::Entity Entity;
  typedef std::vector<const Entity*> EntityVector;
  
  typedef typename TGV::Grid::GlobalIdSet      GlobalIdSet;
  typedef typename TGV::Grid::GlobalIdSet::IdType   IdType;
  typedef std::set<IdType>                           IdSet;
  
  TwoRefs<TGV> gv;
  TwoRefs<GlobalIdSet> gids;
  
  std::map<IdType, int> indices[2];
  int offsetBody[2];
  int offsetInner[2];
  int offsetActive[2];
  int sizeOther[2];

  inline int calcOffset (int body, int idx) const;
public:
  TwoBodyMapper (const TwoRefs<TGV>& _gv,
                 const IdSet active[],
                 const IdSet inactive[],
                 const IdSet other[]);
  
  template<class EntityType> int map (int body, const EntityType& e) const;
  int map (int body, IdType id) const;
  int map (int body, const Element& e, int i, unsigned int cc) const;
  template<class EntityType> int mapInBody (int body, const EntityType& e) const;
  int mapInBody (int body, IdType id) const;
  int mapInBody (int body, const Element& e, int i, unsigned int cc) const;
  int mapInBoundary (int body, IdType id) const;
  int mapInBoundary (int body, const Element& e, int i, unsigned int cc) const;
  int mapInActive (int body, IdType id) const;
  int mapInActive (int body, const Element& e, int i, unsigned int cc) const;
  
  bool contains (int body, const Element& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (int body, const EntityType& e, int& result) const;
  
  int size (int body) const;
  void update (const IdSet* active, const IdSet* inactive, const IdSet* other);
};


template <int codim, class TGV>
TwoBodyMapper<codim, TGV>::TwoBodyMapper (const TwoRefs<TGV>& _gv,
                                          const IdSet* active,
                                          const IdSet* inactive,
                                          const IdSet* other)
: gv (_gv),
  gids (TwoRefs<GlobalIdSet> (_gv[0].grid().globalIdSet(), _gv[1].grid().globalIdSet())),
  offsetBody(), offsetInner(), offsetActive(), sizeOther()
{
  update (active, inactive, other);
}


template <int codim, class TGV>
template<class EntityType>
int TwoBodyMapper<codim, TGV>::map (int body, const EntityType& e) const
{
  return indices[body].at (gids[body].id(e));
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::map (int body, const IdType id) const
{
  return indices[body].at (id);
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::map (int body, const Element& e, int i,
                                    unsigned int cc) const
{
  return indices[body].at (gids[body].subId (e, i, cc));
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInBoundary (int body, const IdType id) const
{
  return indices[body].at (id) - offsetInner[body];
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInBoundary (int body, const Element& e, int i,
                                              unsigned int cc) const
{
  return indices[body].at (gids[body].subId (e, i, cc)) - offsetInner[body];
}

template <int codim, class TGV>
template<class EntityType>
int TwoBodyMapper<codim, TGV>::mapInBody (int body, const EntityType& e) const
{
  return calcOffset (body, indices[body].at (gids[body].id(e)));
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInBody (int body, const IdType id) const
{
  return calcOffset (body, indices[body].at (id));
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInBody (int body, const Element& e, int i,
                                              unsigned int cc) const
{
  return calcOffset (body, indices[body].at (gids[body].subId (e, i, cc)));
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInActive (int body, const IdType id) const
{
  return indices[body].at (id) - offsetActive[body];
}

template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::mapInActive (int body, const Element& e, int i,
                                            unsigned int cc) const
{
  return indices[body].at (gids[body].subId (e, i, cc)) - offsetActive[body];
}

template <int codim, class TGV>
bool TwoBodyMapper<codim, TGV>::contains (int body, const Element& e, int i,
                                          int cc, int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV>
template <class EntityType>
bool TwoBodyMapper<codim, TGV>::contains (int body, const EntityType& e,
                                          int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::size (int body) const
{
  return gv[body].size (codim);
}

#define DOBOTH(x) for (int x=0; x < 2; ++x)

template <int codim, class TGV>
void TwoBodyMapper<codim, TGV>::update (const IdSet* active,
                                        const IdSet* inactive,
                                        const IdSet* other)
{
  int cnt = 0;
  DOBOTH (body) {
    indices[body].clear();
    sizeOther[body] = other[body].size();
    offsetBody[body] = cnt;
    for (auto x : other[body])
      indices[body][x] = cnt++;
  }

  DOBOTH (body) {
    offsetInner[body] = cnt;
    for (auto x : inactive[body])
      indices[body][x] = cnt++;
    offsetActive[body] = cnt;
    
    for (auto x : active[body])
      indices[body][x] = cnt++;
  }
}
template <int codim, class TGV>
int TwoBodyMapper<codim, TGV>::calcOffset (int body, int idx) const
{
    //cout << "idx= " << idx;
  if (idx >= offsetInner[body])
    idx = idx + sizeOther[body] - offsetInner[body];
    //cout << " --> " << idx - offsetBody[body] << " (" << body << ")" << LF;
  return idx - offsetBody[body];
}

/* Quick hack to return to the original mapper interface, for use in the
 Postprocessor. Maybe I should keep the references to this inside the 
 TwoBodyMapper and return them with operator[] there.
 */
template <int codim, class TGV>
class TwoToOneBodyMapper
: public Mapper<typename TGV::Grid, TwoToOneBodyMapper<codim, TGV> >
{
  typedef typename TGV::template Codim<0>::Entity Element;
  typedef typename TGV::template Codim<codim>::Entity Entity;
  typedef std::vector<const Entity*> EntityVector;
  
  typedef typename TGV::Grid::GlobalIdSet      GlobalIdSet;
  typedef typename TGV::Grid::GlobalIdSet::IdType   IdType;
  typedef std::set<IdType>                           IdSet;

  const int body;
  const TwoBodyMapper<codim, TGV>& mapper;

public:
  TwoToOneBodyMapper (int _body, const TwoBodyMapper<codim, TGV>& _mapper)
  : body(_body), mapper(_mapper) {}
  
  template<class EntityType> int map (const EntityType& e) const {
    return mapper.mapInBody (body, e);
  }
  int map (IdType id) const {
    return mapper.mapInBody (body, id);
  }
  int map (int body, const Element& e, int i, unsigned int cc) const {
    return mapper.mapInBody (body, e, i, cc);
  }
  bool contains (const Element& e, int i, int cc, int& result) const {
    return mapper.contains (body, e, i, cc, result);
  }
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const {
    return mapper.contains (body, e, result);
  }
  int size () const {
    return mapper.size (body);
  }
};

#endif /* defined (TWOBODY_ACTIVEINACTIVEMAPPER) */
