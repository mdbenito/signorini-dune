/******************************************************************************
 * activeinactivemapper.hpp                                                   *
 ******************************************************************************/

#ifndef SIGNORINI_ACTIVEINACTIVEMAPPER
#define SIGNORINI_ACTIVEINACTIVEMAPPER

#include <map>
#include <dune/grid/common/mapper.hh>

/*! A mapper to order entities according to their status as active or inactive.
 
 
 Using the specified functor "gap" the entities are divided in two sets 
 corresponding to those which belong to the support (set "S") and those who
 don't (set "N")
 
 The set "S" is further divided in active "A" and inactive "I", according to
 an additional functor "contact"
 
 FIXME!!!!!!!!!!!!!!!!!!!!!
 We use a dirty and very ugly hack which relies on Mapper's bad interface: the
 indices returned by Mapper are platform dependent *signed* ints, so we can use
 the sign bit to store whether a particular index is for an entity in the support
 or not.
 
 This leaves enough room on my 64 bit machine and compiler, but looks like very
 bad practice indeed. However it's fast and easy.
 
 Template parameters:
 TGV: TGridView
 TFN: TFunctor   (for the gap function)
 TAC: TActiveFunctor (to determine active entities in the 
 */
template <int codim, class TGV, class TFN, class TAC>
class ActiveInactiveMapper
: public Mapper<typename TGV::Grid, ActiveInactiveMapper<codim, TGV, TFN, TAC> >
{
  typedef typename TGV::template Codim<0>::Entity Entity;
  
  const TFN&     gap;
  const TAC& contact;
  const TGV&      gv;
  const typename TGV::IndexSet& iset;
  
public:
  ActiveInactiveMapper (const TGV& _gv, const TFN& _gap, const TAC& _contact);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (const Entity& e, int i, unsigned int cc) const;
  
  bool contains (const Entity& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;
  
  int size () const;
  void update ();
  
private:
  int gapSupported;  //<! Number of elements in the support of the gap functor
  int       active;
  int        total;  //<! Total number of elements
  int*     indices;
};


template <int codim, class TGV, class TFN, class TAC>
ActiveInactiveMapper<codim, TGV, TFN, TAC>::ActiveInactiveMapper (const TGV& _gv,
                                                                  const TFN& _gap,
                                                                  const TAC& _contact)
: gv (_gv), gap (_gap), contact (_contact), iset (_gv.indexSet()), indices (NULL)
{
  update ();
}


template <int codim, class TGV, class TFN, class TAC>
template<class EntityType>
int ActiveInactiveMapper<codim, TGV, TFN, TAC>::map (const EntityType& e) const
{
  int idx = indices[iset.index (e)];
  return (idx < 0) ? (-1*idx + gapSupported-1) : idx;
}


template <int codim, class TGV, class TFN, class TAC>
int ActiveInactiveMapper<codim, TGV, TFN, TAC>::map (const Entity& e, int i,
                                                     unsigned int cc) const
{
  int idx = indices[iset.index (e, i, cc)];
  return (idx < 0) ? (-1*idx + gapSupported-1) : idx;
}


template <int codim, class TGV, class TFN, class TAC>
bool ActiveInactiveMapper<codim, TGV, TFN, TAC>::contains (const Entity& e, int i,
                                                           int cc, int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV, class TFN, class TAC>
template <class EntityType>
bool ActiveInactiveMapper<codim, TGV, TFN, TAC>::contains (const EntityType& e,
                                                           int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV, class TFN, class TAC>
int ActiveInactiveMapper<codim, TGV, TFN, TAC>::size () const
{
  return total;
}


template <int codim, class TGV, class TFN, class TAC>
void ActiveInactiveMapper<codim, TGV, TFN, TAC>::update ()
{
  gapSupported = 0;

  total = gv.size (codim);
  delete indices;
  indices = new int[total];
  int cnt = -1;
  
  std::vector<typename TGV::EntityPointer> inGap (total/2);  // z.B...

  for (auto it = gv.template begin<codim>(); it != gv.template end<codim>(); ++it) {
    if (gap.isSupported (it->geometry())) {
      inGap << *it;
      gapSupported++;
    } else {
      indices[iset.index(*it)] = cnt--;
    }
  }
  
  if (gapSupported == 0)
    DUNE_THROW (Exception, "functor had empty support");
  active = 0;
  cnt = 1;
  for (auto& x : inGap) {
    if (contact.isSupported(x.geometry()))
      indices[iset.index(x)] = active++;
    else
      indices[iset.index(x)] = gapSupported-(cnt--);
  }
}


#endif /* defined (SIGNORINI_ACTIVEINACTIVEMAPPER) */
