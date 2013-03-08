/******************************************************************************
 * functorsupportmapper.hpp                                                   *
 ******************************************************************************/

#ifndef SIGNORINI_FUNCTORSUPPORTMAPPER
#define SIGNORINI_FUNCTORSUPPORTMAPPER

#include <map>
#include <dune/grid/common/mapper.hh>

/*! A mapper to order entities according to problem data.
 
 The entities in the GridView will be sorted according to the support of the
 functor: those which DO lie in the support will be sorted first, then the
 rest.
 

 
 FIXME!!!!!!!!!!!!!!!!!!!!!
 We use a dirty and very ugly hack which relies on Mapper's bad interface: the
 indices returned by Mapper are platform dependent *signed* ints, so we can use
 the sign bit to store whether a particular index is for an entity in the support
 or not.
 
 This leaves enough room on my 64 bit machine and compiler, but looks like very
 bad practice indeed. However it's fast and easy.

 Template parameters:
    TGV: TGridView
    TFN: TFunctor
 */
template <int codim, class TGV, class TFN>
class FunctorSupportMapper
  : public Mapper<typename TGV::Grid, FunctorSupportMapper<codim, TGV, TFN> >
{
  typedef typename TGV::template Codim<codim>::Entity Entity;

  const TFN& func;
  const TGV&   gv;
  const typename TGV::IndexSet& iset;

public:
  FunctorSupportMapper (const TGV& _gv, const TFN& _func);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (const Entity& e, int i, unsigned int cc) const;

  bool contains (const Entity& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;

  int size () const;
  void update ();

private:
  int supp;  //<! Number of elements in the support of the functor
  int ntot;  //<! Total number of elements
  int* indices;
};


template <int codim, class TGV, class TFN>
FunctorSupportMapper<codim, TGV, TFN>::FunctorSupportMapper (const TGV& _gv,
                                                             const TFN& _func)
: gv (_gv), func (_func), iset (_gv.indexSet()), indices (NULL)
{
  update ();
}


template <int codim, class TGV, class TFN>
template<class EntityType>
int FunctorSupportMapper<codim, TGV, TFN>::map (const EntityType& e) const
{
  int idx = indices[iset.index (e)];
  return (idx < 0) ? (-1*idx + supp-1) : idx;
}


template <int codim, class TGV, class TFN>
int FunctorSupportMapper<codim, TGV, TFN>::map (const Entity& e, int i,
                                                unsigned int cc) const
{
  int idx = indices[iset.index (e, i, cc)];
  return (idx < 0) ? (-1*idx + supp-1) : idx;
}


template <int codim, class TGV, class TFN>
bool FunctorSupportMapper<codim, TGV, TFN>::contains (const Entity& e, int i,
                                                      int cc, int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV, class TFN>
template <class EntityType>
bool FunctorSupportMapper<codim, TGV, TFN>::contains (const EntityType& e,
                                                      int& result) const
{
  DUNE_THROW (Exception, "not implemented");
}


template <int codim, class TGV, class TFN>
int FunctorSupportMapper<codim, TGV, TFN>::size () const
{
  return supp;
}


template <int codim, class TGV, class TFN>
void FunctorSupportMapper<codim, TGV, TFN>::update ()
{
  supp = 0;
  ntot = gv.size (codim);
  delete indices;
  indices = new int[ntot];
  int cnt = -1;
  for (auto it = gv.template begin<codim>(); it != gv.template end<codim>(); ++it) {
    if (func.isSupported (it->geometry())) indices[iset.index(*it)] = supp++;
    else                                   indices[iset.index(*it)] = cnt--;
  }
  
  if (supp == 0)
    DUNE_THROW (Exception, "functor had empty support");
}


#endif /* defined (SIGNORINI_FUNCTORSUPPORTMAPPER) */
