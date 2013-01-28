/******************************************************************************
 * vertexboundarymapper.hpp                                                   *
 *                                                                            *
 ******************************************************************************/
#ifndef SIGNORINI_VERTEXBOUNDARYMAPPER
#define SIGNORINI_VERTEXBOUNDARYMAPPER

#include <map>
#include <dune/grid/common/mapper.hh>

/*! A mapper to order vertices according to boundary data.
 
 The vertices in the grid will be sorted according to the support of the 
 functor: those which do NOT lie in the support will be sorted first, then the
 rest.
 
 Template parameters:
    TG: TGrid
    TF: TFunctor
 
 */
template <class TG, class TF>
class VertexBoundayMapper : public Mapper<TG, VertexBoundayMapper<TG, TF> >
{
  typedef typename TG::Traits::template Codim<0>::Entity Entity;

  TG&       grid;
  const TF& func;
  
public:
  VertexBoundayMapper (const TG& _grid, const TF& _func);
  
  template<class EntityType> int map (const EntityType& e) const;
  int map (const Entity& e, int i, unsigned int codim) const;

  bool contains (const Entity& e, int i, int cc, int& result) const;
  template<class EntityType>
  bool contains (const EntityType& e, int& result) const;

  int size () const;
  void update ();
};


template <class TG, class TF>
VertexBoundayMapper<TG, TF>::VertexBoundayMapper (const TG& _grid,
                                                  const TF& _func)
: grid (_grid), func (_func)
{
  
}


template <class TG, class TF>
template<class EntityType>
int VertexBoundayMapper<TG, TF>::map (const EntityType& e) const
{
  
}


template <class TG, class TF>
int VertexBoundayMapper<TG, TF>::map (const Entity& e, int i,
                                      unsigned int codim) const
{
  
}


template <class TG, class TF>
bool VertexBoundayMapper<TG, TF>::contains (const Entity& e, int i, int cc,
                                            int& result) const
{
  
}


template <class TG, class TF>
template<class EntityType>
bool VertexBoundayMapper<TG, TF>::contains (const EntityType& e,
                                            int& result) const
{
  
}


template <class TG, class TF>
int VertexBoundayMapper<TG, TF>::size () const
{
  
}


template <class TG, class TF>
void VertexBoundayMapper<TG, TF>::update ()
{
  
}


#endif /* defined (SIGNORINI_VERTEXBOUNDARYMAPPER) */
