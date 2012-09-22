#include "../src/tinyedge.h"
#include "../src/typetraits.h"

#include "test.h"

using namespace MolDB;
using namespace MolDB::TypeTraits;

template<typename T1, typename T2, typename T3>
struct TinyGraph
{
  enum {
    EdgePathUsage
  };

  typedef EmptyType VertexType;
  typedef EmptyType EdgeType;
  typedef unsigned char VertexEdgeSemanticsType;
  typedef unsigned short IndexType;
};

typedef TinyEdge<TinyGraph<EmptyType, EmptyType, EmptyType> > TinyEdgeType;

int main()
{
  TinyVector<TinyEdgeType*, Int2Type<0> > vertices(10);
  TinyEdgeType *edge = 0;
  
  MOLDB_COMPARE(sizeof(TinyEdgeType), 16);
}
