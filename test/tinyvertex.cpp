#include "../src/tinyvertex.h"
#include "../src/typetraits.h"

#include "test.h"

using namespace MolDB;
using namespace MolDB::TypeTraits;

template<typename T1, typename T2, typename T3>
struct TinyGraph
{
  enum {
    VertexEdgesUsage
  };

  typedef EmptyType VertexType;
  typedef EmptyType EdgeType;
  typedef unsigned char VertexEdgeSemanticsType;
  typedef unsigned short IndexType;
};

typedef TinyVertex<TinyGraph<EmptyType, EmptyType, EmptyType> > TinyVertexType;

int main()
{
  TinyVector<TinyVertexType*, Int2Type<0> > vertices(10);
  TinyVertexType *vertex = 0;

  MOLDB_COMPARE(sizeof(TinyVertexType), 16);
}
