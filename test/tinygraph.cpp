#include "../src/tinygraph.h"

#include "test.h"

using namespace MolDB;

typedef TinyGraph<TypeTraits::EmptyType, unsigned char, unsigned short> TinyGraphType;
typedef TinyGraphType::VertexType TinyVertexType;
typedef TinyGraphType::EdgeType TinyEdgeType;

int main()
{
  MOLDB_COMPARE(sizeof(TinyGraphType), 24);

  TinyGraphType g(10, 10);
  MOLDB_COMPARE(g.numVertices(), 0);
  MOLDB_COMPARE(g.numEdges(), 0);

  std::vector<TinyGraphType::IndexType> edges1;
  edges1.push_back(4);
  edges1.push_back(6);
  TinyVertexType *vertex1 = g.addVertex(6, edges1);
  MOLDB_COMPARE(g.numVertices(), 1);
  MOLDB_COMPARE(g.numEdges(), 0);
  MOLDB_COMPARE(g.vertex(0), vertex1);
  MOLDB_COMPARE(vertex1->edges().size(), 2);
  MOLDB_COMPARE(vertex1->edges()[0], 4);
  MOLDB_COMPARE(vertex1->edges()[1], 6);
 
  std::vector<TinyGraphType::IndexType> edges2;
  edges2.push_back(4);
  edges2.push_back(6);
  edges2.push_back(9);
  TinyVertexType *vertex2 = g.addVertex(8, edges2);
  MOLDB_COMPARE(g.numVertices(), 2);
  MOLDB_COMPARE(g.numEdges(), 0);
  MOLDB_COMPARE(g.vertex(1), vertex2);
  MOLDB_COMPARE(vertex2->edges().size(), 3);
  MOLDB_COMPARE(vertex2->edges()[0], 4);
  MOLDB_COMPARE(vertex2->edges()[1], 6);
  MOLDB_COMPARE(vertex2->edges()[2], 9);

  std::vector<TinyGraphType::VertexEdgeSemanticsType> path;
  path.push_back(6);
  path.push_back(7);
  path.push_back(8);
  TinyEdgeType *edge = g.addEdge(vertex1, vertex2, path);
  MOLDB_COMPARE(g.numVertices(), 2);
  MOLDB_COMPARE(g.numEdges(), 1);
  MOLDB_COMPARE(g.edge(0), edge);
  MOLDB_COMPARE(edge->source(), 1);
  MOLDB_COMPARE(edge->target(), 2);
  MOLDB_COMPARE(edge->path().size(), 3);
  MOLDB_COMPARE(edge->path()[0], 6);
  MOLDB_COMPARE(edge->path()[1], 7);
  MOLDB_COMPARE(edge->path()[2], 8);

 

}
