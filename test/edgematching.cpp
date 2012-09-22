#include <memory>

#include "../src/graph.h"
#include "../src/isomorphism.h"

#include "test.h"

using namespace MolDB;

typedef Graph<TypeTraits::EmptyType, unsigned char> GraphType;
typedef GraphType::VertexType VertexType;
typedef GraphType::EdgeType EdgeType;

void testVerticesMatch()
{
  std::auto_ptr<GraphType> query1(GraphType::fromString("{[{6}],[]}"));
  std::auto_ptr<GraphType> queried1(GraphType::fromString("{[{6}],[]}"));
  MOLDB_COMPARE(Impl::Isomorphism::verticesMatch(query1->vertices()[0], queried1->vertices()[0]), true);

  std::auto_ptr<GraphType> query2(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
  std::auto_ptr<GraphType> queried2(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
  MOLDB_COMPARE(Impl::Isomorphism::verticesMatch(query2->vertices()[0], queried2->vertices()[0]), true);

  std::auto_ptr<GraphType> query3(GraphType::fromString("{[{6}],[{1,0,[6,6]},{1,0,[6]}]}"));
  std::auto_ptr<GraphType> queried3(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
  MOLDB_COMPARE(Impl::Isomorphism::verticesMatch(query3->vertices()[0], queried3->vertices()[0]), false);
}

void testEdgesPathMatches()
{
  std::auto_ptr<GraphType> query1(GraphType::fromString("{[{6}],[{1,0,[6,6,8,6]}]}"));
  std::auto_ptr<GraphType> queried1(GraphType::fromString("{[{6}],[{1,0,[6,6,8,6]}]}"));
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query1->edges()[0]->path(), queried1->edges()[0]->path(), 0, true), 1);
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query1->edges()[0]->path(), queried1->edges()[0]->path(), 1, true), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query1->edges()[0]->path(), queried1->edges()[0]->path(), 0, false), 0);
  
  std::auto_ptr<GraphType> queried2(GraphType::fromString("{[{6}],[{1,0,[6,6,6,8,6]}]}"));
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query1->edges()[0]->path(), queried2->edges()[0]->path(), 1, true), 1);
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query1->edges()[0]->path(), queried2->edges()[0]->path(), 1, false), 0);
  
  std::auto_ptr<GraphType> query2(GraphType::fromString("{[{6}],[{1,0,[6,6,6,8,6,6,6]}]}"));
  MOLDB_COMPARE(Impl::Isomorphism::edgesPathMatches(query2->edges()[0]->path(), queried2->edges()[0]->path(), 0, true), -1);
}

void testEdgesMatch()
{
  std::auto_ptr<GraphType> g(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,8,6]},{2,1,[6,8,6,6]},{1,0,[6,6,8,6]}]}"));
  VertexType *v1 = g->vertices()[0];
  VertexType *v2 = g->vertices()[1];
  EdgeType *e1 = g->edges()[0];
  EdgeType *e2 = g->edges()[1];
  EdgeType *e3 = g->edges()[2];

  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e1, e1, 0), 1);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e1, e2, 0), 2);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e2, e1, 0), 2);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e2, e2, 0), 1);

  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v2, e1, e1, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v2, e1, e2, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v2, e2, e1, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v2, e2, e2, 0), 0);

  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v1, e1, e1, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v1, e1, e2, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v1, e2, e1, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v1, e2, e2, 0), 0);

  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v2, e1, e1, 0), 1);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v2, e1, e2, 0), 2);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v2, e2, e1, 0), 2);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v2, v2, e2, e2, 0), 1);
  
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e3, e3, 0), 1);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e1, e3, 0), 0);
  MOLDB_COMPARE(Impl::Isomorphism::edgesMatch(v1, v1, e2, e3, 0), 0);
}

void testCheckUntargetedEdges()
{
  typedef Impl::Isomorphism::State<GraphType> StateType;

  ////////////////////////////////////////////////////////////////////
  //
  // Test matching untargeted query edges
  //
  // figure edgematching1 from paper
  //
  ////////////////////////////////////////////////////////////////////

  // no untargeted query edges, should always match
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[]}]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[1,2,[]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  
  // case 1
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,0,[6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 2
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 3
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 4
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6]},{2,0,[6]},{2,0,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 5a (case 5 from figure)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6]},{1,0,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,1,[6,6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 5b (extra)
  {
    // counter example for case 5a where the two path overlap
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6]},{1,0,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,1,[6,6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 5c (extra)
  {
    // same as 5a but with 2 vertices to cross
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6]},{1,0,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{1,3,[6]},{2,3,[6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 5d (extra)
  {
    // counter example for case 5c where the two path overlap after crossing the vertices
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6,6]},{1,0,[6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{1,3,[6]},{2,3,[6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 6
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6]},{1,0,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]},{1,2,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 7
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6]},{1,0,[6,6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{1,3,[6]},{2,3,[6,6,6]},{2,0,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6]},{1,0,[6,6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{1,3,[6]},{2,3,[6,6,6]},{3,0,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }

  //
  // Extra cases
  //
  // case 8
  {
    // used vertices can't be crossed
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.usedQueriedVertices[1] = 1;
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 9
  {
    // used vertices can't be crossed
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.usedQueriedVertices[2] = 1;
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 10
  {
    // counter example to 8 & 9
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6}],[{1,0,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.usedQueriedVertices[2] = 1;
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 11
  /*
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[6]},{3,0,[7]},{4,0,[8]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[6]},{3,0,[8]},{4,0,[7]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[7]},{3,0,[6]},{4,0,[8]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[7]},{3,0,[8]},{4,0,[6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[8]},{3,0,[6]},{4,0,[7]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[]},{1,3,[]},{1,4,[]},{2,0,[8]},{3,0,[7]},{4,0,[6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    MOLDB_COMPARE(Impl::Isomorphism::checkUntargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  */
}
 
void testCheckTargetedEdges()
{
  typedef Impl::Isomorphism::State<GraphType> StateType;


  ////////////////////////////////////////////////////////////////////
  //
  // Test matching targeted query edges
  //
  // figure edgematching2 from paper
  //
  ////////////////////////////////////////////////////////////////////

  std::cout << "111111111111" << std::endl;
  // case 1
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 1b (wrong mapping)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6,6]},{2,3,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 0);
  }
  // case 2
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 0);
  }
  std::cout << "222222222222" << std::endl;
  // case 3
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 0);
   }
  // case 4a (case 4 from figure)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 2);
    MOLDB_COMPARE(state.usedQueriedVertices[2], 2);
  }
  // case 4b (wrong mapping)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[2], 0);
  }
  // case 4c (4a with longer path)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[6]},{2,3,[6]},{3,4,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[3];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
    MOLDB_COMPARE(state.usedQueriedVertices[0], 0);
    MOLDB_COMPARE(state.usedQueriedVertices[1], 2);
    MOLDB_COMPARE(state.usedQueriedVertices[2], 2);
    MOLDB_COMPARE(state.usedQueriedVertices[3], 2);
  }
  std::cout << "333333333333" << std::endl;
  // case 5a (case 5 from figure)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 5b (variation on 5 with query path to long)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 5c (used queried vertices can't be crossed)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    state.usedQueriedVertices[1] = 1;
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }

  ////////////////////////////////////////////////////////////////////
  //
  // Test matching targeted query edges
  //
  // figure edgematching3 from paper
  //
  ////////////////////////////////////////////////////////////////////

  std::cout << "444444444444" << std::endl;
  // case 1a (case 1 from figure)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]},{1,2,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]},{1,2,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    state.usedQueriedVertices[0] = 1;
    state.usedQueriedVertices[1] = 2;
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  // case 1b (queried has only one edge)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]},{1,2,[6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6}],[{1,2,[6,6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[1];
    state.usedQueriedVertices[0] = 1;
    state.usedQueriedVertices[1] = 2;
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
  // case 2a (case 1 from figure)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]},{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6},{6}],[{1,2,[6]},{1,3,[6]},{2,4,[6]},{3,4,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[3];
    state.usedQueriedVertices[0] = 1;
    state.usedQueriedVertices[3] = 2;
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 1);
  }
  std::cout << "555555555555" << std::endl;
  // case 2b (queried has only one edge)
  {
    std::auto_ptr<GraphType> query(GraphType::fromString("{[{6},{6}],[{1,2,[6,6,6]},{1,2,[6,6,6]}]}"));
    std::auto_ptr<GraphType> queried(GraphType::fromString("{[{6},{6},{6}],[{1,2,[6]},{2,3,[6]}]}"));
    StateType state(query.get(), queried.get());
    state.queryPath.push_back(query->vertices()[1]);
    state.queriedPath.push_back(queried->vertices()[1]);
    state.queryPath.push_back(query->vertices()[0]);
    state.queriedPath.push_back(queried->vertices()[0]);
    state.mapping[0] = queried->vertices()[0];
    state.mapping[1] = queried->vertices()[2];
    state.usedQueriedVertices[0] = 1;
    state.usedQueriedVertices[3] = 2;
    MOLDB_COMPARE(Impl::Isomorphism::checkTargetedEdges(state, query->vertices()[0], queried->vertices()[0]), 0);
  }
 










}
 
int main()
{
  std::cout << "testVerticesMatch" << std::endl;
  testVerticesMatch();
  std::cout << "testEdgesMatch" << std::endl;
  testEdgesMatch();
  std::cout << "testUntargetedEdges" << std::endl;
  testCheckUntargetedEdges();
  std::cout << "testTargetedEdges" << std::endl;
  testCheckTargetedEdges();
}
