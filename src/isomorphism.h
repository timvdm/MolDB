/**
 * Copyright (c) 2012, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MOLDB_ISOMORPHISM_H
#define MOLDB_ISOMORPHISM_H

#include <vector>
#include <iostream>

#include "semantics.h"

#define DEBUG_MOLDB_ISOMORPHISM_H 0

namespace MolDB {
  
  template<typename T>
  void print_vector(const std::string &label, const std::vector<T> &v)
  {
    std::cout << label << ": ";
    for (std::size_t i = 0; i < v.size(); ++i)
      //if (v[i] < 10)
      //  std::cout << " " << v[i] << " ";
      //else
        std::cout << v[i] << " ";

    std::cout << std::endl;
  }

  template<typename C, typename T>
  inline std::size_t indexOf(const C &c, T v)
  {
    return std::find(c.begin(), c.end(), v) - c.begin();
  }


  namespace Impl {

    namespace Isomorphism {
      
      static const char *red    = "\033[1;31m";
      static const char *green  = "\033[1;32m";
      static const char *yellow = "\033[1;33m";
      static const char *blue   = "\033[1;34m";
      static const char *normal = "\033[0m";

      struct UsedEdge
      {
        /**
         * The various use cases for used edges.
         */
        enum Use {
          NotUsed,
          Forward,
          Backward,
          Both 
        };

        UsedEdge() : depth(0), use(NotUsed), count(0)
        {
        }

        UsedEdge(unsigned int depth_, Use use_, unsigned int count_)
            : depth(depth_), use(use_), count(count_)
        {
        }

        /**
         * The depth at which this edge is used. This data member is used for
         * backtracking.
         */
        unsigned int depth;
        /**
         * The use of the edge
         */
        Use use;
        /**
         * The number of path elements used.
         */
        unsigned int count;
      };

      struct CompareUsedEdgeDepth
      {
        CompareUsedEdgeDepth(unsigned int depth_) : depth(depth_)
        {
        }

        bool operator()(const UsedEdge &e)
        {
          return e.depth == depth;
        }
        
        unsigned int depth;
      };


      std::ostream& operator<<(std::ostream &os, const UsedEdge &e)
      {
        os << "[ depth = " << e.depth << ", use = ";
        switch (e.use) {
          case UsedEdge::NotUsed:
            os << "NotUsed, ";
            break;
          case UsedEdge::Forward:
            os << "Forward, ";
            break;
          case UsedEdge::Backward:
            os << "Backward, ";
            break;
          case UsedEdge::Both:
            os << "Both, ";
            break;
        }
        os << "count = " << e.count << " ]";
      }

      template<typename _GraphType>
      struct SharedState
      {
        typedef _GraphType GraphType;
        typedef typename GraphType::VertexType VertexType;

        SharedState(std::size_t querySize, std::size_t queriedSize)
            : queryMapping(querySize, 0), queriedMapping(queriedSize, 0),
            queryDepths(querySize, 0), queriedDepths(queriedSize, 0)
        {
        }

        std::vector<VertexType*> queryMapping;
        std::vector<VertexType*> queriedMapping;
        std::vector<std::size_t> queryDepths;
        std::vector<std::size_t> queriedDepths;
      };

      /**
       * @class State isomorphism.h
       *
       * The state object. This class holds all information needed during the
       * mapping process.
       */
      template<typename _GraphType>
      struct State {
        /**
         * The graph type.
         */
        typedef _GraphType GraphType;
        /**
         * The vertex type.
         */
        typedef typename GraphType::VertexType VertexType;
        /**
         * The edge type.
         */
        typedef typename GraphType::EdgeType EdgeType;
        /**
         * The candidate pair type.
         */
        typedef std::pair<VertexType*, VertexType*> CandidateType;

        typedef typename GraphType::VerticesType GraphVerticesType;
        typedef typename GraphType::EdgesType GraphEdgesType;
        typedef typename VertexType::EdgesType VertexEdgesType;
        typedef typename EdgeType::PathType EdgePathType;
        typedef typename GraphType::IndexType IndexType;

        /**
         * Constructor.
         */
        State(const GraphType *query_, const GraphType *queried_)
          : mappingSize(0), queryDepthsSize(0), queriedDepthsSize(0),
          query(query_), queried(queried_), lastMapped(0, 0),
          sharedState(new SharedState<GraphType>(query->numVertices(), queried->numVertices())),
          ownedSharedState(true)
        {
          // resize and initialze usedQueriedVertices
          usedQueriedVertices.resize(queried->numVertices(), 0);
          usedQueriedEdges.resize(queried->numEdges());
        }

        State(const State &other) : mappingSize(other.mappingSize), queryDepthsSize(other.queryDepthsSize),
            queriedDepthsSize(other.queriedDepthsSize), query(other.query), queried(other.queried),
            lastMapped(0, 0), sharedState(other.sharedState), ownedSharedState(false)
        {
          // FIXME usedQueried...
        }

        ~State()
        {
          if (ownedSharedState)
            delete sharedState;
        }

        std::size_t mappingSize;
        std::size_t queryDepthsSize;
        std::size_t queriedDepthsSize;
        const GraphType *query;
        const GraphType *queried;
        CandidateType lastMapped;
        SharedState<GraphType> *sharedState;
        bool ownedSharedState;

        /**
         * True if a mapping is found.
         */
        //bool match;
        /**
         * The query graph.
         */
        //const GraphType *query;
        /**
         * The queried graph.
         */
        //const GraphType *queried;
        /**
         * The query path.
         */
        //std::vector<VertexType*> queryPath;
        /**
         * The queried path.
         */
        //std::vector<VertexType*> queriedPath;
        /**
         * The query depths of vertices from the queried graph were added.
         */
        //std::vector<unsigned int> queryDepths;
        /**
         * The queried depths of vertices from the queried graph were added.
         */
        //std::vector<unsigned int> queriedDepths;
        /**
         * The mapping from query to queried. Contains queried vertices,
         * indexed by query vertex index.
         */
        //std::vector<VertexType*> mapping;
        /**
         * The used queried vertices. This vector is indexed by queried vertex
         * index and contains the query edge target that crossed or overlapped
         * with it.
         */
        std::vector<std::size_t> usedQueriedVertices;
        std::vector<UsedEdge> usedQueriedEdges;
      };

      /**
       * Compare two vertices' semantic attributes. This function also compares
       * the number of edges in the query and queried vertices. If the query has
       * more vertices than the queried, it follows that there can be no match.
       *
       * @note The runtime complexity for this function is O(1). This can be
       *       made worse if comparing the vertices semantic attributs has a
       *       runtime worse than O(1).
       *
       * @param query The query vertex.
       * @param queried The queried vertex.
       *
       * @return True if the vertices match, false otherwise.
       */
      template<typename VertexType>
      bool verticesMatch(VertexType *query, VertexType *queried)
      {
        // compare semantic attributes
        if (query->semantics() != queried->semantics())
          return false;
        // compare number of edges
        if (query->numEdges() > queried->numEdges())
          return false;
        int numTargetedQueryEdges = 0;
        for (std::size_t i = 0; i < query->numEdges(); ++i)
          numTargetedQueryEdges++;
        int numTargetedQueriedEdges = 0;
        for (std::size_t i = 0; i < queried->numEdges(); ++i)
          numTargetedQueriedEdges++;
        if (numTargetedQueryEdges > numTargetedQueriedEdges)
          return false;
        // the vertices match
        return true;
      }

      /**
       * Compare the semantic attributes stored in and edge's path. The 
       * attributes can be checked in a forward or backwards order. It is also
       * possible to specify the starting index for the query's path.
       *
       * @note This function should not be used directly since there is a
       *       edgesMatch() wrapper function that calls this function with the
       *       correct @p forward argument.
       *
       * @note The runtime complexity for this function is O(n) where n is the
       *       length of the longest path. The runtime complexity is also
       *       influenced by the complexity of comparing path element which is
       *       not guaranteed to be O(1).
       *
       * @param query The query path.
       * @param queried The queried path.
       * @param queryStartIndex The start index in the query path to start
       *        comparing from.
       * @param forward Specify if the matching should be done starting from the
       *        queried path's 1st element. The default (true) is to match in a
       *        forward direction. If false, matching will be done starting with
       *        queried's last element first.
       *
       * @return 0 if there is no match. If the query paths' length is shorter
       *         than (or equal in length) the queried path and all element
       *         match, 1 is returned if the queried edge was compared forwards
       *         and 2 if compared backwards. If the query paths' length is
       *         longer than queried paths' but all elements up to queried
       *         paths' length match, -1 is returned.
       */
      template<typename PathType>
      int edgesPathMatches(const PathType &query, const PathType &queried,
          std::size_t queryStartIndex, bool forward = true)
      {
        if (queried.empty()) // FIXME (check before...)
          return false;
        typename PathType::const_iterator i = query.begin() + queryStartIndex;
        if (forward) {
          // compare forwards
          typename PathType::const_iterator j = queried.begin();
          for (; i != query.end(); ++i, ++j) {
            if (j == queried.end())
              return -1;
            if (!semanticsMatch(*i, *j))
              return 0;
          }
          
          return 1;
        } else {
          // compare backwards
          typename PathType::const_reverse_iterator j = queried.rbegin();
          for (; i != query.end(); ++i, ++j) {
            if (j == queried.rend())
              return -1;
            if (!semanticsMatch(*i, *j))
              return 0;
          }
          
          return 2;
        }
      }


      template<typename VertexType, typename EdgeType>
      inline bool forwardDirection(VertexType *querySource, VertexType *queriedSource, 
          EdgeType *queryEdge, EdgeType *queriedEdge)
      {
        bool forward;
        if (querySource == queryEdge->source())
          forward = (queriedSource == queriedEdge->source()) ? true : false;
        else
          forward = (queriedSource == queriedEdge->source()) ? false : true;
        return forward;
       }
 

      /**
       * Compare the semantic attributes stored in and edge's path. The 
       * attributes can be checked in a forward or backwards order. It is also
       * possible to specify the starting index for the query's path.
       *
       * @note The runtime complexity for this function is O(n) where n is the
       *       length of the longest path. The runtime complexity is also
       *       influenced by the complexity of comparing path element which is
       *       not guaranteed to be O(1).
       *
       * @param query The query path.
       * @param queried The queried path.
       * @param queryStartIndex The start index in the query path to start
       *        comparing from.
       * @param forward Specify if the matching should be done starting from the
       *        queried path's 1st element. The default (true) is to match in a
       *        forward direction. If false, matching will be done starting with
       *        queried's last element first.
       *
       * @return 0 if there is no match. If the query paths' length is shorter
       *         than (or equal in length) the queried path and all element
       *         match, 1 is returned if the queried edge was compared forwards
       *         and 2 if compared backwards. If the query paths' length is longer
       *         than queried paths' but all elements up to queried paths'
       *         length match, -1 is returned.
       */
      template<typename VertexType, typename EdgeType>
      inline int edgesMatch(VertexType *querySource, VertexType *queriedSource, 
          EdgeType *queryEdge, EdgeType *queriedEdge, std::size_t queryStartIndex = 0,
          bool forceBackward = false)
      {
        // if query edge has a target, queried edge must also have a target
        if (queryEdge->target() && !queriedEdge->target())
          return false;
        if (queryEdge->path().empty() && queriedEdge->path().empty())
          return true;
    
        if (forceBackward)
          return edgesPathMatches(queryEdge->path(), queriedEdge->path(), queryStartIndex, false);

        if (queriedEdge->source() == queriedEdge->target()) {
          if (edgesPathMatches(queryEdge->path(), queriedEdge->path(), queryStartIndex, true))
            return 1;
          if (edgesPathMatches(queryEdge->path(), queriedEdge->path(), queryStartIndex, false))
            return 2;
          return 0;
        }

        // determine the direction
        bool forward = forwardDirection(querySource, queriedSource, queryEdge, queriedEdge);
            
        return edgesPathMatches(queryEdge->path(), queriedEdge->path(), queryStartIndex, forward);
      }

                  
      template<bool Targeted, typename StateType>
      bool edgesMatchDFS(StateType &state, 
          typename StateType::EdgeType *queryEdge, 
          typename StateType::VertexType *queryVertex, 
          typename StateType::VertexType *queriedVertex,
          std::size_t queryPathStartIndex)
      {
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::VertexType VertexType;
        typedef typename StateType::VertexEdgesType VertexEdgesType;

        const VertexEdgesType &queriedEdges = queriedVertex->edges();

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << "edgesMatchDFS(queriedVertex = " << queriedVertex->index() << ", index = " << queryPathStartIndex << ")" << std::endl;

        for (std::size_t i = 0; i < queriedVertex->numEdges(); ++i) {
          EdgeType *queriedEdge = state.queried->edge(queriedEdges[i]);
          assert(queriedEdge);


          //print_vector("usedEdges", state.usedQueriedEdges);
          if (!Targeted) {
            // queriedTarget can not be in the path
            //if (std::find(state.queriedPath.begin(), state.queriedPath.end(), queriedEdge->other(queriedVertex)) != state.queriedPath.end())
            //  continue;


            std::size_t queriedEdgeIndex = queriedEdge->index();
            if (state.usedQueriedEdges[queriedEdgeIndex].use != UsedEdge::NotUsed) {
              if (queriedEdge->source() == queriedEdge->target()) {
                if (state.usedQueriedEdges[queriedEdgeIndex].use == UsedEdge::Both)
                  continue;
              } else {
                enum UsedEdge::Use use = forwardDirection(queryVertex, queriedVertex, queryEdge, queriedEdge) ? UsedEdge::Forward : UsedEdge::Backward;
                if (state.usedQueriedEdges[queriedEdgeIndex].use == use)
                  continue;
                if (queriedEdge->pathSize() - state.usedQueriedEdges[queriedEdgeIndex].count < queryEdge->pathSize() - queryPathStartIndex)
                  continue;
              }
            }
          } else {
            // queriedTarget can not be in the path
            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), queriedEdge->other(queriedVertex)) != state.queriedPath.end())
              continue;

            if (queryEdge->pathSize() - queryPathStartIndex < queriedEdge->pathSize())
              continue;
          }
 
          int result = edgesMatch(queryVertex, queriedVertex, queryEdge, queriedEdge, queryPathStartIndex);
          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << queryEdge->toString(state.query) << " ?= " << queriedEdge->toString(state.queried) << "  = " << result << std::endl;

          if (result > 0) {
            if (!Targeted) {
              // queried edge's path is equal in length or longer and matches
              if (state.usedQueriedEdges[queriedEdge->index()].use != UsedEdge::NotUsed)
                state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), UsedEdge::Both, queryEdge->pathSize());
              else {
                enum UsedEdge::Use use = (result == 1) ? UsedEdge::Forward : UsedEdge::Backward;
                state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), use, queryEdge->pathSize() - queryPathStartIndex);
              }
            } else {
              VertexType *mappedQueriedTarget = (queryEdge->source() == queryEdge->target()) ? 
                  queriedVertex : state.mapping[queryEdge->other(queryVertex)->index()];

              //std::cout << "queriedTarget = " << state.queried->vertexIndex(queriedTarget) << std::endl;
              //std::cout << "queriedTarget = " << state.queried->vertexIndex(queriedEdge->other(queriedVertex)) << std::endl;
              if (mappedQueriedTarget != queriedEdge->other(queriedVertex)) {
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << 1 << std::endl;
                continue; 
              }
              if (queryEdge->pathSize() - queryPathStartIndex != queriedEdge->pathSize()) {
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << 2 << std::endl;
                return false;
              }
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "FFFFFFFFFFFFG: " << queriedVertex->index() << " = " << queryEdge->other(queryVertex)->index() + 1 << std::endl;
              state.usedQueriedVertices[queriedVertex->index()] = queryEdge->other(queryVertex)->index() + 1;
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "HHHHHHHHHHhHHH: " << queriedEdge->other(queriedVertex)->index() << " = " << queryEdge->other(queryVertex)->index() + 1 << std::endl;
              state.usedQueriedVertices[queriedEdge->other(queriedVertex)->index()] = queryEdge->other(queryVertex)->index() + 1;
            }
            return true;
          } else if (result == -1) {
            // query edge's path is longer but there is a match up to queried path's length 
            // check if the queried target vertex is a match for query edge's path
            VertexType *queriedTarget = queriedEdge->other(queriedVertex);;
            if (!queriedTarget)
              continue;

            // skip used vertices
            if (state.usedQueriedVertices[queriedTarget->index()]) {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "YYYYYYYY" << std::endl;
              continue;
            }
 
            // check if the semantics of the target vertex match the next
            // element is the query's path
            if (!semanticsMatch(queryEdge->path()[queryPathStartIndex + queriedEdge->pathSize()], queriedTarget->semantics()))
              continue;

            // queriedTarget can not be in the path
            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), queriedTarget) != state.queriedPath.end())
              continue;

            // check if the query path is fully matched
            if (!Targeted && (queryEdge->pathSize() == queryPathStartIndex + queriedEdge->pathSize() + 1))
                return true;

            // recurse
            queryPathStartIndex += queriedEdge->pathSize() + 1;
            if (!Targeted) {
              enum UsedEdge::Use use = forwardDirection(queryVertex, queriedVertex, queryEdge, queriedEdge) ? UsedEdge::Forward : UsedEdge::Backward;
              state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), use, queriedEdge->pathSize());
              state.usedQueriedVertices[queriedVertex->index()] = state.query->numVertices() + queryVertex->index();
            } else {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "GGGGGGGGGGGGGGG: " << queriedVertex->index() << " = " << queryEdge->other(queryVertex)->index() + 1 << std::endl;
              state.usedQueriedVertices[queriedVertex->index()] = queryEdge->other(queryVertex)->index() + 1;
            }

            return edgesMatchDFS<Targeted>(state, queryEdge, queryVertex, queriedTarget, queryPathStartIndex);
          }
        }
      
        return false;
      }

      /**
       * Check if the untargeted edges around @p queryVertex match edges around
       * @p queriedVertex. This function makes two passes over the the queried
       * vertex's edges. This is done to reduce the number of states needing to
       * be explored later on.
       *
       * In the first pass the untargeted query edges are compared with
       * untargeted queried edges. This is a very simple process, either the
       * edges match or they don't. For two edges to match the following must
       * hold:
       * @li queryEdge->pathSize() <= queriedEdge->pathSize()
       * @li for each element in query's edge path: queryEdge->path()[i] == 
       *     queriedEdge->path()[i]
       * During this pass, the used queried edges are tracked to ensure that
       * two query edges are not matched with the same queried edge. Since the
       * vertex edges are ordered by decreasing path length, only one iteration
       * over the queried edges are needed. The matched query edges are also
       * tracked since the second pass can ignore these. If after the first pass
       * all query edges are matched, the second pass can be omitted.
       *
       * In the second pass the remaining unmatched query edges are matched with
       * targted queried edges. In this case there are three possible outcomes:
       * @li The edges don't match (i.e. there is an element in query's path
       *     that doesn't match the corresponding element in queried's edge).
       * @li queryEdge->pathSize() <= queriedEdge->pathSize() and for each
       *     element in the query's edge path: queryEdge->path()[i] == 
       *     queriedEdge->path()[i]
       * @li queryEdge->pathSize() > queriedEdge->pathSize() and for each
       *     element of queried edge's path: queryEdge->path()[i] == 
       *     queriedEdge->path()[i]
       *
       * The first two cases are simple. However, the last case is more
       * complicated. In this case, a depth-first search is initiated to try and
       * match the remaining elements in the query edge with edges in the
       * queried edge's target vertex. In the example below, we are trying to
       * match the edges around query vertex with index 1 to queried vertex
       * with index 5.
       * @code
       * queryEdge: { 1, 0, [ 1, 1, 1, 1, 1, 1 ] }
       * queriedEdge: { 5, 6, [ 1, 1, 1 ] }
       * @endcode
       * Here the query edge's path first 3 (0-2) element match the queried
       * edge's path. Next, the queried target vertex with index 6 is compared
       * with element 3 of the query edge's path. If this comparison is a
       * match, the process continues by doing a recursive depth first search
       * for the remaining 2 elements (4-5) of query edge's path starting from
       * queried vertex with index 6. For example, if queried vertex 6 has an
       * edge like:
       * @code
       * queriedEdge: { 6, 7, [ 1, 1, 1, 1] }
       * @endcode
       * The remaing 2 elements (4-5) from the query edges's path are matched
       * and we conclude that the edges around the query vertex and initial
       * queried vertex match. If this last match was a partial match again, the
       * recursive process would continue from queried vertex 7 untill there is
       * a match or mismatch.
       *
       * @param state The state s.
       * @param queryVertex The query vertex.
       * @param queriedVertex The queried vertex.
       *
       * @return True if the edges around @p queryVertex match edges around @p
       *         queriedVertex.
       */
      template<typename StateType, typename VertexType>
      bool checkUntargetedEdges(StateType &state, VertexType *queryVertex,
          VertexType *queriedVertex)
      {
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::VertexEdgesType VertexEdgesType;

        const VertexEdgesType &queryEdges = queryVertex->edges();
        // make a copy of the queried edges, we need to shuffle this later
        VertexEdgesType queriedEdges = queriedVertex->edges();

        // since there are two passes over the queried edges, we need to keep track of already matched edges
        std::vector<bool> matchedQueryEdges(queryVertex->numEdges(), false);

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << "checkUntargetedEdges()" << std::endl;

        unsigned int numUntargetedQueryEdges = 0;

        //
        // PASS 1: match untargeted query edges with untargeted queried edges
        //         first this reduces the number of states
        //
        // keep track of edges in queried that are already used, this avoid that parrallel
        // edges in the query would (both) match a single edge in queried
        // (this happens for parrallel self-loop and vertex-vertex edges (i.e. rings) )
        // check if each edge in queryVertex has a matching edge in queried
        for (std::size_t i = 0; i < queryVertex->numEdges(); ++i) {
          EdgeType *queryEdge = state.query->edge(queryEdges[i]);

          // check untargeted query edges only
          if (queryEdge->target())
            continue;

          numUntargetedQueryEdges++;

          // look for matching edge in queriedEdges
          for (std::size_t j = 0; j < queriedVertex->numEdges(); ++j) {
            EdgeType *queriedEdge = state.queried->edge(queriedEdges[j]);
            // skip targeted queried edges
            if (queriedEdge->target())
              continue;

            std::size_t queriedEdgeIndex = queriedEdge->index();
            // skip already used queried edges
            if (state.usedQueriedEdges[queriedEdgeIndex].use != UsedEdge::NotUsed)
              continue;
            
            int result = edgesMatch(queryVertex, queriedVertex, queryEdge, queriedEdge);
            if (result > 0) {
              // update matched query edges
              matchedQueryEdges[indexOf(queryEdges, queryEdge->index())] = true;
              state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), UsedEdge::Forward, queryEdge->pathSize());
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << queryEdge->toString(state.query) << " == " << queriedEdge->toString(state.queried) << "  = 1" << std::endl;
              break;
            }
            
            if (DEBUG_MOLDB_ISOMORPHISM_H)
              std::cout << queryEdge->toString(state.query) << " != " << queriedEdge->toString(state.queried) << std::endl;
          }
        }

        // if all edges are matched, we're done
        if (std::count(matchedQueryEdges.begin(), matchedQueryEdges.end(), true) == numUntargetedQueryEdges)
          return true;

        //
        // PASS 2: match untargeted query edges with targeted query edges
        //
        // check if each edge in queryVertex has a matching edge in queried
        for (std::size_t i = 0; i < queriedVertex->numEdges(); ++i) {
          //std::cout << "start" << std::endl;

          for (std::size_t j = 0; j < queryVertex->numEdges(); ++j) {
            EdgeType *queryEdge = state.query->edge(queryEdges[j]);

            // check untargeted query edges only
            if (queryEdge->target())
              continue;

            // look for matching edge in queriedEdges
            bool foundMatchingEdge = false;
            for (std::size_t k = 0; k < queriedVertex->numEdges(); ++k) {
              EdgeType *queriedEdge = state.queried->edge(queriedEdges[k]);
              // skip untargeted queried edges
              if (!queriedEdge->target())
                continue;
              std::size_t queriedEdgeIndex = queriedEdge->index();
              // skip already used queried edges
              bool forceBackward = false;
              if (state.usedQueriedEdges[queriedEdgeIndex].use != UsedEdge::NotUsed) {
                if (queriedEdge->source() == queriedEdge->target()) {
                  if (state.usedQueriedEdges[queriedEdgeIndex].use == UsedEdge::Both)
                    continue;
                  if (state.usedQueriedEdges[queriedEdgeIndex].use == UsedEdge::Forward) {
                    if (queriedEdge->pathSize() - state.usedQueriedEdges[queriedEdgeIndex].count < queryEdge->pathSize())
                      continue;
                    forceBackward = true;
                    if (DEBUG_MOLDB_ISOMORPHISM_H)
                      std::cout << "force backward" << std::endl;
                  }
                } else {
                  enum UsedEdge::Use use = forwardDirection(queryVertex, queriedVertex, queryEdge, queriedEdge) ? UsedEdge::Forward : UsedEdge::Backward;
                  if (state.usedQueriedEdges[queriedEdgeIndex].use == use)
                    continue;
                  if (queriedEdge->pathSize() - state.usedQueriedEdges[queriedEdgeIndex].count < queryEdge->pathSize())
                    continue;
                }
              }
 
              int result = edgesMatch(queryVertex, queriedVertex, queryEdge, queriedEdge, forceBackward);
              if (result > 0) {
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << queryEdge->toString(state.query) << " == " << queriedEdge->toString(state.queried) << "  = 1" << std::endl;
                // queried edge's path is equal in length or longer
                foundMatchingEdge = true;
                if (state.usedQueriedEdges[queriedEdge->index()].use != UsedEdge::NotUsed)
                  state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), UsedEdge::Both, queryEdge->pathSize());
                else {
                  enum UsedEdge::Use use = (result == 1) ? UsedEdge::Forward : UsedEdge::Backward;
                  state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), use, queryEdge->pathSize());
                }
                break;
              } else if (result == -1) {
                // skip used vertices
                VertexType *queriedTarget = queriedEdge->other(queriedVertex);
                if (state.usedQueriedVertices[queriedTarget->index()])
                  continue;
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << queryEdge->toString(state.query) << " >= " << queriedEdge->toString(state.queried) << "  = -1" << std::endl;
                // query edge's path is longer but there is a match up to queried path's length
                // check if the queried target vertex is a match for query edge's path
                if (!semanticsMatch(queryEdge->path()[queriedEdge->pathSize()], queriedTarget->semantics()))
                  continue;
                enum UsedEdge::Use use = forwardDirection(queryVertex, queriedVertex, queryEdge, queriedEdge) ? UsedEdge::Forward : UsedEdge::Backward;
                state.usedQueriedEdges[queriedEdge->index()] = UsedEdge(state.queriedPath.size(), use, queriedEdge->pathSize());
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << "XXXXXXXXXXXXXX: usedQueriedVertices[" << queriedEdge->other(queriedVertex)->index() << "] = " << state.query->numVertices() + queryVertex->index() + 1 << std::endl;
                state.usedQueriedVertices[queriedEdge->other(queriedVertex)->index()] = state.query->numVertices() + queryVertex->index() + 1;
                if (queryEdge->pathSize() == queriedEdge->pathSize() + 1) {
                  foundMatchingEdge = true;
                  break;
                }
                std::size_t queryPathStartIndex = queriedEdge->pathSize() + 1;
                if (edgesMatchDFS<false>(state, queryEdge, queryVertex, queriedTarget, queryPathStartIndex)) {
                  foundMatchingEdge = true;
                  break;
                }
              }

              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << queryEdge->toString(state.query) << " != " << queriedEdge->toString(state.queried) << "  = 0" << std::endl;
            }

            // return false if there is no matching edge
            if (foundMatchingEdge)
              matchedQueryEdges[indexOf(queryEdges, queryEdge->index())] = true;
          }
        
          //if (DEBUG_MOLDB_ISOMORPHISM_H)
          //  print_vector("usedEdges", state.usedQueriedEdges);

          std::replace(state.usedQueriedVertices.begin(), state.usedQueriedVertices.end(),
                       static_cast<unsigned int>(state.query->numVertices() + queryVertex->index() + 1), static_cast<unsigned int>(0));
          std::replace_if(state.usedQueriedEdges.begin(), state.usedQueriedEdges.end(), CompareUsedEdgeDepth(state.queryPath.size()), UsedEdge());

          if (std::count(matchedQueryEdges.begin(), matchedQueryEdges.end(), true) == numUntargetedQueryEdges)
            return true;

          // cycle through queriedEdges
          // TODO only permutations of edges with the same length are needed...
          //      only permutations of targeted  edges are needed
          EdgeType *last = state.queried->edge(queriedEdges.back());
          queriedEdges.pop_back();
          queriedEdges.insert(queriedEdges.begin(), last->index());
        }
        
        
        return (std::count(matchedQueryEdges.begin(), matchedQueryEdges.end(), true) == numUntargetedQueryEdges);
      }


      template<typename StateType, typename VertexType>
      bool checkTargetedEdges(StateType &state, VertexType *queryVertex,
          VertexType *queriedVertex)
      {
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::VertexEdgesType VertexEdgesType;

        const VertexEdgesType &queryEdges = queryVertex->edges();

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << "checkTargetedEdges()" << std::endl;
        std::vector<bool> usedQueriedEdges(queriedVertex->numEdges(), false);

        // check if each edge in queryVertex has a matching edge in queried
        for (std::size_t i = 0; i < queryVertex->numEdges(); ++i) {
          EdgeType *queryEdge = state.query->edge(queryEdges[i]);

          // check targeted edges only
          if (!queryEdge->target())
            continue;

          // if queryEdge->source() == queryEdge->target()
          //     - queryVertex has a self-loop edge
          //     - there may be parallel self-loops
          // else 
          //     - queryEdge conncects two vertices in query, 
          //       we need to find an exact matching edge in queried
          VertexType *queriedTarget = (queryEdge->source() == queryEdge->target()) ? 
            queriedVertex : state.mapping[queryEdge->other(queryVertex)->index()];

          // queriedTarget is not mapped yet, edge will be checked later
          if (!queriedTarget)
            continue;

          // there may be parallel edges
          const VertexEdgesType &queriedEdges = queriedVertex->edges();

          if (DEBUG_MOLDB_ISOMORPHISM_H) {
            std::cout << "querySource = " << queryVertex->index() + 1 << std::endl;
            std::cout << "queryTarget = " << queryEdge->other(queryVertex)->index() + 1 << std::endl;
            std::cout << "queriedVertex = " << queriedVertex->index() + 1 << std::endl;
            std::cout << "queriedTarget = " << queriedTarget->index() + 1 << std::endl;

            std::cout << "queryEdges:" << std::endl;
            for (int k = 0; k < queryVertex->numEdges(); ++k)
              std::cout << "    " << state.query->edge(queryVertex->edges()[k])->toString(state.query) << std::endl;
            std::cout << "queriedEdges:" << std::endl;
            for (int k = 0; k < queriedVertex->numEdges(); ++k)
              std::cout << "    " << state.queried->edge(queriedVertex->edges()[k])->toString(state.queried) << std::endl;
          }


          // look for matching edge in queriedEdges
          bool foundMatchingEdge = false;
          for (std::size_t j = 0; j < queriedVertex->numEdges(); ++j) {
            EdgeType *queriedEdge = state.queried->edge(queriedEdges[j]);
            if (usedQueriedEdges[j])
              continue;
            if (DEBUG_MOLDB_ISOMORPHISM_H)
              std::cout << queryEdge->toString(state.query) << " =? " << queriedEdge->toString(state.queried) << std::endl;
            // skip already used vertices
            /*
            if (state.usedQueriedVertices[state.queried->vertexIndex(queriedEdge->other(queriedVertex))]) {
              std::cout << "skipping used vertex" << std::endl;
              continue;
            }
            */
            int result = edgesMatch(queryVertex, queriedVertex, queryEdge, queriedEdge);
            if (DEBUG_MOLDB_ISOMORPHISM_H)
              std::cout << "result = " << result << std::endl;
            if (result > 0) {
              // make sure the correct mapping is used
              if (queriedTarget != queriedEdge->other(queriedVertex))
                continue; 
              if (queryEdge->pathSize() != queriedEdge->pathSize())
                return false;
              foundMatchingEdge = true;
              state.usedQueriedVertices[queriedEdge->other(queriedVertex)->index()] = queryEdge->other(queryVertex)->index() + 1;
              usedQueriedEdges[j] = true;
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << queryEdge->toString(state.query) << " -> " << queriedEdge->toString(state.queried) << std::endl;
              break;
            } else if (result == -1) {
              // skip used vertices
              VertexType *queriedTarget2 = queriedEdge->other(queriedVertex);
              if (state.usedQueriedVertices[queriedTarget2->index()])
                continue;
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << queryEdge->toString(state.query) << " >= " << queriedEdge->toString(state.queried) << "  = -1" << std::endl;
              // query edge's path is longer but there is a match up to queried path's length
              // check if the queried target vertex is a match for query edge's path
              if (!semanticsMatch(queryEdge->path()[queriedEdge->pathSize()], queriedTarget2->semantics()))
                continue;
              std::size_t queryPathStartIndex = queriedEdge->pathSize() + 1;
              if (edgesMatchDFS<true>(state, queryEdge, queryVertex, queriedTarget2, queryPathStartIndex)) {
                foundMatchingEdge = true;
//              std::cout << "HHHHHHHHHHhHHH: " << state.queried->vertexIndex(queriedVertex) << " = " << state.query->vertexIndex(queryEdge->other(queryVertex)) + 1 << std::endl;
//                state.usedQueriedVertices[state.queried->vertexIndex(queriedEdge->other(queriedVertex))] = state.query->vertexIndex(queryEdge->other(queryVertex)) + 1;
                usedQueriedEdges[j] = true;
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << "edgesMatchDFS successful" << std::endl;
                break;
              }
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "edgesMatchDFS unsuccessful" << std::endl;
            }

          }
          
          // return false if there is no matching edge
          if (!foundMatchingEdge) {
            if (DEBUG_MOLDB_ISOMORPHISM_H)
              std::cout << "return false" << std::endl;
            return false;
          }
        }

        return true;
      }

      template<typename StateType, typename VertexType>
      bool checkEdges(StateType &state, VertexType *queryVertex,
          VertexType *queriedVertex)
      {
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::VertexEdgesType VertexEdgesType;
 
        const VertexEdgesType &queryEdges = queryVertex->edges();
        for (unsigned int i = 0; i < queryVertex->numEdges(); ++i) {
          EdgeType *queryEdge = state.query->edge(queryEdges[i]);
          VertexType *queryNbr = queryEdge->other(queryVertex);

          if (state.sharedState->queryMapping[queryNbr->index()]) {
            VertexType *queriedNbr = state.sharedState->queryMapping[queryNbr->index()];
            EdgeType *queriedEdge = queriedVertex->edge(state.queried, queriedNbr);

            if (!queriedEdge)
              return false;

            if (queryEdge->pathSize() != queriedEdge->pathSize())
              return false;
            if (!queryEdge->path().empty())
              if (queryEdge->path()[0] != queriedEdge->path()[0])
                return false;
          }
        }

        return true;
      }

      /**
       * The VF2 feasibility function F(s, n, m). The n and m vertex are not
       * needed as arguments since the @p state already contains them. Therefor
       * the feasibility function simply becomes F(s). If the state is not
       * feasible, the state is restored.
       *
       * A state s is feasible if and only if the following equalities hold:
       * @code
       * size( T1(s) ) <= size( T2(s) )
       * size( N1 - M1(s) - T1(s) ) <= size( N2 - M2 - T2(s) )
       * @endcode
       * 
       * Here, T1(s) and T2(s) are the terminal sets for the query and queried
       * graph respectivily. N1 and N2 are the sets of all vertices (again) for
       * the query and queried graph respectivily. M1(s) and M2(s) are the sets
       * of vertices that are already matched in the state s. The function size()
       * returns the size of the set and the operation '-' means set difference.
       *
       * @note Runtime complexity O(2n) ~ O(n) where n is the number of vertices
       *       already mapped in state s.
       *
       * @param state The state s.
       *
       * @return True if the state s is feasible.
       */
      /*
      template<typename StateType>
      bool isStateFeasible(StateType &state)
      {
        //
        // Feasibility rules for the VF2 algorithm:
        //
        //  size( T1(s) ) <= size( T2(s) )
        //
        //  size( N1 - M1(s) - T1(s) ) <= size( N2 - M2(s) - T2(s) )
        //

        // compute T1(s) size
        unsigned int numT1 = 0;
        for (std::size_t i = 0; i < state.query->numVertices(); ++i)
          if (isInTerminalSet(state.queried, state.queryDepths, state.queryPath, i)) // O(n)  n = # mapped vertices in state s
            numT1++;
        // compute T2(s) size
        unsigned int numT2 = 0;
        for (std::size_t i = 0; i < state.queried->numVertices(); ++i)
          if (isInTerminalSet(state.query, state.queriedDepths, state.queriedPath, i)) // O(n)  n = # mapped vertices in state s
            numT2++;

        // T1(s) > T2(s)
        if (numT1 > numT2) {
          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << "not feasible: T1(s) > T2(s)" << std::endl;
          backtrack(state);
          return false;
        }
        //  N1 - M1(s) - T1(s) > N2 - M2(s) - T2(s)
        if ((state.query->numVertices() - state.queryPath.size() - numT1) > (state.queried->numVertices() - state.queriedPath.size() - numT2)) {
          if (DEBUG_MOLDB_ISOMORPHISM_H) {
            std::cout << "not feasible: N1 - M1(s) - T1(s) > N2 - M2(s) - T2(s)" << std::endl;
            std::cout << "not feasible: " << state.query->numVertices() << " - " << state.queryPath.size() << " - " << numT1 << " > "
                      << state.queried->numVertices() << " - " << state.queriedPath.size() << " - " << numT2 << std::endl;
          }
          backtrack(state);
          return false;
        }

        return true;
      }
      */

      template<typename StateType, typename CandidateType>
      bool isStateFeasible(const StateType &state, const CandidateType &candidate)
      {
        typedef typename StateType::VertexType VertexType;

        VertexType *queryVertex = candidate.first;
        VertexType *queriedVertex = candidate.second;
    
        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << blue << "isStateFeasible(" << queryVertex->index() << " -> " << queriedVertex->index() << ")" << normal << std::endl;

        if (!verticesMatch(queryVertex, queriedVertex))
          return false;

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << red << "    vertices match..." << normal << std::endl;

        std::size_t queryTerminalNbrCount = 0;
        std::size_t queriedTerminalNbrCount = 0;
        std::size_t queryNewNbrCount = 0;
        std::size_t queriedNewNbrCount = 0;

        for (std::size_t i = 0; i < queryVertex->numEdges(); ++i) {
          VertexType *queryNbr = state.query->edge(queryVertex->edges()[i])->other(queryVertex);
          if (state.sharedState->queryMapping[queryNbr->index()])
            continue;
          if (state.sharedState->queryDepths[queryNbr->index()]) {
            ++queryTerminalNbrCount;
          } else {
            ++queryNewNbrCount;
          }
        }

        for (std::size_t i = 0; i < queriedVertex->numEdges(); ++i) {
          VertexType *queriedNbr = state.queried->edge(queriedVertex->edges()[i])->other(queriedVertex);
          if (state.sharedState->queriedMapping[queriedNbr->index()])
            continue;
          if (state.sharedState->queriedDepths[queriedNbr->index()]) {
            ++queriedTerminalNbrCount;
          } else {
            ++queriedNewNbrCount;
          }
        }

        return (queryTerminalNbrCount <= queriedTerminalNbrCount) && (queryNewNbrCount <= queriedNewNbrCount);
      }

      /**
       * Compute the state s' by adding the pair (n, m) to state s.
       *
       * @note Runtime complexity O(nv + mw) where n is the number of vertices
       *       in the query, m is the number of vertices in the queried, v is
       *       the number of edges in the query vertex and w is the number of
       *       edges in the queried vertex.
       *
       * @param state The state s, will be updated to state s'.
       * @param queryVertex The vertex n from the query.
       * @param queriedVertex The vertex m from the queried.
       */
      template<typename StateType, typename CandidateType>
      void updateState(StateType &state, const CandidateType &candidate)
      {
        typedef typename StateType::VertexType VertexType;
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::VertexEdgesType VertexEdgesType;

        ++state.mappingSize;
        state.lastMapped = candidate;
          
        VertexType *queryVertex = candidate.first;
        VertexType *queriedVertex = candidate.second;

        // update mappings
        state.sharedState->queryMapping[queryVertex->index()] = queriedVertex;
        state.sharedState->queriedMapping[queriedVertex->index()] = queryVertex;

        //
        // update queryDepths
        //
        // add queryVertex's depth if not already set
        std::size_t queryVertexIndex = queryVertex->index();
        if (!state.sharedState->queryDepths[queryVertexIndex])
          state.sharedState->queryDepths[queryVertexIndex] = state.mappingSize;

        // add queryVertex's neighbors' depths if not already set
        const VertexEdgesType &queryEdges = queryVertex->edges(); 
        for (unsigned int i = 0; i < queryVertex->numEdges(); ++i) { // O(nm)  m = # edges in query vertex
          EdgeType *queryEdge = state.query->edge(queryEdges[i]);
          if (!queryEdge->target())
            continue;
          std::size_t index = queryEdge->other(queryVertex)->index(); // O(n)  n = # vertices in query
          if (!state.sharedState->queryDepths[index]) {
            state.sharedState->queryDepths[index] = state.mappingSize;
            ++state.queryDepthsSize;
          }
        }

        //
        // update queriedDepths
        //
        // add queriedVertex's depth if not already set
        std::size_t queriedVertexIndex = queriedVertex->index();
        if (!state.sharedState->queriedDepths[queriedVertexIndex])
          state.sharedState->queriedDepths[queriedVertexIndex] = state.mappingSize;

        // add queriedVertex's neighbors' depths if not already set
        const VertexEdgesType &queriedEdges = queriedVertex->edges();
        for (unsigned int i = 0; i < queriedVertex->numEdges(); ++i) { // O(nm)  m = # vertices in queried vertex
          EdgeType *queriedEdge = state.queried->edge(queriedEdges[i]);
          if (!queriedEdge->target())
            continue;
          std::size_t index = queriedEdge->other(queriedVertex)->index(); // O(n)  n = # vertices in queried
          if (!state.sharedState->queriedDepths[index])
            state.sharedState->queriedDepths[index] = state.mappingSize;
        }

        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          std::cout << "    queryDepths:   ";
          for (unsigned int i = 0; i < state.query->numVertices(); ++i)
            std::cout << state.sharedState->queryDepths[i] << " ";
          std::cout << std::endl;
          std::cout << "    queriedDepths: ";
          for (unsigned int i = 0; i < state.queried->numVertices(); ++i)
            std::cout << state.sharedState->queriedDepths[i] << " ";
          std::cout << std::endl;
          std::cout << "    queryMapping: ";
          for (unsigned int i = 0; i < state.sharedState->queryMapping.size(); ++i)
            if (state.sharedState->queryMapping[i])
              std::cout << state.sharedState->queryMapping[i]->index() << " ";
            else
              std::cout << "0 ";
          std::cout << std::endl;
          std::cout << "    queriedMapping: ";
          for (unsigned int i = 0; i < state.sharedState->queriedMapping.size(); ++i)
            if (state.sharedState->queriedMapping[i])
              std::cout << state.sharedState->queriedMapping[i]->index() << " ";
            else
              std::cout << "0 ";
          std::cout << std::endl;
        }
 
      }

      /**
       * Check if the vertex with index @p i is in the terminal set. A vertex
       * is in the terminal set if the depth is non-zero and the vertex is not
       * already mapped (not in the Mx(s) set).
       *
       * @note Runtime complexity O(n) where n is the number of vertices already
       *       mapped (size( Mx(s))
       *
       * @param state The state s.
       * @param g The graph.
       * @param depths The std::vector with depths.
       * @param path The path of mapped vertices.
       * @param i The vertex index.
       *
       * @return True if vertex with index @p i is in the terminal set Tx(s).
       */
      template<typename GraphType, typename VertexType>
      bool isInTerminalSet(const GraphType *g, 
          const std::vector<unsigned int> &depths,
          const std::vector<VertexType*> &path, std::size_t i)
      {
        if (!depths[i])
          return false;

        // NOTE: this O(n) search can be avoided at the cost of O(n) memory
        if (std::find(path.begin(), path.end(), g->vertices()[i]) != path.end())
          return false;

        return true;
      }


      /**
       * Restore the @p state data structure.
       *
       * @note Runtime cemplexity O(2n + m) ~ O(n + m) where n is the number of
       *       vertices in the query and m is the number of vertices in the
       *       queried graph.
       *
       * @param state The state s to restore.
       */
      template<typename StateType>
      void backtrack(StateType &state)
      {
        typedef typename StateType::VertexType VertexType;
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::CandidateType CandidateType;

        VertexType *queryVertex = state.lastMapped.first;
        VertexType *queriedVertex = state.lastMapped.second;

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << red << "backtrack... " << normal << queryVertex->index() << std::endl;

        if (state.sharedState->queryDepths[queryVertex->index()] == state.mappingSize)
          state.sharedState->queryDepths[queryVertex->index()] = 0;

        for (std::size_t i = 0; i < queryVertex->numEdges(); ++i) {
          EdgeType *queryEdge = state.query->edge(queryVertex->edges()[i]);
          VertexType *queryNbr = queryEdge->other(queryVertex);
          if (queryEdge->target()) {
            std::replace(state.usedQueriedVertices.begin(), state.usedQueriedVertices.end(), 
                         static_cast<unsigned int>(queryNbr->index()), 
                         static_cast<unsigned int>(0));
            if (state.sharedState->queryDepths[queryNbr->index()] == state.mappingSize)
              state.sharedState->queryDepths[queryNbr->index()] = 0;
          }
        }

        if (state.sharedState->queriedDepths[queriedVertex->index()] == state.mappingSize)
          state.sharedState->queriedDepths[queriedVertex->index()] = 0;

        for (std::size_t i = 0; i < queriedVertex->numEdges(); ++i) {
          EdgeType *queriedEdge = state.queried->edge(queriedVertex->edges()[i]);
          VertexType *queriedNbr = queriedEdge->other(queriedVertex);
          if (queriedEdge->target()) {
            if (state.sharedState->queriedDepths[queriedNbr->index()] == state.mappingSize)
              state.sharedState->queriedDepths[queriedNbr->index()] = 0;
          }
        }

        state.sharedState->queryMapping[queryVertex->index()] = 0;
        state.sharedState->queriedMapping[queriedVertex->index()] = 0;

        --state.mappingSize;
        state.lastMapped = CandidateType(0, 0);
      }


      /**
       * Compute the candidates to consider to match next. This is the P(s)
       * function from the VF2 algorithm.
       *
       * The candidates are computed using this formula:
       * @code
       * P(s) = T2(s) x { min T1(s) }
       * @endcode
       *
       * Here, T1(s) and T2(s) are the terminal sets for the query and queried
       * graph respectivily. The min operator takes the smalles element from a
       * set, any ordering will do and in this implementation the index is
       * used. The { } around min T1(s) makes the element a set again. The 'x'
       * operator means cartesian product of two sets (i.e. {1, 2} x {3, 4} =
       * { {1, 3}, {1, 4}, {2, 3}, {2, 4} })
       *
       * If the function above returns an empty set of pairs, the following
       * function is used to obtain candidate pairs.
       * @code
       * P(s) = (N2 - M2(s)) x { min (N1 - M1(s)) }
       * @endcode
       *
       * Here, N1 and N2 are the sets of all vertices for the query and queried
       * respectivily. M1(s) and M2(s) are the sets of the already mapped
       * vertices in state s. The '-' operation means set difference (i.e.
       * {1, 2, 3, 4} - {2, 3} = [1, 4}).
       *
       * Before a pair is added to the candidates list, a check is performed to
       * ensure the vertex semantic attributes match.
       * 
       * @note Runtime complexity O(2np + 2mp) ~ O(np + mp) usually O(mp)
       *       though. Here n is the number of vertices in the query, m is the
       *       number of vertices in the queried and p is the number of already
       *       mapped vertices in state s.
       *
       * @param state The state s.
       * @param candidates The std::vector in which to store the found pairs.
       */
      /*
      template<typename StateType, typename CandidateType>
      void computeCandidates(StateType &state,
          std::vector<CandidateType> &candidates)
      {
        typedef typename StateType::VertexType VertexType;

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << "computeCandidates()" << std::endl;

        // The terminology used in these comments is taken from the VF2 paper.
        // Since molecules are undirected graphs, vertices have only one set of
        // edges.
        //
        //  G1 : The query graph (i.e. state.query)
        //  G2 : The queried graph (i.e. state.queried)
        //
        //  N1 : The set of query vertices in G1
        //  N2 : The set of queried vertices in G2
        //
        //  M1(s) : The set of already mapped query vertices (i.e. state.queryPath)
        //  M2(s) : The set of already mapped queried vertices (i.e. state.queriedPath)
        //
        // Variables from the original C++ implementation referenced in paper:
        //
        //  core_1 : M1(s) or state.queryPath
        //  core_2 : M2(s) or state.queriedPath
        //
        //  in_1, out_1 : state.queryDepths
        //  in_2, out_2 : state.queriedDepths
        //
        //
        //  P(s) = T1(s) x { min T2(s) }
        //
        //  P(s) : The set of mapping candidates to consider for state s
        //
        //  T1(s) : The set of "terminal" query vertices for state s
        //  T2(s) : The set of "terminal" queried vertices for state s
        //
        //  min T2(s) : The element from T2(s) with the lowest index. The { } brackets make the element a set again.
        //
        //  x : means cartesian product (e.g. {1, 2} x {3, 4} = {{1, 3}, {1, 4}, {2, 3}, {2, 4}}
        //
        // T1(s) is computed using state.queryDepths and the current mapping (i.e. state.queryPath).
        // A query vertex with index n is said to be in the T1(s) set if state.queryDepths[n] is non-zero
        // and n is not part of the current mapping (i.e. n is not in state.queryPath). T2(s) is computed in
        // a similar way.

        // compute min T2(s)
        VertexType *minT2 = 0;
        for (std::size_t i = 0; i < state.queried->numVertices(); ++i) { // O(np)  n = # vertices in query
          if (isInTerminalSet(state.queried, state.queriedDepths, state.queriedPath, i)) { // O(p)  p = # vertices mapped in state
            minT2 = state.queried->vertices()[i];
            break;
          }
        }
        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << "P(s) = T1(s) x { min T2(s) } = { ";
        // compute T1(s)
        if (minT2) {
          for (std::size_t i = 0; i < state.query->numVertices(); ++i) { // O(mp)  n = # vertices in query
            VertexType *T1_i = state.query->vertices()[i];
            if (!isInTerminalSet(state.query, state.queryDepths, state.queryPath, i)) // O(p)  p = # vertices mapped in state
              continue;
            if (verticesMatch(T1_i, minT2)) {
              candidates.push_back(CandidateType(T1_i, minT2));
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << T1_i->index() << " ";
            }
          }
        }
        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          if (minT2)
            std::cout << "} x { " << minT2->index() << " }" << std::endl;
          else
            std::cout << "} x { }" << std::endl;
        }

        // If the P(s) set from above is empty, use this set:
        //
        //  P(s) = (N1 - M1(s)) x { min (N2 - M2(s)) }
        //
        //  N1, N2 : all vertices in the query and queried respectivly
        //
        //  M1, M2 : all vertices in the mapping for the query and queried respectivly
        //
        //  - : simply means minus (e.g. {1, 2, 3, 4} - {2, 3} = {1, 4})
        //
        // N1 - M1(s) can simply be computed by ensuring that the quey vertex is not in
        // state.queryPath

        // compute min(N2 - M2(s))
        if (candidates.empty()) {
          VertexType *minNM2 = 0;
          for (std::size_t i = 0; i < state.queried->numVertices(); ++i) { // O(np)  n = # vertices in query
            if (std::find(state.queriedPath.begin(), state.queriedPath.end(), state.queried->vertices()[i]) == state.queriedPath.end()) { // O(p)  p = # vertices mapped in state
              minNM2 = state.queried->vertices()[i];
              break;
            }
          }
          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << "P(s) = (N1 - M1(s)) x { min (N2 - M2(s)) } = { ";
          // compute (N1 - M1(s))
          if (minNM2) {
            for (std::size_t i = 0; i < state.query->numVertices(); ++i) { // O(mp)  m = # vertices in query
              VertexType *NM1_i = state.query->vertices()[i];
              if (std::find(state.queryPath.begin(), state.queryPath.end(), NM1_i) != state.queryPath.end()) // O(p)  p = # vertices mapped in state
                continue;
              if (verticesMatch(NM1_i, minNM2)) {
                candidates.push_back(CandidateType(NM1_i, minNM2));
                if (DEBUG_MOLDB_ISOMORPHISM_H)
                  std::cout << NM1_i->index() << " ";
              }
            }
          }
          if (DEBUG_MOLDB_ISOMORPHISM_H) {
            if (minNM2)
              std::cout << "} x { " << minNM2->index() << " }" << std::endl;
            else
              std::cout << "} x { }" << std::endl;
          }
        }
 
      }
      */

      template<typename StateType, typename CandidateType>
      CandidateType nextCandidate(StateType &state, const CandidateType &lastCandidate)
      {
        std::size_t lastQueryVertex = lastCandidate.first ? lastCandidate.first->index() : 0;
        std::size_t lastQueriedVertex = lastCandidate.second ? lastCandidate.second->index() + 1: 0;

        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          std::cout << blue << "nextCandidate()" << normal << std::endl;
          std::cout << "    queryDepths:   ";
          for (unsigned int i = 0; i < state.query->numVertices(); ++i)
            std::cout << state.sharedState->queryDepths[i] << " ";
          std::cout << std::endl;
          std::cout << "    queriedDepths: ";
          for (unsigned int i = 0; i < state.queried->numVertices(); ++i)
            std::cout << state.sharedState->queriedDepths[i] << " ";
          std::cout << std::endl;
          std::cout << "    queryMapping: ";
          for (unsigned int i = 0; i < state.sharedState->queryMapping.size(); ++i)
            if (state.sharedState->queryMapping[i])
              std::cout << state.sharedState->queryMapping[i]->index() << " ";
            else
              std::cout << "0 ";
          std::cout << std::endl;
          std::cout << "    queriedMapping: ";
          for (unsigned int i = 0; i < state.sharedState->queriedMapping.size(); ++i)
            if (state.sharedState->queriedMapping[i])
              std::cout << state.sharedState->queriedMapping[i]->index() << " ";
            else
              std::cout << "0 ";
          std::cout << std::endl;
        }
 
        std::size_t querySize = state.query->numVertices();
        std::size_t queriedSize = state.queried->numVertices();

        if (state.queryDepthsSize > state.mappingSize && state.queriedDepthsSize > state.mappingSize) {
          while (lastQueryVertex < querySize && (state.sharedState->queryMapping[lastQueryVertex] || !state.sharedState->queryDepths[lastQueryVertex])) {
            lastQueryVertex++;
            lastQueriedVertex = 0;
          }
        } else {
          while(lastQueryVertex < querySize && state.sharedState->queryMapping[lastQueryVertex]) {
            lastQueryVertex++;
            lastQueriedVertex = 0;
          }
        }

        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          std::cout << blue << "lastQueryVertex: " << lastQueryVertex << normal << std::endl;
          std::cout << blue << "lastQueriedVertex: " << lastQueriedVertex << normal << std::endl;
        }

        if (state.queryDepthsSize > state.mappingSize && state.queriedDepthsSize > state.mappingSize) {
          while (lastQueriedVertex < queriedSize && (state.sharedState->queriedMapping[lastQueriedVertex] || !state.sharedState->queriedDepths[lastQueriedVertex])) 
            lastQueriedVertex++;
        } else {
          while(lastQueriedVertex < queriedSize && state.sharedState->queriedMapping[lastQueriedVertex])
            lastQueriedVertex++;
        }

        if (DEBUG_MOLDB_ISOMORPHISM_H)
          std::cout << blue << "lastQueriedVertex: " << lastQueriedVertex << normal << std::endl;

        if (lastQueryVertex < querySize && lastQueriedVertex < queriedSize)
          return CandidateType(state.query->vertex(lastQueryVertex), state.queried->vertex(lastQueriedVertex));

        return CandidateType();
      }

      template<bool ReducedGraph, typename StateType>
      bool checkForMapping(StateType &state)
      {
        // Check if there is a mapping found
        bool match = (state.mappingSize == state.query->numVertices());

        // if we already found a match, return
        /* FIXME
        if (ReducedGraph && match) {
          bool allUntargetedEdgesMatch = true;
          for (std::size_t i = 0; i < state.queryPath.size(); ++i) {
            if (!checkUntargetedEdges(state, state.queryPath[i], state.queriedPath[i]))
              allUntargetedEdgesMatch = false;
            state.usedQueriedEdges.clear();
            state.usedQueriedEdges.resize(state.queried->numEdges(), UsedEdge());
          }
          if (allUntargetedEdgesMatch) {
            if (DEBUG_MOLDB_ISOMORPHISM_H)
              std::cout << green << "-----------------> MATCH" << normal << std::endl;
            return true;
          }

          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << "    untargeted edges do not match..." << std::endl;

          return false;
        }
        */

        if (DEBUG_MOLDB_ISOMORPHISM_H && match)
          std::cout << green << "-----------------> MATCH" << normal << std::endl;
        return match;
      }

      /**
       * This is the main VF2 algorithm loop.
       *
       * @code
       * PROCEDURE Match(s)
       *     IF M(s) covers all the nodes of G2 THEN
       *         OUTPUT M(s)
       *     ELSE
       *         Compute the set P(s) of the pairs candidate for inclusion in M(s)
       *         FOREACH (n, m) in P(s)
       *             Compute the state s' obtained by adding (n, m) to M(s)
       *             Check if the edges around n match edgess in m
       *             IF F(s) THEN
       *                 CALL Match(s')
       *             END IF
       *         END FOREACH
       *         Restore data structures
       *     END IF
       * END PROCEDURE
       * @endcode
       *
       * In the pseudo-code above, Match(s) is the isomorphismDFS() function.
       * P(s) is computed using the computeCandidates() function.
       * F(s) is computed using the isStateFeasible() function.
       * The state s' is obtained from state s using the updateState() function.
       * The data structures are restored using the backtrack() function().
       *
       * @note Runtime complexity O(2^n) where n is the number of vertices in
       *       the query. This is a recursive function. See the
       *       GraphIsomorphismVF2 class for a detailed discussion of
       *       the actual runtime.
       *  
       * @param state The state s.
       */
      template<bool ReducedGraph, typename StateType>
      bool isomorphismDFS(StateType &state)
      {
        typedef typename StateType::GraphType GraphType;
        typedef typename StateType::VertexType VertexType;
        typedef typename StateType::EdgeType EdgeType;
        typedef typename StateType::CandidateType CandidateType;

        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          if (state.mappingSize)
            std::cout << "isomorphismDFS(" << state.lastMapped.first->index() << ", " << state.lastMapped.second->index() << ")" << std::endl;
          else
            std::cout << "isomorphismDFS(-1, -1)" << std::endl;
        }

        if (checkForMapping<ReducedGraph>(state))
          return true;

        CandidateType lastCandidate = std::make_pair<VertexType*, VertexType*>(0, 0);

        bool match = false;
        while (!match) {
         if (DEBUG_MOLDB_ISOMORPHISM_H && lastCandidate.first)
            std::cout << green << "lastCandidate: " << lastCandidate.first->index() << " -> " << lastCandidate.second->index() << normal << std::endl;

         CandidateType candidate = nextCandidate(state, lastCandidate);

          if (!candidate.first)
            return false;

          lastCandidate = candidate;

          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << yellow << "candidate: " << candidate.first->index() << " -> " << candidate.second->index() << normal << std::endl;

         
          // check feasibility function F(s)
          bool feasible = isStateFeasible(state, candidate);
          //if (ReducedGraph && feasible) {
            /*
            // check if the edges from queryVertex to already matched queried vertices match
            // alse check if the queryVertex self-loop and non-targeted edges match
            if (!checkTargetedEdges(state, candidate.first, candidate.second)) {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "    targeted edges do not match..." << std::endl;
              feasible = false;
            }
            */
          //} else {
            if (!checkEdges(state, candidate.first, candidate.second)) {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "    edges do not match..." << std::endl;
              feasible = false;
            }
          //}

          if (feasible) {
            StateType nextState(state);
            // compute state s' obtained by adding (n, m) to state s
            updateState(nextState, candidate);
            match = isomorphismDFS<ReducedGraph>(nextState);
            backtrack(nextState);
          }
        }

        return match;


        // compute the candidates using the P(s) function
        /*
        std::vector<CandidateType> candidates;
        computeCandidates(state, candidates);


        if (DEBUG_MOLDB_ISOMORPHISM_H) {
          std::cout << "Candidates:" << std::endl;
          for (std::size_t i = 0; i < candidates.size(); ++i)
            std::cout << "        " << candidates[i].first->index() << " -> " << candidates[i].second->index() << std::endl;
        }

        // do the mapping by checking the candidates
        while (candidates.size()) {
          if (state.match)
            return;

          CandidateType candidate(candidates.back());
          candidates.pop_back();

          if (DEBUG_MOLDB_ISOMORPHISM_H)
            std::cout << yellow << "candidate: " << candidate.first->index() << " -> " 
              << candidate.second->index() << normal << std::endl;

          // compute state s' obtained by adding (n, m) to state s
          updateState(state, candidate.first, candidate.second);
          
          // check feasibility function F(s)
          if (!ReducedGraph)
            if(!isStateFeasible(state))
              continue;

          if (ReducedGraph) {
            // check if the edges from queryVertex to already matched queried vertices match
            // alse check if the queryVertex self-loop and non-targeted edges match
            if (!checkTargetedEdges(state, candidate.first, candidate.second)) {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "    targeted edges do not match..." << std::endl;
              backtrack(state);
              continue;;
            }

            //if (!checkUntargetedEdges(state, candidate.first, candidate.second)) {
            //  if (DEBUG_MOLDB_ISOMORPHISM_H)
            //    std::cout << "    untargeted edges do not match..." << std::endl;
            //  state.usedQueriedEdges.clear();
            //  state.usedQueriedEdges.resize(state.queried->numEdges(), UsedEdge());
            //  backtrack(state);
            //  continue;
            //}
          } else {
            if (!checkEdges(state, candidate.first, candidate.second)) {
              if (DEBUG_MOLDB_ISOMORPHISM_H)
                std::cout << "    targeted edges do not match..." << std::endl;
              backtrack(state);
              continue;;
            }
          }


          // recurse
          isomorphismDFS<ReducedGraph>(state);
        }

        // backtrack
        if (!state.match)
          backtrack(state);
        */
      }

    } // namespace Isomorphism

  } // namespace Impl

  /**
   * @class GraphIsomorphismVF2 isomorphism.h
   *
   * Class for finding the isomorphism between two reduced graphs
   * (Graph). This class implements the VF2 algorithm.
   */
  template<bool ReducedGraph, typename GraphType>
  class GraphIsomorphismVF2
  {
    public:
      typedef typename GraphType::VertexType VertexType;
      typedef typename GraphType::EdgeType EdgeType;

      bool isSubgraph(const GraphType *query, const GraphType *queried)
      {
        // if query is empty, return true
        if (query->vertices().empty())
          return true;
        // if the queried graph is empty, there can be no match
        if (queried->vertices().empty())
          return false;


        /*
        for (std::size_t i = 0; i < query->vertices().size(); ++i) {
          bool vertexMatched = false;
          for (std::size_t j = 0; j < queried->vertices().size(); ++j)
            if (Impl::Isomorphism::verticesMatch(query->vertices()[i], queried->vertices()[j])) {
              vertexMatched = true;
              break;
            }
          if (!vertexMatched)
            return false;
        }
        */

        // create the state object
        Impl::Isomorphism::State<GraphType> state(query, queried);
        // start the recursive algorithm
        return Impl::Isomorphism::isomorphismDFS<ReducedGraph>(state);
      }
  };

}

#endif
