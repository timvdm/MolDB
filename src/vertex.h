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

#ifndef MOLDB_VERTEX_H
#define MOLDB_VERTEX_H

#include <vector>
#include <algorithm>

#include "memory.h"

namespace MolDB {

  namespace Impl {

    /**
     * @class GraphHelper vertex.h
     *
     * This class is a helper class for the Graph class to enable the graph
     * class to call Vertex's private addEdge function. This class is a friend
     * of the Vertex class.
     */
    template<typename VertexType, typename EdgeType>
    struct GraphHelper
    {
      /**
       * Add edge @p e to vertex @v.
       */
      void addEdge(VertexType *v, EdgeType *e)
      {
        v->addEdge(e);
      }

      void setIndex(VertexType *v, unsigned int index)
      {
        v->setIndex(index);
      }
      void setIndex(EdgeType *e, unsigned int index)
      {
        e->setIndex(index);
      }
    };

    template<typename GraphType>
    struct CompareEdgePathSize
    {
      typedef typename GraphType::IndexType IndexType;

      CompareEdgePathSize(GraphType *graph) : m_graph(graph)
      {
      }

      bool operator()(IndexType e1, IndexType e2)
      {
        return m_graph->edge(e1)->pathSize() > m_graph->edge(e2)->pathSize();
      }

      GraphType *m_graph;
    };

  }

  template<typename _GraphType>
  class Vertex
  {
    public:
      /**
       * The graph type.
       */
      typedef _GraphType GraphType;
      /**
       * The vertex type.
       */
      typedef typename GraphType::VertexType VertexType;
      typedef typename GraphType::EdgeType EdgeType;
      typedef typename GraphType::VertexEdgeSemanticsType SemanticsType;
      typedef typename GraphType::IndexType IndexType;
      typedef std::vector<IndexType> EdgesType;

      /**
       * Constructor.
       *
       * @param semantics The semantic attributes for this vertex.
       */
      Vertex(SemanticsType semantics) : m_semantics(semantics)
      {
      }

      /**
       * Get the semantic attributes for this vertex.
       *
       * @note Runtime complexity O(1).
       *
       * @return The semantic attributes for this vertex.
       */
      const SemanticsType& semantics() const
      {
        return m_semantics;
      }

      IndexType numEdges() const
      {
        return m_edges.size();
      }

      /**
       * Get the edges for this vertex.
       *
       * @note Runtime complexity O(1).
       *
       * @return The edges for this vertex.
       */
      const EdgesType& edges() const
      {
        return m_edges;
      }

      EdgeType* edge(const GraphType *graph, const VertexType *other) const
      {
        for (std::size_t i = 0; i < m_edges.size(); ++i) {
          if (graph->edge(m_edges[i])->source() == this && graph->edge(m_edges[i])->target() == other)
            return graph->edge(m_edges[i]);
          else if (graph->edge(m_edges[i])->source() == other && graph->edge(m_edges[i])->target() == this)
            return graph->edge(m_edges[i]);
        }
        
        return 0;
      }


      /**
       * Get the edges connecting this vertex to vertex @p other. This function
       * returns a list of edges since there may be parallel edges.
       *
       * @note Runtime complexity O(n) where n is the number of edges for this
       *       vertex.
       *
       * @param other The other vertex for the edges to which this vertex is
       *        connected.
       *
       * @return A std::vector containing the edges connecting this vertex to
       *         vertex @p other.
       */
      EdgesType edges(const GraphType *graph, const VertexType *other) const
      {
        EdgesType e;
        for (std::size_t i = 0; i < m_edges.size(); ++i) {
          if (graph->edge(m_edges[i])->source() == this && graph->edge(m_edges[i])->target() == other)
            e.push_back(m_edges[i]);
          else if (graph->edge(m_edges[i])->source() == other && graph->edge(m_edges[i])->target() == this)
            e.push_back(m_edges[i]);
        }
        
        return e;
      }

      /**
       * Get the memory size of this vertex object in bytes.
       *
       * @note Runtime complexity O(1).
       *
       * @return The memory size of this vertex object in bytes.
       */
      std::size_t memory_size() const
      {
        return sizeof(VertexType) + sizeof(EdgeType*) * m_edges.size();
      }

      IndexType index(GraphType *graph = 0) const
      {
        return m_index;
      }
 
      /**
       * Friend classes.
       */
      template<typename, typename> friend class Impl::GraphHelper; // FIXME
      template<typename, typename> friend class Graph;;

    private:
      /**
       * Add @p edge to this vertex.
       *
       * @note Runtime complexity O(n log(n)) to O(n^2) depending on the
       *       implementation of std::sort. Here, n is the number of edges
       *       for this vertex after adding edge @p edge.
       */
      void addEdge(GraphType *graph, EdgeType *edge)
      {
        if (std::find(m_edges.begin(), m_edges.end(), edge->index()) != m_edges.end())
          return;
        m_edges.push_back(edge->index());
        std::sort(m_edges.begin(), m_edges.end(), Impl::CompareEdgePathSize<GraphType>(graph));
      }

      void setIndex(IndexType index)
      {
        m_index = index;
      }

      /**
       * The vertex semantic attributes.
       */
      SemanticsType m_semantics;
      /**
       * The vertex index.
       */
      IndexType m_index;
      /**
       * The vertex edges.
       */
      EdgesType m_edges;
  };

}

#endif
