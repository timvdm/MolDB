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

#ifndef MOLDB_TINYVERTEX_H
#define MOLDB_TINYVERTEX_H

#include <vector>
#include <algorithm>

#include "tinyvector.h"
#include "typetraits.h"

namespace MolDB {

  namespace Impl {

    template<typename EdgeType>
    inline bool compareTinyEdgePathLength(const EdgeType *e, const EdgeType *f)
    {
      return e->path().size() > f->path().size();
    }

  }

  template<typename _GraphType>
  class TinyVertex
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
      typedef TinyVector<IndexType, TypeTraits::Int2Type<GraphType::VertexEdgesUsage> > EdgesType;

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
        return m_numEdges;
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

      EdgeType* edge(IndexType index) const
      {
        return m_edges[index];
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
      std::vector<EdgeType*> edges(VertexType *other)
      {
        std::vector<EdgeType*> e;
        for (std::size_t i = 0; i < m_edges.size(); ++i) {
          if (m_edges[i]->source() == this && m_edges[i]->target() == other)
            e.push_back(m_edges[i]);
          else if (m_edges[i]->source() == other && m_edges[i]->target() == this)
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

      IndexType index() const
      {
        return m_index;
      }
 
      /**
       * Friend classes.
       */
      template<typename, typename, typename> friend class TinyGraph;
      template<typename, typename> friend class TinyGraphVertexEdgeMemoryPool;

    private:
      /**
       * Constructor.
       *
       * @param semantics The semantic attributes for this vertex.
       */
      TinyVertex(SemanticsType semantics, const std::vector<IndexType> &edges) : m_semantics(semantics), m_edges(edges.size())
      {
        for (std::size_t i = 0; i < edges.size(); ++i)
          m_edges.push_back(edges[i]);
        m_numEdges = edges.size();
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
       * The number of edges.
       */
      IndexType m_numEdges;
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
