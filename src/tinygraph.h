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

#ifndef MOLDB_TINYGRAPH_H
#define MOLDB_TINYGRAPH_H

#include <vector>
#include <cassert>
#include <map>
#include <sstream>


#include <iostream> // FIXME

#include "typetraits.h"
#include "tinyvertex.h"
#include "tinyedge.h"
#include "tinyvector.h"
#include "exceptions.h"
#include "semantics.h"

namespace MolDB {


  template<typename T, typename UsageType>
  class TinyGraphVertexEdgeMemoryPool
  {
      typedef TinyGraphVertexEdgeMemoryPool<T, UsageType> MemoryPool;
    public:
      std::size_t size() const
      {
        return m_pool.size();
      }

      std::size_t capacity() const
      {
        return m_pool.capacity();
      }

      template<typename T1, typename T2, typename T3>
      T* allocate(const T1 &t1, const T2 &t2, const T3 &t3)
      {
        resizePool(t1, t2, t3);
        T *p = &m_pool[m_pool.size() - 1];
        return p;
      }

      void deallocate()
      {
        m_pool.clear();
      }

      static MemoryPool& instance(std::size_t size = 0)
      {
        static MemoryPool *p = 0;
        if (!p)
          p = new MemoryPool;
        if (size)
          p->m_pool.reserve(size);
        return *p;
      }

    private:
      TinyGraphVertexEdgeMemoryPool()
      {
      }

      template<typename T1, typename T2, typename T3>
      inline void resizePool(const T1 &t1, const T2 &t2, const T3 &t3, TypeTraits::Int2Type<true>)
      {
        m_pool.resize(m_pool.size() + 1, T(t1, t2));
      }

      template<typename T1, typename T2, typename T3>
      inline void resizePool(const T1 t1, const T2 t2, const T3 &t3, TypeTraits::Int2Type<false>)
      {
        m_pool.resize(m_pool.size() + 1, T(t1, t2, t3));
      }

      template<typename T1, typename T2, typename T3>
      inline void resizePool(const T1 &t1, const T2 &t2, const T3 &t3)
      {
        resizePool(t1, t2, t3, TypeTraits::Int2Type<TypeTraits::IsSame<T, typename T::GraphType>::result || TypeTraits::IsSame<T, typename T::VertexType>::result>());
      }

      std::vector<T> m_pool;
  };

  template<typename _GraphSemanticsType, typename _VertexEdgeSemanticsType, typename _IndexType>
  class TinyGraph
  {
    public:
      enum TinyUsages {
        GraphUsage,
        VertexUsage,
        EdgeUsage,
        GraphVerticesUsage,
        GraphEdgesUsage,
        VertexEdgesUsage,
        EdgePathUsage
      };

      typedef TinyGraph<_GraphSemanticsType, _VertexEdgeSemanticsType, _IndexType> GraphType;
      typedef TinyVertex<GraphType> VertexType;
      typedef TinyEdge<GraphType> EdgeType;
      typedef TinyVector<VertexType*, TypeTraits::Int2Type<GraphVerticesUsage> > VerticesType;
      typedef TinyVector<EdgeType*, TypeTraits::Int2Type<GraphEdgesUsage> > EdgesType;
      typedef _GraphSemanticsType GraphSemanticsType;
      typedef _VertexEdgeSemanticsType VertexEdgeSemanticsType;
      typedef _IndexType IndexType;
      typedef TinyGraphVertexEdgeMemoryPool<VertexType, TypeTraits::Int2Type<VertexUsage> > VertexMemoryPool;
      typedef TinyGraphVertexEdgeMemoryPool<EdgeType, TypeTraits::Int2Type<EdgeUsage> > EdgeMemoryPool;

      TinyGraph(std::size_t numVertices, std::size_t numEdges) 
          : m_vertices(numVertices), m_edges(numEdges)
      {
      }

      template<typename OtherGraphType>
      static GraphType* copy(OtherGraphType *other)
      {
        GraphType *g = new GraphType(other->numVertices(), other->numEdges());

        // copy the vertices
        for (std::size_t i = 0; i < other.numVertices(); ++i)
          g->addVertex(other->vertex(i)->semantics(), other->vertex(i)->edges());

        // copy edges
        for (std::size_t i = 0; i < other.numEdges(); ++i)
          g->addVertex(other->edges(i)->source(), other->edge(i)->target(), other->edge(i)->path());

        return g;
     }

      IndexType numVertices() const
      {
        return m_vertices.size();
      }

      const VerticesType& vertices() const
      {
        return m_vertices;
      }

      VertexType* vertex(IndexType index) const
      {
        return m_vertices[index];
      }

      VertexType* addVertex(VertexEdgeSemanticsType semantics, const std::vector<IndexType> &edges)
      {
        VertexType *vertex = VertexMemoryPool::instance().allocate(semantics, edges, 0);
        std::size_t index = m_vertices.size();
        vertex->setIndex(index);
        m_vertices.push_back(vertex);
        return m_vertices[index];
      }

      IndexType numEdges() const
      {
        return m_edges.size();
      }

      const std::vector<EdgeType*>& edges() const
      {
        return m_edges;
      }

      EdgeType* edge(IndexType index) const
      {
        return m_edges[index];
      }

      EdgeType* addEdge(VertexType *source, VertexType *target, const std::vector<VertexEdgeSemanticsType> &path)
      {
        EdgeType *edge = EdgeMemoryPool::instance().allocate(source, target, path);
        std::size_t index = m_edges.size();
        edge->setIndex(index);
        m_edges.push_back(edge);
        return m_edges[index];
      }

      std::size_t memory_size() const
      {
        std::size_t size = sizeof(GraphType);

        for (std::size_t i = 0; i < m_vertices.size(); ++i)
          size += m_vertices[i]->memory_size();
        
        for (std::size_t i = 0; i < m_edges.size(); ++i)
          size += m_edges[i]->memory_size();

        return size;
      }

      std::string toString(const std::string &ws = std::string(), const std::string &el = std::string())
      {
        // store vertex and edge indices in maps
        std::map<VertexType*, std::size_t> vertexIndices;
        std::map<EdgeType*, std::size_t> edgeIndices;

        for (std::size_t i = 0; i < m_vertices.size(); ++i)
          vertexIndices[m_vertices[i]] = i + 1;
        for (std::size_t i = 0; i < m_edges.size(); ++i)
          edgeIndices[m_edges[i]] = i + 1;

        // construct the string
        std::stringstream ss;
        ss << "{" << ws << "[" << ws << el;
        for (std::size_t i = 0; i < m_vertices.size(); ++i) {
          ss << ws << ws << ws << ws << "{" << ws << semanticsToString(m_vertices[i]->semantics());
          /*
          ss << "," << ws << "[" << ws;
          for (std::size_t j = 0; j < m_vertices[i]->edges().size(); ++j) {
            ss << edgeIndices[m_vertices[i]->edges()[j]];
            if (j + 1 < m_vertices[i]->edges().size())
              ss << "," << ws;
          }
          ss << ws << "]";
          */
          ss << ws << "}";
          if (i + 1 < m_vertices.size())
            ss << "," << el;
        }
        ss << el << "]" << ws << "," << ws << "[" << el;
        for (std::size_t i = 0; i < m_edges.size(); ++i) {
          ss << ws << ws << ws << ws << "{" << ws << vertexIndices[m_edges[i]->source()] << "," << ws;
          if (m_edges[i]->target())
            ss << vertexIndices[m_edges[i]->target()];
          else
            ss << 0;
          ss << ","<< ws << "[";
          for (std::size_t j = 0; j < m_edges[i]->path().size(); ++j) {
            ss << semanticsToString(m_edges[i]->path()[j]);  // FIXME
            if (j + 1 < m_edges[i]->path().size())
              ss << "," << ws;
          }
          ss << "]" << ws << "}";
          if (i + 1 < m_edges.size())
            ss << "," << el;
        }
        ss << el << "]" << ws << "}" << el;

        return ss.str();
      }
      
      template<typename, typename> friend class TinyGraphVertexEdgeMemoryPool;

    private:
      VerticesType m_vertices;
      EdgesType m_edges;
      GraphSemanticsType m_semantics;
  };

}

#endif
