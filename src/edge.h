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

#ifndef MOLDB_EDGE_H
#define MOLDB_EDGE_H

#include <vector>

#include "typetraits.h"
#include "semantics.h"

namespace MolDB {

  template<typename _GraphType>
  class Edge
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
      /**
       * The edge type.
       */
      typedef typename GraphType::EdgeType EdgeType;
      /**
       * The semantics type.
       */
      typedef typename GraphType::VertexEdgeSemanticsType SemanticsType;
      /**
       * The path type. If the semantics type (SemanticsType) is an integer
       * type this type will be std::vector<SemanticsType>. In all other cases
       * this type becomes std::vector<SemanticsType*> to avoid costly copying.
       */
      typedef std::vector<typename 
          TypeTraits::Select<TypeTraits::IsInteger<SemanticsType>::result, 
          SemanticsType, SemanticsType*>::result> PathType;

      typedef typename GraphType::IndexType IndexType;

      /**
       * Constructor.
       *
       * @param source The source vertex.
       * @param target The target vertex.
       * @param path The path containing the semantic attributes for this edge.
       */
      Edge(VertexType *source, VertexType *target, const PathType &path)
          : m_source(source), m_target(target), m_path(path)
      {
      }

      /**
       * Get the source vertex for this edge.
       *
       * @note Runtime complexity O(1).
       *
       * @return The source vertex for this edge.
       */
      VertexType* source() const
      {
        return m_source;
      }

      /**
       * Get the target vertex for this edge.
       *
       * @note Runtime complexity O(1).
       *
       * @return The target vertex for this edge.
       */
      VertexType* target() const
      {
        return m_target;
      }

      /**
       * Get the other vertex given vertex @p v. If @p v matches this edge's
       * source vertex, the target vertex is returned. Otherwise the source
       * vertex is returned.
       *
       * @note Runtime complexity O(1).
       *
       * @return The other vertex for this edge.
       */
      VertexType* other(VertexType *v)
      {
        return (v == m_source) ? m_target : m_source;
      }

      IndexType pathSize() const
      {
        return m_path.size();
      }

      /**
       * Get the path of semantic attributes for this vertex.
       *
       * @note Runtime complexity O(1).
       *
       * @return A std::vector containing pointers to the semantic attributes.
       */
      const PathType& path() const
      {
        return m_path;
      }

      /**
       * Get the memory size for this edge object in bytes.
       *
       * @note Runtime complexity O(1).
       *
       * @return The memory size for this edge object in bytes.
       */
      std::size_t memory_size() const
      {
        return sizeof(EdgeType) + sizeof(SemanticsType) * m_path.size();
      }

      std::string toString(const GraphType *g, const std::string &ws = std::string(), const std::string &el = std::string())
      {
        // store vertex and edge indices in maps
        std::map<VertexType*, std::size_t> vertexIndices;

        for (std::size_t i = 0; i < g->vertices().size(); ++i)
          vertexIndices[g->vertices()[i]] = i + 1;

        // construct the string
        std::stringstream ss;
        ss << "{" << ws << vertexIndices[m_source] << "," << ws;
        if (m_target)
          ss << vertexIndices[m_target];
        else
          ss << 0;
        ss << ","<< ws << "[";
        for (std::size_t j = 0; j < m_path.size(); ++j) {
          ss << semanticsToString(m_path[j]);
          if (j + 1 < m_path.size())
            ss << "," << ws;
        }
        ss << "]" << ws << "}";
        
        return ss.str();
      }

      IndexType index() const
      {
        return m_index;
      }
      
      template<typename, typename> friend class Graph;

    private:
      void setIndex(IndexType index)
      {
        m_index = index;
      }

      /**
       * The source vertex.
       */
      VertexType *m_source;
      /**
       * The target vertex.
       */
      VertexType *m_target;
      /**
       * The edge index.
       */
      IndexType m_index;
      /**
       * The path containing semantic attributes.
       */
      PathType m_path;
  };

}

#endif
