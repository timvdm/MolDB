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

#ifndef MOLDB_GRAPH_H
#define MOLDB_GRAPH_H

#include <vector>
#include <cassert>
#include <map>
#include <sstream>


#include <iostream> // FIXME

#include "vertex.h"
#include "edge.h"
#include "exceptions.h"
#include "semantics.h"

namespace MolDB {

  namespace Impl {
  
    // findEnd(3, "   { { { } }, { }, }", "{", "}")
    inline std::size_t findEnd(std::size_t pos, const std::string &s, const std::string &open, const std::string &close)
    {
      int indent = 1;
      std::size_t currpos = pos;

      while (indent) {
        std::size_t openpos = s.find(open, currpos);
        std::size_t closepos = s.find(close, currpos);

        if (closepos == std::string::npos) {
          std::stringstream ss;
          ss << "No matching '" << close << "' for '" << open << "' at position " << pos << " in string: \"" << s << "\"";
          throw Exceptions::GraphStringSyntaxError(ss.str(), s, pos);
        }

        if (openpos == std::string::npos) {
          indent--;
          currpos = closepos + 1;
        } else if (openpos > closepos) {
          indent--;
          currpos = closepos + 1;
        } else {
          indent++;
          currpos = openpos + 1;
        }
      }

      return currpos - 1;
    }

    inline std::vector<std::string> tokenize(const std::string &s, const std::string &delimiter)
    {
      std::vector<std::string> tokens;
      std::size_t currpos = 0, nextpos = 0;

      while ((nextpos = s.find(delimiter, currpos)) != std::string::npos) {
        tokens.push_back(s.substr(currpos, nextpos - currpos));
        currpos = nextpos + 1;
      }
      tokens.push_back(s.substr(currpos, s.length() - currpos));

      return tokens;
    }

  }

  template<typename _GraphSemanticsType, typename _VertexEdgeSemanticsType>
  class Graph
  {
    public:
      typedef Graph<_GraphSemanticsType, _VertexEdgeSemanticsType> GraphType;
      typedef Vertex<GraphType> VertexType;
      typedef Edge<GraphType> EdgeType;
      typedef _GraphSemanticsType GraphSemanticsType;
      typedef _VertexEdgeSemanticsType VertexEdgeSemanticsType;
      typedef unsigned int IndexType;
      typedef std::vector<VertexType*> VerticesType;
      typedef std::vector<EdgeType*> EdgesType;

      Graph()
      {
      }

      ~Graph()
      {
        for (std::size_t i = 0; i < m_vertices.size(); ++i)
          delete m_vertices[i];
        
        for (std::size_t i = 0; i < m_edges.size(); ++i)
          delete m_edges[i];
      }

      IndexType numVertices() const
      {
        return m_vertices.size();
      }

      const std::vector<VertexType*>& vertices() const
      {
        return m_vertices;
      }

      VertexType* vertex(IndexType index) const
      {
        return m_vertices[index];
      }

      void addVertex(VertexType *vertex)
      {
        Impl::GraphHelper<VertexType, EdgeType>().setIndex(vertex, m_vertices.size());
        m_vertices.push_back(vertex);
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

      void addEdge(EdgeType *edge)
      {
        edge->setIndex(m_edges.size());
        assert(edge->source());
        m_edges.push_back(edge);
        edge->source()->addEdge(this, edge);
        if (edge->target())
          edge->target()->addEdge(this, edge);
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

      // { [
      //     { 1, [ 1, 2, 3 ] }
      // ] , [
      //     { 1, 2, [ 6, 6 ] }
      // ] }

      static GraphType* fromString(const std::string &s)
      {
        GraphType *g = new GraphType;
        std::size_t pos;

        std::size_t startpos = s.find("{");
        if (pos == std::string::npos) {
          std::stringstream ss;
          ss << "Cannot find opening '{' for graph object in string: " << s;
          throw Exceptions::GraphStringSyntaxError(ss.str(), s, 0);
        }

        std::size_t endpos = Impl::findEnd(startpos + 1, s, "{", "}");

        std::size_t atomsStartPos = s.find("[");
        if (atomsStartPos != std::string::npos) {
          std::size_t atomsEndPos = Impl::findEnd(atomsStartPos + 1, s, "[", "]");
          std::size_t atomStartPos = atomsStartPos + 1;

          //std::cout << "atoms = " << s.substr(atomsStartPos, atomsEndPos - atomsStartPos + 1) << std::endl;

          while ((atomStartPos = s.find("{", atomStartPos)) < atomsEndPos) {
 
            std::size_t atomEndPos = Impl::findEnd(atomStartPos + 1, s, "{", "}");
            std::string element = s.substr(atomStartPos + 1, atomEndPos - atomStartPos - 1);
            //std::cout << "atom = " << s.substr(atomStartPos, atomEndPos - atomStartPos + 1) << std::endl;
            /*
            std::size_t delimiterPos = s.find(",", atomStartPos + 1);
            if (delimiterPos == std::string::npos || delimiterPos > atomEndPos) {
              std::stringstream ss;
              ss << "Cannot find ',' seperating element from edge list in vertex ";
              ss << s.substr(atomStartPos, s.length() - atomEndPos - 1);
              ss << " starting at position " << atomStartPos << " in string: " << s;
              throw GraphStringSyntaxError(ss.str(), s, atomStartPos);
            }

            std::string element = s.substr(atomStartPos + 1, atomEndPos - delimiterPos - 1);
            
            std::size_t edgesStartPos = s.find("[");
            if (edgeStartPos == std::string::npos || edgeStartPos > atomEndPos) {
              std::stringstream ss;
              ss << "Cannot find '[' starting the edge list in vertex ";
              ss << s.substr(atomStartPos, s.length() - atomEndPos - 1);
              ss << " starting at position " << atomStartPos << " in string: " << s;
              throw GraphStringSyntaxError(ss.str(), s, atomStartPos);
            }
            
            std::size_t edgesEndPos = Impl::findEnd(edgesStartPos + 1, s, "[", "]");
            if (edgesEndPos > atomEndPos) {
              std::stringstream ss;
              ss << "Cannot find ']' ending the edge list in vertex ";
              ss << s.substr(atomStartPos, s.length() - atomEndPos - 1);
              ss << " starting at position " << atomStartPos << " in string: " << s;
              throw GraphStringSyntaxError(ss.str(), s, atomStartPos);
            }

            std::vector<std::string> edges = Impl::tokenize(s.substr(edgesStartPos, s.length() - edgesEndPos - 1, ","));
            */

            std::stringstream ss(element);
            int elementInt;
            try {
              ss >> elementInt;
            } catch (...) {
              ss.str() = "";
              ss << "Cannot interpret element \"" << element << "\"";
              throw Exceptions::GraphStringSyntaxError(ss.str(), element, 0);
            }

            g->addVertex(new VertexType(elementInt));
            atomStartPos = atomEndPos + 1;
          }

          std::size_t delimiterVertexEdge = s.find(",", atomsEndPos);
          // todo handle npos
          std::size_t edgesStartPos = s.find("[", atomsEndPos);
          // handle npos & delimiterVertexEdge > edgesStartPos

          std::size_t edgesEndPos = Impl::findEnd(edgesStartPos + 1, s, "[", "]");
          //std::cout << "edges = " << s.substr(edgesStartPos, edgesEndPos - edgesStartPos + 1) << std::endl;
          std::size_t edgeStartPos = s.find("{", edgesStartPos + 1);
          while ((edgeStartPos = s.find("{", edgeStartPos)) < edgesEndPos) {
            std::size_t edgeEndPos = Impl::findEnd(edgeStartPos + 1, s, "{", "}");

            int source, target;
            std::size_t delimiter1 = s.find(",", edgeStartPos + 1);
            std::stringstream ss1(s.substr(edgeStartPos + 1, delimiter1 - edgeStartPos - 1));
            ss1 >> source;
            
            std::size_t delimiter2 = s.find(",", delimiter1 + 1);
            std::stringstream ss2(s.substr(delimiter1 + 1, delimiter2 - delimiter1 - 1));
            ss2 >> target;

            //std::cout << "edge = " << s.substr(edgeStartPos, edgeEndPos - edgeStartPos + 1) << std::endl;

            typename EdgeType::PathType path;
            std::size_t vertexStartPos = s.find("[", delimiter2 + 1);
            std::size_t vertexEndPos = s.find("]", vertexStartPos + 1);

            std::vector<std::string> tokens = Impl::tokenize(s.substr(vertexStartPos + 1, vertexEndPos - vertexStartPos - 1), ",");

            for (std::size_t i = 0; i < tokens.size(); ++i) {
              int p;
              std::stringstream ss(tokens[i]);
              ss >> p;
              path.push_back(p);
            }

            if (target)
              g->addEdge(new EdgeType(g->vertices()[source - 1], g->vertices()[target - 1], path));
            else
              g->addEdge(new EdgeType(g->vertices()[source - 1], 0, path));

            edgeStartPos = edgeEndPos + 1;
          }


        }

        return g;
      }





    private:
      std::vector<VertexType*> m_vertices;
      std::vector<EdgeType*> m_edges;
      GraphSemanticsType m_semantics;
  };

}

#endif
