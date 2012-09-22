#ifndef MOLDB_OPENBABEL_H
#define MOLDB_OPENBABEL_H

#include <string>
#include <fstream>
#include <bitset>
#include <cassert>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/fingerprint.h>

#include "exceptions.h"
#include "graph.h"

#define DEBUG_MOLDB_OPENBABEL_H 0

namespace MolDB {

  using namespace OpenBabel;
  
  enum PathTypes {
    // 0 - 199 reserved for element descriptors (atomic numbers and aromatic atoms)
    // single bonds are implicit
    AromaticBond = 200,
    DoubleBond = 201,
    TripleBond = 202
  };


  namespace Impl {

    inline int elementDescriptor(int atomicNumber, bool aromatic)
    {
      if (aromatic)
        switch (atomicNumber) {
          case 6: // C
            return 128;
          case 7: // N
            return 129;
          case 8: // O
            return 130;
          case 15: // P
            return 131;
          case 16: // S
            return 132;
          case 33: // As
            return 133;
          case 34: // Se
            return 134;
          default:
            break;
        }

      return atomicNumber;
    }


    template<typename VertexType>
    VertexType* createVertex(OpenBabel::OBAtom *atom)
    {
      VertexType *v = new VertexType(elementDescriptor(atom->GetAtomicNum(), atom->IsAromatic()));
      return v;
    }

    template<typename GraphType>
    void graphDFS(OpenBabel::OBAtom *from, OpenBabel::OBAtom *to, 
        typename GraphType::VertexType *fromVertex, 
        std::vector<bool> &visitedAtoms, std::vector<bool> &visitedBonds,
        std::map<OpenBabel::OBAtom*, typename GraphType::VertexType*> &vertices, 
        GraphType *g, typename GraphType::EdgeType::PathType &path)
    {
      typedef typename GraphType::VertexType VertexType;
      typedef typename GraphType::EdgeType EdgeType;
      typedef typename EdgeType::PathType PathType;

      //std::cout << "graphDFS(" << to->GetIndex() << ")" << std::endl;


      // mark atom as visited
      visitedAtoms[to->GetIndex()] = true;


      // add bond to path if needed (i.e. not single)
      if (from) {
        OpenBabel::OBBond *bond = from->GetBond(to);
        visitedBonds[bond->GetIdx()] = true;
        if (bond->IsAromatic())
          path.push_back(AromaticBond);
        else if (bond->IsDouble())
          path.push_back(DoubleBond);
        else if (bond->IsTriple())
          path.push_back(TripleBond);
      }

      if (to->GetValence() == 1 && fromVertex) {
        // we've reached the end of a path that has no target vertex
        path.push_back(elementDescriptor(to->GetAtomicNum(), to->IsAromatic()));
        EdgeType *edge = new EdgeType(fromVertex, 0, path);
        if (DEBUG_MOLDB_OPENBABEL_H)
          std::cout << "    addEdge_1(" << fromVertex->index() + 1 << ", 0)" << std::endl;
        g->addEdge(edge);
        path.clear();
      } else if (to->GetValence() == 2) {
        // close rings where there are no vertices with degree > 2
        if (vertices.find(to) != vertices.end()) {
          EdgeType *edge = new EdgeType(fromVertex, vertices[to], path);
          if (DEBUG_MOLDB_OPENBABEL_H)
            std::cout << "    addEdge_2(" << fromVertex->index() + 1 << ", " << vertices[to]->index() + 1 << ")" << std::endl;
          g->addEdge(edge);
          return;
        }

        // continue the path (until we reach a new vertex or end)
        path.push_back(elementDescriptor(to->GetAtomicNum(), to->IsAromatic()));
        FOR_NBORS_OF_ATOM (nbr, to)
          if (from != &*nbr) {
            graphDFS(to, &*nbr, fromVertex, visitedAtoms, visitedBonds, vertices, g, path);
            break;
          }
      } else if (to->GetValence() > 2) {
        // to->GetValence() > 2

        // is the vertex already constructed? if so create an Edge and return
        if (vertices.find(to) != vertices.end()) {
          EdgeType *edge = new EdgeType(fromVertex, vertices[to], path);
          if (DEBUG_MOLDB_OPENBABEL_H)
            std::cout << "    addEdge_3(" << fromVertex->index() + 1 << ", " << vertices[to]->index() + 1 << ")" << std::endl;
          g->addEdge(edge);
          return;
        }

        // create new vertex in the graph
        VertexType *v = createVertex<VertexType>(to);
        vertices[to] = v;
        g->addVertex(v);
        if (fromVertex) {
          EdgeType *edge = new EdgeType(fromVertex, v, path);
          if (DEBUG_MOLDB_OPENBABEL_H)
            std::cout << "    addEdge_4(" << fromVertex->index() + 1 << ", " << v->index() + 1 << ")" << std::endl;
          g->addEdge(edge);
        }

        FOR_NBORS_OF_ATOM (nbr, to) {
          if (visitedAtoms[nbr->GetIndex()]) {
            if (fromVertex) {
              OpenBabel::OBBond *bond = nbr->GetBond(to);
              if (!visitedBonds[bond->GetIdx()]) {
                visitedBonds[bond->GetIdx()] = true;
                PathType path;
                if (bond->IsAromatic())
                  path.push_back(AromaticBond);
                else if (bond->IsDouble())
                  path.push_back(DoubleBond);
                else if (bond->IsTriple())
                  path.push_back(TripleBond);
                EdgeType *edge = new EdgeType(v, vertices[&*nbr], path);
                g->addEdge(edge);
                if (DEBUG_MOLDB_OPENBABEL_H)
                  std::cout << "    addEdge_5(" << v->index() + 1 << ", " << vertices[&*nbr]->index() + 1 << ")" << std::endl;
              }
            }
            continue;
          }
          // start a new path
          PathType newPath;
          graphDFS(to, &*nbr, v, visitedAtoms, visitedBonds, vertices, g, newPath);
        }
      }
    }

  }


  struct OB
  {
      
    class File
    {
      public:
        File(const std::string &filename)
        {
          m_ifstream.open(filename.c_str());

          if (!m_ifstream)
            throw Exceptions::CouldNotOpenFile();

          OpenBabel::OBFormat *format = m_conv.FormatFromExt(filename);

          if (!format || !m_conv.SetInFormat(format))
            throw Exceptions::CouldNotFindFormat();
        }

        OpenBabel::OBMol* next()
        {
          assert(m_ifstream);
          OpenBabel::OBMol *mol = new OpenBabel::OBMol;

          if (!m_conv.Read(mol, &m_ifstream)) {
            delete mol;
            return 0;
          }

          return mol;
        }

        std::ios_base::streampos tellg()
        {
          return m_ifstream.tellg();
        }

        void seekg(std::ios_base::streampos pos)
        {
          m_ifstream.seekg(pos);
        }

        static OpenBabel::OBMol* readFile(const std::string &filename)
        {
          File file(filename);
          return file.next();
        }

        static OpenBabel::OBMol* readSmiles(const std::string &smiles)
        {
          OpenBabel::OBConversion conv;
          if (!conv.SetInFormat("smi"))
            throw Exceptions::CouldNotFindFormat();

          OpenBabel::OBMol *mol = new OpenBabel::OBMol;
          if (!conv.ReadString(mol, smiles)) {
            delete mol;
            return 0;
          }

          return mol;
        }

        static bool isValidFileFormat(const std::string &filename)
        {
          OpenBabel::OBConversion conv;
          return conv.FormatFromExt(filename);
        }

        static std::string writeSmiles(OpenBabel::OBMol *mol)
        {
          OpenBabel::OBConversion conv;
          if (!conv.SetOutFormat("can"))
            throw Exceptions::CouldNotFindFormat();

         std::string smiles = conv.WriteString(mol, true);
         if (smiles.find(" ") != std::string::npos)
           smiles = smiles.substr(0, smiles.find(" "));
         if (smiles.find("\t") != std::string::npos)
           smiles = smiles.substr(0, smiles.find("\t"));
         
         return smiles;
        }


      private:
        std::ifstream m_ifstream;
        OpenBabel::OBConversion m_conv;
    };

    template<int Size>
    static std::bitset<Size> fingerprint(OpenBabel::OBMol *mol, const std::string &type)
    {
      OpenBabel::OBFingerprint *fp = OpenBabel::OBFingerprint::FindFingerprint(type.c_str());

      if (!fp)
        throw Exceptions::CouldNotFindFingerprint();

      std::vector<unsigned int> fpbits;
      fp->GetFingerprint(mol, fpbits, Size);

      std::bitset<Size> fpbitset;
      for (std::size_t i = 0; i < fpbits.size(); ++i)
        fpbitset = fpbitset | (std::bitset<Size>(fpbits[i]) << (8 * sizeof(unsigned int) * i));

      return fpbitset;
    }

    template<typename GraphType>
    static GraphType* reducedGraph(OpenBabel::OBMol *mol)
    {
      typedef typename GraphType::VertexType VertexType;
      typedef typename GraphType::EdgeType::PathType PathType;

      GraphType *g = new GraphType;

      if (!mol->NumAtoms())
        return g;

      // find a node atom (has at least 3 bonds)
      OpenBabel::OBAtom *node = 0;
      FOR_ATOMS_OF_MOL (atom, mol)
        if (atom->GetValence() > 2) {
          node = &*atom;
          break;
        }
      
      VertexType *fromVertex = 0;
      std::map<OpenBabel::OBAtom*, VertexType*> vertices;
      

      // special case: linear molecules with no degree(atom) > 2
      if (!node) {
        FOR_ATOMS_OF_MOL (atom, mol)
          if (atom->GetValence() == 1) {
            node = &*atom;
            break;
          }
        
        if (!node)
          node = mol->GetAtom(1);

        fromVertex = Impl::createVertex<VertexType>(node);
        vertices[node] = fromVertex;
      }

      PathType path;
      std::vector<bool> visitedAtoms(mol->NumAtoms(), false);
      std::vector<bool> visitedBonds(mol->NumBonds(), false);

      Impl::graphDFS(static_cast<OpenBabel::OBAtom*>(0), node, fromVertex, visitedAtoms, visitedBonds, vertices, g, path);
 
      return g;
    }

    template<typename GraphType>
    static GraphType* graph(OpenBabel::OBMol *mol)
    {
      typedef typename GraphType::VertexType VertexType;
      typedef typename GraphType::EdgeType EdgeType;
      typedef typename GraphType::EdgeType::PathType PathType;

      GraphType *g = new GraphType;

      if (!mol->NumAtoms())
        return g;

      FOR_ATOMS_OF_MOL (atom, mol) {
        VertexType *v = Impl::createVertex<VertexType>(&*atom);
        g->addVertex(v);
      }

      FOR_BONDS_OF_MOL (bond, mol) {
        PathType path;
        if (bond->IsAromatic())
          path.push_back(AromaticBond);
        else if (bond->IsDouble())
          path.push_back(DoubleBond);
        else if (bond->IsTriple())
          path.push_back(TripleBond);

        EdgeType *edge = new EdgeType(g->vertices()[bond->GetBeginAtom()->GetIndex()], g->vertices()[bond->GetEndAtom()->GetIndex()], path);
        g->addEdge(edge);
      }
 
      return g;
    }


    template<typename GraphType>
    static GraphType* reducedGraph(const std::string &smiles)
    {
      OpenBabel::OBMol *mol = File::readSmiles(smiles);
      GraphType *g = reducedGraph<GraphType>(mol);
      delete mol;
      return g;
    }
 
  };

}

#endif
