#include "database.h"
#include "openbabel.h"
#include "stringvector.h"
#include "timer.h"
#include "util.h"
#include "isomorphism.h"
#include "messagehandler.h"

#include "smilesindex.h"
#include "fingerprintindex.h"

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

// OpenBabel
#include <openbabel/fingerprint.h>
#include <openbabel/obconversion.h>
// ChemKit
#ifdef HAVE_CHEMKIT
#include <chemkit/molecule.h>
#include <chemkit/fingerprint.h>
#include <chemkit/substructurequery.h>
#endif
// RDKit
#ifdef HAVE_RDKIT
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#endif

#ifdef _MSC_VER
#define FUNCTION_SIGNATURE __FUNCSIG__
#else
#define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#endif


namespace MolDB {

  typedef Graph<TypeTraits::EmptyType, unsigned char> GraphType;
  typedef GraphType::VertexType VertexType;
  typedef GraphType::EdgeType EdgeType;

  typedef TinyGraph<TypeTraits::EmptyType, unsigned char, unsigned short> TinyGraphType;
  typedef TinyGraphType::VertexType TinyVertexType;
  typedef TinyGraphType::EdgeType TinyEdgeType;

  typedef TinyGraphVertexEdgeMemoryPool<TinyGraphType, TypeTraits::Int2Type<TinyGraphType::GraphUsage> > TinyGraphMemoryPool;
  typedef Impl::TinyVectorMemoryPool<TinyVertexType*, TypeTraits::Int2Type<TinyGraphType::GraphVerticesUsage> > TinyGraphVerticesMemoryPool;
  typedef Impl::TinyVectorMemoryPool<TinyEdgeType*, TypeTraits::Int2Type<TinyGraphType::GraphEdgesUsage> > TinyGraphEdgesMemoryPool;
  typedef TinyGraphVertexEdgeMemoryPool<TinyVertexType, TypeTraits::Int2Type<TinyGraphType::VertexUsage> > TinyVertexMemoryPool;
  typedef Impl::TinyVectorMemoryPool<TinyGraphType::IndexType, TypeTraits::Int2Type<TinyGraphType::VertexEdgesUsage> > TinyVertexEdgesMemoryPool;
  typedef TinyGraphVertexEdgeMemoryPool<TinyEdgeType, TypeTraits::Int2Type<TinyGraphType::EdgeUsage> > TinyEdgeMemoryPool;
  typedef Impl::TinyVectorMemoryPool<TinyGraphType::VertexEdgeSemanticsType, TypeTraits::Int2Type<TinyGraphType::EdgePathUsage> > TinyEdgePathMemoryPool;

  // fingerprint storage
  typedef InMemoryFingerprintStorage FingerprintStorageType;
  //typedef MemoryMappedFingerprintStorage FingerprintStorageType;

  // fingerprint index
  typedef BruteForceFingerprintIndex<FingerprintStorageType> FingerprintIndexType;
  //typedef InvertedFingerprintIndex<FingerprintStorageType> FingerprintIndexType;
  //ypedef BitSelectionBinarySearchTreeFingerprintIndex<FingerprintStorageType, 14, 1> FingerprintIndexType;
  //typedef KDBinarySearchTreeFingerprintIndex<FingerprintStorageType, 5> FingerprintIndexType;
  //typedef CompressedBinarySearchTreeFingerprintIndex<FingerprintStorageType, false> FingerprintIndexType;
  //typedef KDGridFingerprintIndex<FingerprintStorageType, 1> FingerprintIndexType;
  //typedef MixedFingerprintIndex<FingerprintStorageType, BitSelectionBinarySearchTreeFingerprintIndex<FingerprintStorageType, 16, 1>, KDGridFingerprintIndex<FingerprintStorageType, 3> > FingerprintIndexType;

  // smiles index
  typedef InMemoryMapSmilesIndex<> SmilesIndexType;


  struct DatabasePrivate
  {
    DatabasePrivate() : graphCount(0), vertexCount(0), edgeCount(0),
        graphVerticesCount(0), graphEdgesCount(0),vertexEdgesCount(0), 
        edgePathCount(0), isCreated(false), isLoaded(false)
    {
    }

    ~DatabasePrivate()
    {
      delete smilesIndex;
      for (std::size_t i = 0; i < fingerprintIndexes.size(); ++i) {
        delete fingerprintStorages[i];
        delete fingerprintIndexes[i];
      }
    }

    /**
     * The database filename.
     */
    std::string filename;
    /**
     * The custom data functor.
     */
    Database::CustomDataFunctor *functor;
    /**
     * The message handler.
     */
    MessageHandler *msg;
    /**
     * True if the message handler is owned by the database.
     */
    bool ownedMessageHandler;
    /**
     * The maximum number of graphs in memory.
     */
    unsigned long maxGraphsInMemory;
    std::ios_base::streampos endOfFile;
    bool isCreated;
    bool isLoaded;


    std::vector<std::string> fingerprintNames;
    std::vector<unsigned long> positions;

    std::vector<FingerprintStorageType*> fingerprintStorages;
    std::vector<FingerprintIndex*> fingerprintIndexes;

    SmilesIndexType *smilesIndex;

    std::vector<TinyGraphType*> graphs;
 
    unsigned long graphCount;
    unsigned long vertexCount;
    unsigned long edgeCount;
    unsigned long graphVerticesCount;
    unsigned long graphEdgesCount;
    unsigned long vertexEdgesCount;
    unsigned long edgePathCount;

    boost::mutex mutex;
  };

  namespace Impl {
    
#ifdef HAVE_CHEMKIT
    chemkit::Molecule* chemKitMolecule(OpenBabel::OBMol *mol)
    {
      chemkit::Molecule *molecule = new chemkit::Molecule;
      std::map<OpenBabel::OBAtom*, chemkit::Atom*> atoms;

      FOR_ATOMS_OF_MOL (atom, mol)
        atoms[&*atom] = molecule->addAtom(atom->GetAtomicNum());

      FOR_BONDS_OF_MOL (bond, mol)
        molecule->addBond(atoms[bond->GetBeginAtom()], atoms[bond->GetEndAtom()], bond->GetBO());

      return molecule;
    }
#endif

    Fingerprint fingerprint(const std::string &fingerprintName, OpenBabel::OBMol *molecule);
  
    std::size_t fingerprintSize(const std::string &fingerprintName)
    {
      OpenBabel::OBMol *molecule = OB::File::readSmiles("CCC");
      std::size_t result = fingerprint(fingerprintName, molecule).size();
      return result;
    }


    Fingerprint fingerprint(const std::string &fingerprintName, OpenBabel::OBMol *molecule)
    {
      if (fingerprintName.find("MolDB::") != std::string::npos) {
        std::string filename;
        if (fingerprintName == "MolDB::MolDB")
          filename = "moldb.def";
        if (fingerprintName == "MolDB::Andrew")
          filename = "andrew.def";
        if (fingerprintName == "MolDB::SubStruct")
          filename = "substruct.def";

        FingerprintGenerator *generator = FingerprintGenerator::instance(filename);

        return generator->fingerprint(molecule);
      }
      
      if (fingerprintName.find("OpenBabel::") != std::string::npos) {
        OpenBabel::OBFingerprint *fp = OpenBabel::OBFingerprint::FindFingerprint(fingerprintName.substr(11, fingerprintName.size() - 11).c_str());
        if (!fp)
          return Fingerprint();

        std::vector<unsigned int> fpbits;
        fp->GetFingerprint(molecule, fpbits);

        Fingerprint::size_type numBits = fpbits.size() * fp->Getbitsperint();
        Fingerprint fpbitset(numBits);
        for (std::size_t i = 0; i < fpbits.size(); ++i)
          fpbitset |= (Fingerprint(numBits, fpbits[i]) << (fp->Getbitsperint() * i));

        return fpbitset;
      }

#ifdef HAVE_CHEMKIT
      if (fingerprintName.find("ChemKit::") != std::string::npos) {
        chemkit::Fingerprint *fp = chemkit::Fingerprint::create(fingerprintName.substr(9, fingerprintName.size() - 9));
        assert(fp);
        if (!fp)
          return Fingerprint();
        chemkit::Molecule chemkitMolecule(OB::File::writeSmiles(molecule), "smiles");
        Fingerprint fpbitset = fp->value(&chemkitMolecule);

        delete fp;
        return fpbitset;
      }
#endif

#ifdef HAVE_RDKIT
      if (fingerprintName.find("RDKit::") != std::string::npos) {
        RDKit::RWMol *rdkitMolecule = 0;
        try {
          rdkitMolecule = RDKit::SmilesToMol(OB::File::writeSmiles(molecule));
        } catch (...) {
          return Fingerprint(fingerprintSize(fingerprintName)).flip(); 
        }

        if (fingerprintName == "RDKit::Topological") {
          ExplicitBitVect *fpbits = RDKit::RDKFingerprintMol(*rdkitMolecule);
          Fingerprint fpbitset = *fpbits->dp_bits;
          delete fpbits;
          delete rdkitMolecule;
          return fpbitset;
        }
        
        if (fingerprintName == "RDKit::AtomPairs") {
          ExplicitBitVect *fpbits = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(*rdkitMolecule);
          Fingerprint fpbitset = *fpbits->dp_bits;
          delete rdkitMolecule;
          delete fpbits;
          return fpbitset;
        }

        if (fingerprintName == "RDKit::TopologicalTorsions") {
          ExplicitBitVect *fpbits = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(*rdkitMolecule);
          Fingerprint fpbitset = *fpbits->dp_bits;
          delete rdkitMolecule;
          delete fpbits;
          return fpbitset;
        }

        if (fingerprintName == "RDKit::Layered") {
          ExplicitBitVect *fpbits = RDKit::LayeredFingerprintMol(*rdkitMolecule);
          Fingerprint fpbitset = *fpbits->dp_bits;
          delete rdkitMolecule;
          delete fpbits;
          return fpbitset;
        }  
      }
#endif

      return Fingerprint();
    }

    struct DatabaseHeader
    {
      unsigned int magicNumber;
      unsigned long endOfFile;
      unsigned long graphCount;
      unsigned long vertexCount;
      unsigned long edgeCount;
      unsigned long graphVerticesCount;
      unsigned long graphEdgesCount;
      unsigned long vertexEdgesCount;
      unsigned long edgePathCount;
    };
 

  }


  struct DatabaseImpl
  {
    /**
     * Get the database filename (including path) withouth the .moldb extension.
     */
    static std::string baseFilename(const std::string &filename)
    {
      return filename.substr(0, filename.find(".moldb"));
    }
   
    /**
     * Set various data members in the d data pointer.
     */
    static void construct(DatabasePrivate *d, const std::string &filename, Database::CustomDataFunctor *functor,
        MessageHandler *messageHandler, bool ownedMessageHandler)
    {
      // make sure OpenBabel plugins are loaded
      OpenBabel::OBConversion();
      d->maxGraphsInMemory = 1000000;
      d->filename = filename;
      d->functor = functor;
      d->msg = messageHandler;
      d->ownedMessageHandler = ownedMessageHandler;

      d->smilesIndex = 0;

      // ensure the database filename has a .moldb extension
      std::size_t pos = filename.find(".moldb");
      if (pos == std::string::npos || pos != filename.size() - 6)
        d->msg->report(Error, make_string("File ", filename, " does not have a .moldb extension."));
    }

    static bool isLoaded(const Database *db)
    {
      if (!db->d->isLoaded) {
        db->d->msg->report(Warning, make_string("Databse ", db->d->filename, " is not loaded."));
        return false;
      }

      return true;
    }

    static std::string algorithm2string(enum Database::SubStructureSearchAlgorithm algorithm)
    {
      switch (algorithm) {
        case Database::MolDBRegularGraphVF2:
          return "MolDBRegularGraphVF2";
        case Database::MolDBReducedGraphVF2:
          return "MolDBReducedGraphVF2";
        case Database::OpenBabelVF2:
          return "OpenBabelVF2";
        case Database::OpenBabelSMARTS:
          return "OpenBabelSMARTS";
        case Database::ChemKitVF2:
          return "ChemKitVF2";
        case Database::RDKitVF2:
          return "RDKitVF2";
      }

      return std::string();
    }

    static void writeHeader(Database *db, std::fstream &ofs, std::ios_base::streampos endOfFile)
    {
      Impl::DatabaseHeader header;
      header.magicNumber = 0xABCDEFAB;
      header.endOfFile = endOfFile;
      header.graphCount = db->d->graphCount;
      header.vertexCount = db->d->vertexCount;
      header.edgeCount = db->d->edgeCount;
      header.graphVerticesCount = db->d->graphVerticesCount;
      header.graphEdgesCount = db->d->graphEdgesCount;
      header.vertexEdgesCount = db->d->vertexEdgesCount;
      header.edgePathCount = db->d->edgePathCount;

      ofs.write(reinterpret_cast<char*>(&header), sizeof(Impl::DatabaseHeader));
    }

    static bool readHeader(Database *db, std::ifstream &ifs)
    {
      Impl::DatabaseHeader header;
      ifs.read(reinterpret_cast<char*>(&header), sizeof(Impl::DatabaseHeader));
      if (header.magicNumber != 0xABCDEFAB) {
        db->d->msg->report(Error, make_string("File ", db->d->filename, " is not a MolDB database."));
        return false;
      }

      db->d->endOfFile = header.endOfFile;
      db->d->graphCount = header.graphCount;
      db->d->vertexCount = header.vertexCount;
      db->d->edgeCount = header.edgeCount;
      db->d->graphVerticesCount = header.graphVerticesCount;
      db->d->graphEdgesCount = header.graphEdgesCount;
      db->d->vertexEdgesCount = header.vertexEdgesCount;
      db->d->edgePathCount = header.edgePathCount;

      db->d->msg->report(Information, make_string(db->d->graphCount, " molecules in the databse."));
      return true;
    }

    static void writeFingerprintNames(Database *db, std::fstream &ofs)
    {
      // write the number of fingerprints
      unsigned int numFingerprints = db->d->fingerprintNames.size();
      ofs.write(reinterpret_cast<char*>(&numFingerprints), sizeof(unsigned int));
      // write the fingerprint names
      for (unsigned int i = 0; i < numFingerprints; ++i) {
        unsigned int fingerprintNameSize = db->d->fingerprintNames[i].size();
        ofs.write(reinterpret_cast<char*>(&fingerprintNameSize), sizeof(unsigned int));
        ofs.write(db->d->fingerprintNames[i].c_str(), fingerprintNameSize);
      }
    }

    static void readFingerprintNames(Database *db, std::ifstream &ifs)
    {
      // read the number of fingerprints
      unsigned int numFingerprints;
      ifs.read(reinterpret_cast<char*>(&numFingerprints), sizeof(unsigned int));
      // read the fingerprint names
      for (unsigned int i = 0; i < numFingerprints; ++i) {
        unsigned int fingerprintNameSize;
        ifs.read(reinterpret_cast<char*>(&fingerprintNameSize), sizeof(unsigned int));
        char *buffer = new char[fingerprintNameSize];
        ifs.read(buffer, fingerprintNameSize);
        db->d->fingerprintNames.push_back(std::string(buffer, fingerprintNameSize));
        delete [] buffer;
      }
    }

    static void writeGraph(Database *db, OpenBabel::OBMol *mol, std::fstream &ofs)
    {
      // create the graph
      GraphType *graph = OB::reducedGraph<GraphType>(mol);

      // update the count values
      db->d->graphCount++;
      db->d->vertexCount += graph->numVertices();
      db->d->edgeCount += graph->numEdges();
      db->d->graphVerticesCount += graph->numVertices() + 1;
      db->d->graphEdgesCount += graph->numEdges() + 1;
      for (std::size_t i = 0; i < graph->numVertices(); ++i)
        db->d->vertexEdgesCount += graph->vertex(i)->numEdges() + 1;
      for (std::size_t i = 0; i < graph->numEdges(); ++i)
        db->d->edgePathCount += graph->edge(i)->pathSize() + 1;

      //
      // write the vertices
      //
      // write the number of vertices
      unsigned int numVertices = graph->numVertices();
      ofs.write(reinterpret_cast<char*>(&numVertices), sizeof(unsigned int));
      // write the vertices
      for (std::size_t i = 0; i < numVertices; ++i) {
        // write the vertex semantics
        GraphType::VertexEdgeSemanticsType semantics = graph->vertex(i)->semantics();
        ofs.write(reinterpret_cast<char*>(&semantics), sizeof(GraphType::VertexEdgeSemanticsType));
      }
 
      //
      // write the edges
      //
      // write the number of edges
      unsigned int numEdges = graph->numEdges();
      ofs.write(reinterpret_cast<char*>(&numEdges), sizeof(unsigned int));
      // write the edges
      for (std::size_t i = 0; i < graph->numEdges(); ++i) {
        // write the source vertex index
        unsigned int source = graph->edge(i)->source()->index() + 1;
        ofs.write(reinterpret_cast<char*>(&source), sizeof(unsigned int));
        // write the target vertex index
        unsigned int target = graph->edge(i)->target() ? graph->edge(i)->target()->index() + 1 : 0;
        ofs.write(reinterpret_cast<char*>(&target), sizeof(unsigned int));
        // write the edge path size
        unsigned int pathSize = graph->edge(i)->pathSize();
        ofs.write(reinterpret_cast<char*>(&pathSize), sizeof(unsigned int));
        // write the edge path
        for (std::size_t j = 0; j < graph->edge(i)->pathSize(); ++j) {
          GraphType::VertexEdgeSemanticsType semantics = graph->edge(i)->path()[j];
          ofs.write(reinterpret_cast<char*>(&semantics), sizeof(GraphType::VertexEdgeSemanticsType));
        }
      }

      delete graph;
    }

    static void skipGraph(std::ifstream &ifs)
    {
      // read the number of vertices
      unsigned int numVertices;
      ifs.read(reinterpret_cast<char*>(&numVertices), sizeof(unsigned int));
      // skip the vertices
      ifs.seekg(numVertices * sizeof(GraphType::VertexEdgeSemanticsType), std::ios_base::cur);
      // read the number of edges
      unsigned int numEdges;
      ifs.read(reinterpret_cast<char*>(&numEdges), sizeof(unsigned int));
      // skip the edges
      for (std::size_t i = 0; i < numEdges; ++i) {
        // skip the source and target vertex index
        ifs.seekg(2 * sizeof(unsigned int), std::ios_base::cur);
        // read the edge path size
        unsigned int pathSize;
        ifs.read(reinterpret_cast<char*>(&pathSize), sizeof(unsigned int));
        // skip the edge path
        ifs.seekg(pathSize * sizeof(GraphType::VertexEdgeSemanticsType), std::ios_base::cur);
      }
    }
    
    static std::string memorySize(unsigned long bytes)
    {
      if (bytes > 1000000000)
        return make_string(static_cast<double>(bytes) / 1000000000, " GB");
      if (bytes > 1000000)
        return make_string(static_cast<double>(bytes) / 1000000, " MB");
      if (bytes > 1000)
        return make_string(static_cast<double>(bytes) / 1000, " KB");
      return make_string(bytes, " bytes");
    }

    static std::vector<std::size_t> fingerprintScreen(const Database *db, OpenBabel::OBMol *query, 
        const std::vector<std::string> &fingerprints, std::size_t maxResults, std::size_t startIndex,
        std::size_t stopIndex)
    {
      if (stopIndex > db->d->graphCount)
        stopIndex = db->d->graphCount;

      std::vector<std::size_t> fingerprintIndices;
      std::vector<Fingerprint> queryFingerprints;
      for (std::size_t i = 0; i < fingerprints.size(); ++i) {
        std::vector<std::string>::const_iterator pos = std::find(db->d->fingerprintNames.begin(), db->d->fingerprintNames.end(), fingerprints[i]);
        if (pos == db->d->fingerprintNames.end()) {
          return std::vector<std::size_t>();
        } else {
          fingerprintIndices.push_back(pos - db->d->fingerprintNames.begin());
          db->d->mutex.lock();
          queryFingerprints.push_back(Impl::fingerprint(fingerprints[i], query));
          db->d->mutex.unlock();
        }
      }

      if (!maxResults)
        maxResults = db->d->graphCount;
       
      std::vector<std::vector<std::size_t> > hits(fingerprintIndices.size());
      for (std::size_t j = 0; j < fingerprintIndices.size(); ++j) {
        std::vector<std::size_t> fphits = db->d->fingerprintIndexes[fingerprintIndices[j]]->contains(queryFingerprints[j]);
        std::sort(fphits.begin(), fphits.end());
        std::copy(fphits.begin(), fphits.end(), std::back_inserter(hits[j]));
        if (hits.size() >= maxResults && hits.size() == 1)
          return hits[0];
        maxResults -= fphits.size();
      }

      if (hits.size() == 1)
        return hits[0];

      // join the hits
      // select shortest hit list
      int shortest = 0;
      for (std::size_t i = 1; i < hits.size(); ++i)
        if (hits[i].size() < hits[shortest].size())
          shortest = i;
      

      std::vector<std::size_t> joined;
      for (std::size_t i = 0; i < hits[shortest].size(); ++i) {
        bool inAllHits = true;

        for (std::size_t j = 0; j < hits.size(); ++j) {
          if (i == shortest)
            continue;
          if (!std::binary_search(hits[j].begin(), hits[j].end(), hits[shortest][i])) {
            inAllHits = false;
            break;
          }
        }

        if (inAllHits)
          joined.push_back(hits[shortest][i]);
        if (joined.size() >= maxResults)
          return joined;
      }

      return joined;
    }
 
    static std::vector<std::pair<std::size_t, double> > similaritySearch(const Database *db, OpenBabel::OBMol *query, 
        double tanimotoThreshold, const std::vector<std::string> &fingerprints, std::size_t maxResults, std::size_t startIndex,
        std::size_t stopIndex)
    {
      if (stopIndex > db->d->graphCount)
        stopIndex = db->d->graphCount;

      std::vector<std::pair<std::size_t, double> > hits;
      std::vector<std::size_t> fingerprintIndices;
      std::vector<Fingerprint> queryFingerprints;
      for (std::size_t i = 0; i < fingerprints.size(); ++i) {
        std::vector<std::string>::const_iterator pos = std::find(db->d->fingerprintNames.begin(), db->d->fingerprintNames.end(), fingerprints[i]);
        if (pos == db->d->fingerprintNames.end()) {
          return hits;
        } else {
          fingerprintIndices.push_back(pos - db->d->fingerprintNames.begin());
          db->d->mutex.lock();
          queryFingerprints.push_back(Impl::fingerprint(fingerprints[i], query));
          db->d->mutex.unlock();
        }
      }

      if (!maxResults)
        maxResults = db->d->graphCount;
       
      for (std::size_t j = 0; j < fingerprintIndices.size(); ++j) {
        std::vector<std::pair<std::size_t, double> > fphits = db->d->fingerprintIndexes[fingerprintIndices[j]]->tanimoto(queryFingerprints[j], tanimotoThreshold);
        std::copy(fphits.begin(), fphits.end(), std::back_inserter(hits));
        if (hits.size() >= maxResults)
          return hits;
        maxResults -= fphits.size();
      }

      return hits;
    }


    class SubStructureSearchThread
    {
      public:
        SubStructureSearchThread(const Database *db, OpenBabel::OBMol *query, 
            enum Database::SubStructureSearchAlgorithm algorithm, std::size_t maxResults,
            std::size_t startIndex, std::size_t stopIndex, std::vector<std::size_t> &hits)
            : m_db(db), m_query(query), m_algorithm(algorithm), m_maxResults(maxResults),
            m_startIndex(startIndex), m_stopIndex(stopIndex), m_hits(hits)
        {
        }

        void operator()()
        {
          m_hits = m_db->fullSubStructureSearch(m_query, m_algorithm, m_maxResults, m_startIndex, m_stopIndex);
        }

      private:
        const Database *m_db;
        OpenBabel::OBMol *m_query;
        enum Database::SubStructureSearchAlgorithm m_algorithm;
        std::size_t m_maxResults;
        std::size_t m_startIndex;
        std::size_t m_stopIndex;
        std::vector<std::size_t> &m_hits;
    };

  };



  Database::Database(const std::string &filename) : d(new DatabasePrivate)
  {
    DatabaseImpl::construct(d, filename, 0, new StdOutMessageHandler, true);
  }


  Database::Database(const std::string &filename, CustomDataFunctor *customDataFunctor)
      : d(new DatabasePrivate)
  {
    DatabaseImpl::construct(d, filename, customDataFunctor, new StdOutMessageHandler, true);
  }

  Database::Database(const std::string &filename, MessageHandler *messageHandler)
      : d(new DatabasePrivate)
  {
    DatabaseImpl::construct(d, filename, 0, messageHandler, false);
  }

  Database::Database(const std::string &filename, CustomDataFunctor *customDataFunctor, 
      MessageHandler *messageHandler) : d(new DatabasePrivate)
  {
    DatabaseImpl::construct(d, filename, customDataFunctor, messageHandler, false);
  }

  Database::~Database()
  {
    if (d->msg && d->ownedMessageHandler)
      delete d->msg;
    delete d;
  }

  std::vector<std::string> Database::availableFingerprints()
  {
    OpenBabel::OBConversion();
    std::vector<std::string> fingerprints;

    // MolDB fingerprints
    fingerprints.push_back("MolDB::MolDB");
    fingerprints.push_back("MolDB::Andrew");
    fingerprints.push_back("MolDB::SubStruct");
    
    // OpenBabel fingerprints
    std::vector<std::string> obfingerprints;
    OpenBabel::OBFingerprint::ListAsVector("fingerprints", 0, obfingerprints);
    for (std::size_t i = 0; i < obfingerprints.size(); ++i) {
      std::string fp = "OpenBabel::" + obfingerprints[i].substr(0, obfingerprints[i].find(' '));
      fingerprints.push_back(fp);
    }

#ifdef HAVE_CHEMKIT
    // ChemKit fingerprints
    std::vector<std::string> chemkitfingerprints = chemkit::Fingerprint::fingerprints();
    for (std::size_t i = 0; i < chemkitfingerprints.size(); ++i) {
      std::string fp = "ChemKit::" + chemkitfingerprints[i].substr(0, chemkitfingerprints[i].find(' '));
      fingerprints.push_back(fp);
    }
#endif

#ifdef HAVE_RDKIT
    // RDKit fingerprints
    fingerprints.push_back("RDKit::Topological");
    fingerprints.push_back("RDKit::Layered");
    fingerprints.push_back("RDKit::AtomPairs");
    fingerprints.push_back("RDKit::TopologicalTorsions");
#endif

    return fingerprints;
  }

  void Database::create(const std::vector<std::string> &fingerprints)
  {
    std::ifstream ifs(d->filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (ifs) {
      unsigned int magicNumber;
      ifs.read(reinterpret_cast<char*>(&magicNumber), sizeof(unsigned int));
      if (magicNumber == 0xABCDEFAB) {
        d->msg->report(Error, "Database already exists, delete the database file first.");
        return;
      }
    }
    ifs.close();

    d->fingerprintNames = fingerprints;

    std::fstream ofs(d->filename.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    if (!ofs) {
      d->msg->report(Error, make_string("Could not open database ", d->filename, "."));
      return;
    }
    
    d->msg->report(Information, make_string("Creating database ", d->filename, "..."));

    // write the fingerprint names
    ofs.seekp(sizeof(Impl::DatabaseHeader));
    DatabaseImpl::writeFingerprintNames(this, ofs);
    for (std::size_t i = 0; i < 10000; ++i) {
      char data = 0x00;
      ofs.write(&data, 1);
    }
    std::ios_base::streampos endOfFile = ofs.tellp();
    // write the header
    ofs.seekp(0);
    DatabaseImpl::writeHeader(this, ofs, endOfFile);
  
    // create the smiles index
    d->msg->report(Information, make_string("Creating smiles index..."));
    d->smilesIndex = new SmilesIndexType(DatabaseImpl::baseFilename(d->filename));
    d->smilesIndex->create(d->msg);


    d->msg->report(Information, make_string("Creating fingerprint index..."));
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
      d->fingerprintStorages.push_back(new FingerprintStorageType(DatabaseImpl::baseFilename(d->filename), d->fingerprintNames[i], Impl::fingerprintSize(d->fingerprintNames[i])));
      d->fingerprintStorages.back()->create(d->msg);
      d->fingerprintIndexes.push_back(new FingerprintIndexType(DatabaseImpl::baseFilename(d->filename), d->fingerprintNames[i], Impl::fingerprintSize(d->fingerprintNames[i]), *d->fingerprintStorages[i]));
      d->fingerprintIndexes.back()->create(d->msg);
    }

    d->msg->report(Information, make_string("Database ", d->filename, " successfuly created."));
 
    d->isCreated = true;
  }
      
  void Database::addFingerprint(const std::string &fingerprint)
  { 
    FingerprintStorageType *storage;
    FingerprintIndex *index;

    std::vector<std::string>::iterator nameIter = std::find(d->fingerprintNames.begin(), d->fingerprintNames.end(), fingerprint);
    if (nameIter == d->fingerprintNames.end()) {
      d->msg->report(Information, make_string("Adding fingerprint ", fingerprint, "."));

      d->fingerprintNames.push_back(fingerprint);

      d->fingerprintStorages.push_back(new FingerprintStorageType(DatabaseImpl::baseFilename(d->filename), fingerprint, Impl::fingerprintSize(fingerprint)));
      d->fingerprintIndexes.push_back(new FingerprintIndexType(DatabaseImpl::baseFilename(d->filename), fingerprint, Impl::fingerprintSize(fingerprint), *d->fingerprintStorages.back()));

      storage = d->fingerprintStorages.back();
      index = d->fingerprintIndexes.back();
    } else {
      d->msg->report(Information, make_string("Updating fingerprint ", fingerprint, "."));

      int nameIndex = nameIter - d->fingerprintNames.begin();
      delete d->fingerprintStorages[nameIndex];
      delete d->fingerprintIndexes[nameIndex];
      d->fingerprintStorages[nameIndex] = new FingerprintStorageType(DatabaseImpl::baseFilename(d->filename), fingerprint, Impl::fingerprintSize(fingerprint));
      d->fingerprintIndexes[nameIndex] = new FingerprintIndexType(DatabaseImpl::baseFilename(d->filename), fingerprint, Impl::fingerprintSize(fingerprint), *d->fingerprintStorages.back());

      storage = d->fingerprintStorages[nameIndex];
      index = d->fingerprintIndexes[nameIndex];
      storage->clear();
      index->clear();
    }
 
    storage->create(d->msg);
    index->create(d->msg);
    storage->startInsert(d->msg);
    index->startInsert(d->msg);
   
    for (std::size_t i = 0; i < d->smilesIndex->size(); ++i) {
      boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(d->smilesIndex->get(i)));
      Fingerprint fp = Impl::fingerprint(fingerprint, mol.get());
      storage->insert(d->msg, fp);
      index->insert(d->msg, fp);
      if (i % 10000 == 0)
        d->msg->report(Information, make_string(i, "/", d->smilesIndex->size()));
    }

    storage->commit();
    index->commit();
    storage->stopInsert(d->msg);
    index->stopInsert(d->msg);


    std::fstream ofs(d->filename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    if (!ofs) {
      d->msg->report(Error, make_string("Could not open database ", d->filename, "."));
      return;
    }
    
    // write the fingerprint names
    ofs.seekp(sizeof(Impl::DatabaseHeader));
    DatabaseImpl::writeFingerprintNames(this, ofs);
  }
  
  void Database::create(const std::string &fingerprint)
  {
    std::vector<std::string> fingerprints(1 , fingerprint);
    create(fingerprints);
  }
 
  bool Database::load(const std::vector<std::string> &fingerprints)
  {
    Timer timer;

    // clear any existing records
    d->graphCount = 0;
    d->vertexCount = 0;
    d->edgeCount = 0;
    d->graphVerticesCount = 0;
    d->graphEdgesCount = 0;
    d->vertexEdgesCount = 0;
    d->edgePathCount = 0;

    d->fingerprintNames.clear();
    d->positions.clear();
    d->graphs.clear();

    TinyGraphMemoryPool::instance().deallocate();
    TinyVertexMemoryPool::instance().deallocate();
    TinyEdgeMemoryPool::instance().deallocate();
    TinyGraphVerticesMemoryPool::instance().deallocate();
    TinyGraphEdgesMemoryPool::instance().deallocate();
    TinyVertexEdgesMemoryPool::instance().deallocate();
    TinyEdgePathMemoryPool::instance().deallocate();

    // open the file
    std::ifstream ifs(d->filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!ifs) {
      d->msg->report(Error, make_string("Could not open database ", d->filename, "."));
      return false;
    }

    // read the header
    if (!DatabaseImpl::readHeader(this, ifs))
      return false;
    // read the fingerprint names
    DatabaseImpl::readFingerprintNames(this, ifs);
    // skip insert resume & free space
    ifs.seekg(10000, std::ios_base::cur);

    if (fingerprints.size()) {
      std::vector<std::string> dbFingerprints = d->fingerprintNames;
      d->fingerprintNames.clear();
      for (std::size_t i = 0; i < dbFingerprints.size(); ++i)
        if (std::find(fingerprints.begin(), fingerprints.end(), dbFingerprints[i]) != fingerprints.end())
          d->fingerprintNames.push_back(dbFingerprints[i]);
    }

    // create the storages
    if (d->fingerprintStorages.size() != d->fingerprintNames.size())
      for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i)
        d->fingerprintStorages.push_back(new FingerprintStorageType(DatabaseImpl::baseFilename(d->filename), d->fingerprintNames[i], Impl::fingerprintSize(d->fingerprintNames[i])));

    // create the indexes
    if (d->fingerprintIndexes.size() != d->fingerprintNames.size())
      for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i)
        d->fingerprintIndexes.push_back(new FingerprintIndexType(DatabaseImpl::baseFilename(d->filename), d->fingerprintNames[i], Impl::fingerprintSize(d->fingerprintNames[i]), *d->fingerprintStorages[i]));
    if (!d->smilesIndex)
      d->smilesIndex = new SmilesIndexType(DatabaseImpl::baseFilename(d->filename));

    // reserve memory to avoid unneeded copying
    d->positions.reserve(d->graphCount + 50000);
    d->smilesIndex->reserve(d->graphCount + 50000, 60); // assume average smiles length of 60 bytes
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
      d->fingerprintStorages[i]->reserve(d->graphCount + 50000);
      d->fingerprintIndexes[i]->reserve(d->graphCount + 50000);
    }
    d->graphs.reserve(d->graphCount + 500000);
    // reserve memory for the various objects
    unsigned long factor = (d->maxGraphsInMemory < d->graphCount) ? d->graphCount / d->maxGraphsInMemory : 1;
    TinyGraphMemoryPool::instance(d->graphCount / factor + 500000);
    TinyVertexMemoryPool::instance(d->graphVerticesCount / factor + 500000);
    TinyEdgeMemoryPool::instance(d->graphEdgesCount / factor + 500000);
    TinyGraphVerticesMemoryPool::instance(d->vertexCount / factor + 500000);
    TinyGraphEdgesMemoryPool::instance(d->vertexEdgesCount / factor + 500000);
    TinyVertexEdgesMemoryPool::instance(d->edgeCount / factor + 500000);
    TinyEdgePathMemoryPool::instance(d->edgePathCount / factor + 500000);

    std::size_t twoPercent = d->graphCount / 50;

    // start reading the records
    for (std::size_t i = 0; i < d->graphCount; ++i) {
      // create a new record
      //d->graphs.push_back(0);

//      if (d->graphs.size() % 100 == 0)
//        d->msg->report(Information, make_string("Loading molecule # ", i + 1), "Database::load()");

      // record the graph position
      //d->positions.push_back(ifs.tellg());
      // skip the graph
      //DatabaseImpl::skipGraph(ifs);

      if (i % twoPercent == 0) {
        d->msg->report(Information, progressBar(i, d->graphCount));
      }
    }

    d->msg->report(Information, "Loading smiles index...");
    d->smilesIndex->load(d->msg);
    d->msg->report(Information, "Loading fingerprint indexes...");
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
      d->fingerprintStorages[i]->load(d->msg);
      d->fingerprintIndexes[i]->load(d->msg);
    }
 
    d->msg->report(Information, make_string("Loaded ", d->graphCount, " molecules in ", timer.elapsed(), " seconds."));
    d->isCreated = true;
    d->isLoaded = true;
    return true;
  }

  void Database::insert(const std::string &filename, bool unique, int resumeInterval)
  {
    if (!DatabaseImpl::isLoaded(this))
      return;

    Timer timer;
    unsigned long startIndex = 0, startPos;

    // open the database file for reading
    std::ifstream ifs(d->filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!ifs) {
      d->msg->report(Error, make_string("Could not open database file ", d->filename, "."));
      return;
    }

    // get the end of file position
    if (!DatabaseImpl::readHeader(this, ifs)) {
      d->msg->report(Error, make_string("Database ", d->filename, " does not exist, create it first."));
      return;
    }
        
    d->msg->report(Information, "Inserting molecules...");

    d->fingerprintNames.clear();
    DatabaseImpl::readFingerprintNames(this, ifs);
    std::ios_base::streampos endOfFile = d->endOfFile;
    std::ios_base::streampos position = ifs.tellg();
    // check if we need to resume
    unsigned int filenameLength;
    ifs.read(reinterpret_cast<char*>(&filenameLength), sizeof(unsigned int));
    if (filenameLength) {
      char *buffer = new char[filenameLength];
      ifs.read(buffer, filenameLength);
      std::string resumeFilename(buffer, filenameLength);
      //std::cout << resumeFilename << std::endl;
      delete [] buffer;
      if (resumeFilename == d->filename) {
        // read the start index to resume from
        ifs.read(reinterpret_cast<char*>(&startIndex), sizeof(unsigned long));
        // read the insert file stream position to resume from
        ifs.read(reinterpret_cast<char*>(&startPos), sizeof(unsigned long));
        d->msg->report(Information, make_string("Resuming inserting from index ", startIndex, "."));
      }
    }
    ifs.close();

    // open the database file for writing
    std::fstream ofs(d->filename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    if (!ofs) {
      d->msg->report(Error, make_string("Could not open database ", d->filename, "."));
      return;
    }
    ofs.seekp(endOfFile);


    OB::File file(filename.c_str());
    if (startIndex)
      file.seekg(startPos);

    d->smilesIndex->startInsert(d->msg);
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
      d->fingerprintStorages[i]->startInsert(d->msg);
      d->fingerprintIndexes[i]->startInsert(d->msg);
    }
 
    unsigned long graphCount = startIndex;
    OpenBabel::OBMol *mol;
    while ((mol = file.next())) {
      graphCount++;

      d->msg->report(Information, make_string("Inserting molecule # ", graphCount));
      d->msg->report(Information, make_string("    title: ", mol->GetTitle()));

      if (unique) {
        std::string smiles = OB::File::writeSmiles(mol);
        if (d->smilesIndex->contains(smiles) != maxIndex()) {
          d->msg->report(Information, "    not unique");
          continue;
        }
      }

      
      if (graphCount % (resumeInterval - 1) == 0)
        startPos = file.tellg();
      if (graphCount % resumeInterval == 0) {
        // write resume data
        endOfFile = ofs.tellp();
        ofs.seekp(position);
        // write filename size
        filenameLength = d->filename.size();
        ofs.write(reinterpret_cast<char*>(&filenameLength), sizeof(unsigned int));
        // write filename
        ofs.write(d->filename.c_str(), filenameLength);
        // write start index
        startIndex = graphCount - 1;
        ofs.write(reinterpret_cast<char*>(&startIndex), sizeof(unsigned long));
        // write the start position in the insert file
        ofs.write(reinterpret_cast<char*>(&startPos), sizeof(unsigned long));
        // write the header
        ofs.seekp(0);
        DatabaseImpl::writeHeader(this, ofs, endOfFile);
        ofs.flush();
        ofs.seekp(endOfFile);
        d->smilesIndex->commit();
        for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
          d->fingerprintStorages[i]->commit();
          d->fingerprintIndexes[i]->commit();
        }
      }

      // create a new record
      d->graphs.push_back(0);

      // record the graph position in the file      
      d->positions.push_back(ofs.tellp());
      // write the graph for the molecule
      DatabaseImpl::writeGraph(this, mol, ofs);
      // insert the smiles in the smiles index
      d->smilesIndex->insert(d->msg, OB::File::writeSmiles(mol));
      // insert the fingerprints in the fingerprint indexes
      for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
        Fingerprint fp = Impl::fingerprint(d->fingerprintNames[i], mol);
        d->fingerprintStorages[i]->insert(d->msg, fp);
        d->fingerprintIndexes[i]->insert(d->msg, fp);
      }
 
      delete mol;
    }

    endOfFile = ofs.tellp();
    // write the header
    ofs.seekp(0);
    DatabaseImpl::writeHeader(this, ofs, endOfFile);

    d->smilesIndex->commit();
    d->smilesIndex->stopInsert(d->msg);
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
      d->fingerprintStorages[i]->commit();
      d->fingerprintStorages[i]->stopInsert(d->msg);
      d->fingerprintIndexes[i]->commit();
      d->fingerprintIndexes[i]->stopInsert(d->msg);
    }
 
    // clear the resume data
    ofs.seekp(position);
    filenameLength = 0;
    ofs.write(reinterpret_cast<char*>(&filenameLength), sizeof(unsigned int));
 
    d->msg->report(Information, make_string("Inserted ", graphCount, " molecules in ", timer.elapsed(), " seconds."));
  }
  
  std::size_t Database::numRecords() const
  {
    if (!DatabaseImpl::isLoaded(this))
      return 0;
    return d->graphCount;
  }
  
  Database::Record Database::record(std::size_t index) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return Database::Record();

    Database::Record record;
    record.position = d->positions[index];
    for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i)
      record.fingerprints.push_back((*d->fingerprintStorages[i])[index]);
 
    record.graph = d->graphs[index];
    record.canonicalSmiles = d->smilesIndex->get(index);
    return record;
  }

  std::vector<std::string> Database::databaseFingerprints() const
  {
    DatabaseImpl::isLoaded(this);
    return d->fingerprintNames;
  }
  
  void Database::databaseInformation() const
  {
    DatabaseImpl::isLoaded(this);

    // get the database file size
    std::ifstream ifs(d->filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!ifs) {
      d->msg->report(Error, make_string("Could not open databse ", d->filename, "."));
      return;
    }
    ifs.seekg(0, std::ios_base::end);
    unsigned long databaseFileSize = ifs.tellg();
    ifs.close();

    d->msg->report(Information, "------------------------------------------------------------");
    // filename
    d->msg->report(Information, make_string("Database: ", d->filename));
    // database size
    d->msg->report(Information, make_string("File Size: ", DatabaseImpl::memorySize(databaseFileSize)));
    // number of molecules
    d->msg->report(Information, make_string("Number of Molecules: ", d->graphCount));
    // memory usage
    d->msg->report(Information, "Memory Usage:");
    unsigned long memory = sizeof(Database) + sizeof(DatabasePrivate);
    if (d->graphCount) {
      // positions
      memory += d->graphCount * sizeof(unsigned long);
      d->msg->report(Information, make_string("    Positions: ", DatabaseImpl::memorySize(d->graphCount * sizeof(unsigned long))));
      // fingerprints
      d->msg->report(Information, "    Fingerprints: ");
      d->msg->report(Information, "        Storage: ");
      unsigned long fingerprintSize = 0;
      for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
        fingerprintSize += d->fingerprintStorages[i]->memoryUsage();
        d->msg->report(Information, make_string("            ", d->fingerprintNames[i], " (", d->fingerprintStorages[i]->numBits(), " bits): ", 
              DatabaseImpl::memorySize(d->fingerprintStorages[i]->memoryUsage())));
      }
      d->msg->report(Information, "        Indexes: ");
      for (std::size_t i = 0; i < d->fingerprintNames.size(); ++i) {
        d->msg->report(Information, make_string("            ", d->fingerprintIndexes[i]->name()));
        fingerprintSize += d->fingerprintIndexes[i]->memoryUsage();
        d->msg->report(Information, make_string("            ", d->fingerprintNames[i], " (", d->fingerprintStorages[i]->numBits(), " bits): ", 
              DatabaseImpl::memorySize(d->fingerprintIndexes[i]->memoryUsage())));
      }
      memory += fingerprintSize;
      d->msg->report(Information, "");
      d->msg->report(Information, make_string("        Total: ", DatabaseImpl::memorySize(fingerprintSize)));
      // graph pointers
      memory += d->graphCount * sizeof(TinyGraphType*);
      d->msg->report(Information, make_string("    Pointers to Graphs: ", DatabaseImpl::memorySize(d->graphCount * sizeof(TinyGraphType*))));
      // canonical smiles
      memory += d->smilesIndex->memoryUsage();
      d->msg->report(Information, make_string("    Canonical Smiles Index: ", DatabaseImpl::memorySize(d->smilesIndex->memoryUsage())));
      d->msg->report(Information, "");
      d->msg->report(Information, make_string("    Total: ", DatabaseImpl::memorySize(memory)));
    }

    d->msg->report(Information, "------------------------------------------------------------");
  }
 
  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Exact Structure Searching
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  std::size_t Database::exactStructureSearch(OpenBabel::OBMol *query) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return false;

    Timer timer;
    d->msg->report(Information, "Exact structure search");

    std::string smiles = OB::File::writeSmiles(query);
    std::size_t index = d->smilesIndex->contains(smiles);

    d->msg->report(Information, make_string("Searched ", d->smilesIndex->size(), " molecules in ", timer.elapsed(), " seconds."));
    if (index != maxIndex())
      d->msg->report(Information, make_string("    Found: ", index));
    else
      d->msg->report(Information, "    Not Found");

    return index;
  }

  std::size_t Database::exactStructureSearch(const std::string &smiles) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return exactStructureSearch(mol.get());
  }
 
  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Fingerprint Screening
  //
  //
  //////////////////////////////////////////////////////////////////////////////
  
  std::vector<std::size_t> Database::fingerprintScreen(OpenBabel::OBMol *query,
      const std::vector<std::string> &fingerprints, std::size_t maxResults) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return std::vector<std::size_t>();

    Timer timer;
    std::stringstream ss;
    for (std::size_t i = 0; i < fingerprints.size(); ++i) {
      ss << fingerprints[i];
      if (i + 1 < fingerprints.size())
        ss << " ";
    }
    d->msg->report(Information, make_string("Fingerprint screening using fingerprint(s) ", ss.str(), "."));

    for (std::size_t i = 0; i < fingerprints.size(); ++i) {
      if (std::find(d->fingerprintNames.begin(), d->fingerprintNames.end(), fingerprints[i]) == d->fingerprintNames.end()) {
        d->msg->report(Warning, make_string("Fingerprint ", fingerprints[i], " not in databse."));
        return std::vector<std::size_t>();
      } 
    }

    std::vector<std::size_t> hits = DatabaseImpl::fingerprintScreen(this, query, fingerprints, maxResults, 0, maxIndex());

    d->msg->report(Information, make_string("Screened ", d->graphCount * fingerprints.size(), " fingerprints in ", timer.elapsed(), " seconds."));
    d->msg->report(Information, make_string("    Hits: ", hits.size()));

    return hits;   
  }
  
  std::vector<std::size_t> Database::fingerprintScreen(OpenBabel::OBMol *query,
      const std::string &fingerprint, std::size_t maxResults) const
  {
    std::vector<std::string> fingerprints(1, fingerprint);
    return fingerprintScreen(query, fingerprints, maxResults);
  }
      
  std::vector<std::size_t> Database::fingerprintScreen(const std::string &smiles,
      const std::vector<std::string> &fingerprints, std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return fingerprintScreen(mol.get(), fingerprints, maxResults);
  }
      
  std::vector<std::size_t> Database::fingerprintScreen(const std::string &smiles,
      const std::string &fingerprint, std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return fingerprintScreen(mol.get(), fingerprint, maxResults);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Similarity Searching
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<std::size_t, double> > Database::similaritySearch(OpenBabel::OBMol *query,
      double tanimotoThreshold, const std::vector<std::string> &fingerprints,
      std::size_t maxResults) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return std::vector<std::pair<std::size_t, double> >();

    Timer timer;
    std::stringstream ss;
    for (std::size_t i = 0; i < fingerprints.size(); ++i) {
      ss << fingerprints[i];
      if (i + 1 < fingerprints.size())
        ss << " ";
    }
    d->msg->report(Information, make_string("Similarity search (>", tanimotoThreshold, ") using fingerprint(s) ", ss.str(), "."));

    for (std::size_t i = 0; i < fingerprints.size(); ++i) {
      if (std::find(d->fingerprintNames.begin(), d->fingerprintNames.end(), fingerprints[i]) == d->fingerprintNames.end()) {
        d->msg->report(Warning, make_string("Fingerprint ", fingerprints[i], " not in databse."));
        return std::vector<std::pair<std::size_t, double> >();
      } 
    }

    std::vector<std::pair<std::size_t,double> > hits = DatabaseImpl::similaritySearch(this, query, tanimotoThreshold, fingerprints, maxResults, 0, maxIndex());

    d->msg->report(Information, make_string("Searched ", d->graphCount * fingerprints.size(), " fingerprints in ", timer.elapsed(), " seconds."));
    d->msg->report(Information, make_string("    Hits: ", hits.size()));

    return hits;
  }

  std::vector<std::pair<std::size_t, double> > Database::similaritySearch(OpenBabel::OBMol *query,
      double tanimotoThreshold, const std::string &fingerprint,
      std::size_t maxResults) const
  {
    std::vector<std::string> fingerprints(1, fingerprint);
    return similaritySearch(query, tanimotoThreshold, fingerprints, maxResults);
  }

  std::vector<std::pair<std::size_t, double> > Database::similaritySearch(const std::string &smiles,
      double tanimotoThreshold, const std::vector<std::string> &fingerprints,
      std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return similaritySearch(mol.get(), tanimotoThreshold, fingerprints, maxResults);
  }

  std::vector<std::pair<std::size_t, double> > Database::similaritySearch(const std::string &smiles,
      double tanimotoThreshold, const std::string &fingerprint,
      std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return similaritySearch(mol.get(), tanimotoThreshold, fingerprint, maxResults);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Full SubStructure Searching
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  std::vector<std::size_t> Database::fullSubStructureSearch(OpenBabel::OBMol *query,
      enum Database::SubStructureSearchAlgorithm algorithm, std::size_t maxResults,
      std::size_t startIndex, std::size_t stopIndex) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return std::vector<std::size_t>();

    std::vector<std::size_t> hits;

    std::string algorithmString = DatabaseImpl::algorithm2string(algorithm);
    if (algorithmString.empty()) {
      d->msg->report(Warning, "SubStructure search algorithm not supported.");
      return hits;
    }

    d->msg->report(Information, make_string("SubStructure search using ", algorithmString, " algorithm."));

    Timer timer;
    unsigned long numSearched = 0;

    if (!maxResults)
      maxResults = std::numeric_limits<std::size_t>::max();

    if (algorithm == MolDBRegularGraphVF2) {
      // create the query graph
      d->mutex.lock();
      boost::shared_ptr<GraphType> queryGraph(OB::graph<GraphType>(query));
      d->mutex.unlock();

      if (stopIndex > d->graphCount)
        stopIndex = d->graphCount;

      for (std::size_t i = startIndex; i < stopIndex; ++i) {
        ++numSearched;
        // create the queried graph
        d->mutex.lock();
        boost::shared_ptr<OpenBabel::OBMol> queried(OB::File::readSmiles(d->smilesIndex->get(i)));
        d->mutex.unlock();
        boost::shared_ptr<GraphType> queriedGraph(OB::graph<GraphType>(queried.get()));
        // run the isomorphism algorithm
        GraphIsomorphismVF2<false, GraphType> iso;
        if (iso.isSubgraph(queryGraph.get(), queriedGraph.get()))
          hits.push_back(i);
        
        if (numSearched % 100 == 0)
          d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }

    if (algorithm == OpenBabelSMARTS) {
      OpenBabel::OBSmartsPattern sp;
      sp.Init(OB::File::writeSmiles(query));

      for (std::size_t i = 0; i < d->graphCount; ++i) {
        ++numSearched;
        // create the queried molecule
        boost::shared_ptr<OpenBabel::OBMol> queried(OB::File::readSmiles(d->smilesIndex->get(i)));
        // run the isomorphism algorithm
        if (sp.Match(*queried.get()))
          hits.push_back(i);
        
        if (numSearched % 100 == 0)
          d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }

#ifdef HAVE_CHEMKIT
    if (algorithm == ChemKitVF2) {
      boost::shared_ptr<chemkit::Molecule> chemkitQuery(Impl::chemKitMolecule(query));
      chemkit::SubstructureQuery q(chemkitQuery);

      for (std::size_t i = 0; i < d->graphCount; ++i) {
        ++numSearched;
        // create the queried molecule
        boost::shared_ptr<chemkit::Molecule> queried(Impl::chemKitMolecule(OB::File::readSmiles(d->smilesIndex->get(i))));
        // run the isomorphism algorithm
        if (q.matches(queried.get()))
          hits.push_back(i);
        
        if (numSearched % 100 == 0)
          d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }
#endif

  }
      
  std::vector<std::size_t> Database::fullSubStructureSearch(const std::string &smiles, 
      enum Database::SubStructureSearchAlgorithm algorithm, std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return fullSubStructureSearch(mol.get(), algorithm, maxResults);
  }
  
  std::vector<std::size_t> Database::fullSubStructureSearchThreaded(OpenBabel::OBMol *query,
      enum SubStructureSearchAlgorithm algorithm, int numThreads, std::size_t maxResults) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return std::vector<std::size_t>();

    boost::thread_group threads;
    std::size_t part = d->graphCount / numThreads;
    std::vector<std::vector<std::size_t> > threadHits(numThreads);


    for (int i = 0; i < numThreads; ++i) {
      std::size_t startIndex = i * part;
      std::size_t stopIndex = (i + 1 == numThreads) ? maxIndex() : (i + 1) * part - 1;
      threads.create_thread(DatabaseImpl::SubStructureSearchThread(this, query, algorithm, maxResults, startIndex, stopIndex, threadHits[i]));
    }

    threads.join_all();

    std::vector<std::size_t> hits;
    for (int i = 0; i < numThreads; ++i)
      std::copy(threadHits[i].begin(), threadHits[i].end(), std::back_inserter(hits));

    return hits;
  }
 
  std::vector<std::size_t> Database::fullSubStructureSearchThreaded(const std::string &smiles,
      enum SubStructureSearchAlgorithm algorithm, int numThreads, std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return fullSubStructureSearchThreaded(mol.get(), algorithm, numThreads, maxResults);
  }
 
  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Quick SubStructure Searching
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  std::vector<std::size_t> Database::quickSubStructureSearch(OpenBabel::OBMol *query,
      const std::vector<std::string> &fingerprints, enum Database::SubStructureSearchAlgorithm algorithm,
      std::size_t maxResults) const
  {
    if (!DatabaseImpl::isLoaded(this))
      return std::vector<std::size_t>();

    std::vector<std::size_t> hits;

    std::string algorithmString = DatabaseImpl::algorithm2string(algorithm);
    if (algorithmString.empty()) {
      d->msg->report(Warning, "SubStructure search algorithm not supported.");
      return hits;
    }

    d->msg->report(Information, make_string("SubStructure search using ", algorithmString, " algorithm."));

    Timer timer;
    unsigned long numSearched = 0;

    std::vector<std::size_t> screened = DatabaseImpl::fingerprintScreen(this, query, fingerprints, maxResults, 0, maxIndex());

    if (!maxResults)
      maxResults = std::numeric_limits<std::size_t>::max();


    /*
    if (algorithm == MolDBRegularGraphVF2) {
      // create the query graph
      d->mutex.lock();
      boost::shared_ptr<GraphType> queryGraph(OB::graph<GraphType>(query));
      d->mutex.unlock();

      if (stopIndex > d->graphCount)
        stopIndex = d->graphCount;

      for (std::size_t i = startIndex; i < stopIndex; ++i) {
        ++numSearched;
        // create the queried graph
        d->mutex.lock();
        boost::shared_ptr<OpenBabel::OBMol> queried(OB::File::readSmiles(d->smilesIndex->get(i)));
        d->mutex.unlock();
        boost::shared_ptr<GraphType> queriedGraph(OB::graph<GraphType>(queried.get()));
        // run the isomorphism algorithm
        GraphIsomorphismVF2<false, GraphType> iso;
        if (iso.isSubgraph(queryGraph.get(), queriedGraph.get()))
          hits.push_back(i);
        
        if (numSearched % 100 == 0)
          d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }
    */

    if (algorithm == OpenBabelSMARTS) {
      OpenBabel::OBSmartsPattern sp;
      sp.Init(OB::File::writeSmiles(query));

      for (std::size_t i = 0; i < screened.size(); ++i) {
        ++numSearched;
        // create the queried molecule
        boost::shared_ptr<OpenBabel::OBMol> queried(OB::File::readSmiles(d->smilesIndex->get(screened[i])));
        // run the isomorphism algorithm
        if (sp.Match(*queried.get()))
          hits.push_back(screened[i]);
        
        //if (numSearched % 1000 == 0)
        //  d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }

#ifdef HAVE_CHEMKIT
    if (algorithm == ChemKitVF2) {
      boost::shared_ptr<chemkit::Molecule> chemkitQuery(Impl::chemKitMolecule(query));
      chemkit::SubstructureQuery q(chemkitQuery);

      for (std::size_t i = 0; i < screened.size(); ++i) {
        ++numSearched;
        // create the queried molecule
        boost::shared_ptr<chemkit::Molecule> queried(Impl::chemKitMolecule(OB::File::readSmiles(d->smilesIndex->get(screened[i]))));
        // run the isomorphism algorithm
        if (q.matches(queried.get()))
          hits.push_back(screened[i]);
        
        if (numSearched % 1000 == 0)
          d->msg->report(Information, make_string("    ", hits.size(), " / ", numSearched, " hits"));
        
        // return the results if maxResults is reached
        if (hits.size() == maxResults)
          return hits;
      }
    
      d->msg->report(Information, make_string("Searched ", numSearched, " molecules in ", timer.elapsed(), " seconds."));
      d->msg->report(Information, make_string("    Hits: ", hits.size()));
        
      return hits;
    }
#endif


    return hits;
  }
      
  std::vector<std::size_t> Database::quickSubStructureSearch(OpenBabel::OBMol *query,
      const std::string &fingerprint, enum Database::SubStructureSearchAlgorithm algorithm,
      std::size_t maxResults) const
  {
    std::vector<std::string> fingerprints(1, fingerprint);
    return quickSubStructureSearch(query, fingerprints, algorithm, maxResults);
  }
  
  std::vector<std::size_t> Database::quickSubStructureSearch(const std::string &smiles,
      const std::vector<std::string> &fingerprints, enum Database::SubStructureSearchAlgorithm algorithm,
      std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return quickSubStructureSearch(mol.get(), fingerprints, algorithm, maxResults);
  }

  std::vector<std::size_t> Database::quickSubStructureSearch(const std::string &smiles,
      const std::string &fingerprint, enum Database::SubStructureSearchAlgorithm algorithm,
      std::size_t maxResults) const
  {
    boost::shared_ptr<OpenBabel::OBMol> mol(OB::File::readSmiles(smiles));
    return quickSubStructureSearch(mol.get(), fingerprint, algorithm, maxResults);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  //
  // Custom Data Functors
  //
  //
  //////////////////////////////////////////////////////////////////////////////

  void TitleCustomDataFunctor::write(std::fstream &ofs, OpenBabel::OBMol *molecule)
  {
    std::string title = molecule->GetTitle();
    unsigned int titleSize = title.size();
    ofs.write(reinterpret_cast<char*>(&titleSize), sizeof(unsigned int));
    ofs.write(title.c_str(), titleSize);
  }

  void TitleCustomDataFunctor::read(std::ifstream &ifs)
  {
    unsigned int titleSize;
    ifs.read(reinterpret_cast<char*>(&titleSize), sizeof(unsigned int));
    char *buffer = new char[titleSize];
    ifs.read(buffer, titleSize);
    m_titles.push_back(std::string(buffer, titleSize));
    delete [] buffer; 
  }
  
  std::string TitleCustomDataFunctor::title(std::size_t index) const
  {
    return std::string(m_titles[index]);
  }

  void MolecularWeightCustomDataFunctor::write(std::fstream &ofs, OpenBabel::OBMol *molecule)
  {
    double weight = molecule->GetMolWt();
    ofs.write(reinterpret_cast<char*>(&weight), sizeof(double));
  }
  
  void MolecularWeightCustomDataFunctor::read(std::ifstream &ifs)
  {
    double weight;
    ifs.read(reinterpret_cast<char*>(&weight), sizeof(double));
    m_weights.push_back(weight);
  }
  
  double MolecularWeightCustomDataFunctor::weight(std::size_t index) const
  {
    return m_weights[index];
  }

  CombineCustomDataFunctor::CombineCustomDataFunctor(const std::vector<Database::CustomDataFunctor*> &functors) : m_functors(functors)
  {
  }

  void CombineCustomDataFunctor::write(std::fstream &ofs, OpenBabel::OBMol *molecule)
  {
    for (std::size_t i = 0; i < m_functors.size(); ++i)
      m_functors[i]->write(ofs, molecule);
  }

  void CombineCustomDataFunctor::read(std::ifstream &ifs)
  {
    for (std::size_t i = 0; i < m_functors.size(); ++i)
      m_functors[i]->read(ifs);
  }
          

}
