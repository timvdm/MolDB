#ifndef MOLDB_DATABASE_H
#define MOLDB_DATABASE_H

#include "typetraits.h"
#include "graph.h"
#include "tinygraph.h"
#include "stringvector.h"
#include "fingerprint.h"

#include <openbabel/mol.h>

#include <boost/dynamic_bitset.hpp>

#include <string>

namespace MolDB {

  struct DatabasePrivate;
  class MessageHandler;

  /**
   * @class Database database.h
   *
   * The Database class is the main class in the MolDB library and most users
   * will only need this class. 
   *
   * @section Messages
   * Messages are handled by Database::MessageHandler subclasses. There are a
   * number of these subclasses already provided. These include the
   * StdOutMessageHandler, StdStreamMessageHandler and RecordingMessageHandler
   * classes. Adding new message handlers is easy and does not require to
   * recompile the library.
   *
   * @section Custom Data
   * By default, a MolDB database contains only data to perform searches. To
   * incorporate additional data Database::CustomDataFunctor subclasses can be
   * used. 
   */
  class Database
  {
      typedef TinyGraph<TypeTraits::EmptyType, unsigned char, unsigned short> TinyGraphType;
    public:
      enum SubStructureSearchAlgorithm
      {
        MolDBRegularGraphVF2,
        MolDBReducedGraphVF2,
        OpenBabelVF2,
        OpenBabelSMARTS,
        ChemKitVF2,
        RDKitVF2,
        InvalidSubStructureSearchAlgorithm
      };

      /**
       * A database record.
       */
      struct Record
      {
        /**
         * Constructor.
         */
        Record() : position(-1), graph(0)
        {
        }

        /**
         * The stream position where the molecule graph is located.
         */
        std::ios_base::streampos position;
        /**
         * The fingerprints for the molecule.
         */
        std::vector<Fingerprint> fingerprints;
        /**
         * The fingerprint bit counts.
         */
        std::vector<unsigned int> bitcounts;
        /**
         * The pointer to the graph. If this pointer is 0, the molecule
         * is not loaded in memory yet.
         */
        TinyGraphType* graph;
        /**
         * The canonical smiles.
         */
        std::string canonicalSmiles;
      };

      class CustomDataFunctor
      {
        public:
          virtual void write(std::fstream &ofs, OpenBabel::OBMol *molecule) = 0;
          virtual void read(std::ifstream &ifs) = 0;
          virtual ~CustomDataFunctor()
          {
          }
      };

      static std::size_t maxIndex()
      {
        return std::numeric_limits<std::size_t>::max();
      }



      /**
       * Construcor.
       *
       * @param filename The database filename.
       */
      Database(const std::string &filename);
      Database(const std::string &filename, CustomDataFunctor *customDataFunctor);
      Database(const std::string &filename, MessageHandler *messageHandler);
      Database(const std::string &filename, CustomDataFunctor *customDataFunctor, MessageHandler *messageHandler);

      ~Database();

      /**
       * Load the database. A database must always be loaded before it can be
       * used.
       */
      bool load(const std::vector<std::string> &fingerprints = std::vector<std::string>());

      /**
       * Get a list of all available fingerprints.
       */
      static std::vector<std::string> availableFingerprints();

      /**
       * Create a new database. A database must always be created before
       * molecules can be inserted.
       *
       * @param fingerprints The fingerprints to use in the database.
       */
      void create(const std::vector<std::string> &fingerprints);
      /**
       * @overload
       */
      void create(const std::string &fingerprint);
      
      void addFingerprint(const std::string &fingerprint);

      /**
       * Insert molecules in the database.
       */
      void insert(const std::string &filename, bool unique = true, int resumeInterval = 1000);


      /**
       * Get the number of records in the database.
       *
       * @return The number of record.
       */
      std::size_t numRecords() const;

      /**
       * Get the record with index @p index.
       *
       * @return The record with index @p index.
       */
      Record record(std::size_t index) const;


      /**
       * Get the list of fingerprints used in the database.
       */
      std::vector<std::string> databaseFingerprints() const;

      void databaseInformation() const;

      /**
       * Do an exact structure search. This is done by computing the canonical
       * smiles of the query molecule and comparing this to the canonical smiles
       * of the molecules in the database. This makes an exact structure search
       * much faster than a substructure search.
       *
       * @param query The query molecule to search for.
       *
       * @return True if the exact structure is found in the database.
       */
      std::size_t exactStructureSearch(OpenBabel::OBMol *query) const;
      /**
       * @overload
       */
      std::size_t exactStructureSearch(const std::string &smiles) const;
    
      /**
       * Do a fingerprint screen for the query molecule in the database. This
       * method is useful for evaluating the efficiency of the various
       * fingerprints. The lower the number of false positives (i.e. molecules
       * identified as match by the fingerprint but do not actually contain the
       * query as confirmed by an substructure isomorphism search), the better
       * the fingerprint is.
       *
       * @param query The query molecule.
       * @param fingerprints The fingerprints to use. These must be in the
       *        database.
       *
       * @return The record indices for all molecules matching the query.
       */
      std::vector<std::size_t> fingerprintScreen(OpenBabel::OBMol *query,
          const std::vector<std::string> &fingerprints, std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> fingerprintScreen(OpenBabel::OBMol *query,
          const std::string &fingerprint, std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> fingerprintScreen(const std::string &smiles,
          const std::vector<std::string> &fingerprints, std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> fingerprintScreen(const std::string &smiles,
          const std::string &fingerprint, std::size_t maxResults = 0) const;

      std::vector<std::pair<std::size_t, double> > similaritySearch(OpenBabel::OBMol *query,
          double tanimotoThreshold, const std::vector<std::string> &fingerprints,
          std::size_t maxResults = 0) const;
      std::vector<std::pair<std::size_t, double> > similaritySearch(OpenBabel::OBMol *query,
          double tanimotoThreshold, const std::string &fingerprint,
          std::size_t maxResults = 0) const;
      std::vector<std::pair<std::size_t, double> > similaritySearch(const std::string &smiles,
          double tanimotoThreshold, const std::vector<std::string> &fingerprints,
          std::size_t maxResults = 0) const;
      std::vector<std::pair<std::size_t, double> > similaritySearch(const std::string &smiles,
          double tanimotoThreshold, const std::string &fingerprints,
          std::size_t maxResults = 0) const;




 
      /**
       * Do a substructure isomorphism search for the entire database. This can
       * be very slow and is only useful for small databases and testing.
       *
       * @param query The query molecule.
       *
       * @return The record indices for all molecules matching the query.
       */
      std::vector<std::size_t> fullSubStructureSearch(OpenBabel::OBMol *query,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0, std::size_t startIndex = 0, 
          std::size_t stopIndex = Database::maxIndex()) const;
      /**
       * @overload
       */
      std::vector<std::size_t> fullSubStructureSearch(const std::string &smiles,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0) const;

      std::vector<std::size_t> fullSubStructureSearchThreaded(OpenBabel::OBMol *query,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2, int numThreads = 1,
          std::size_t maxResults = 0) const;
      std::vector<std::size_t> fullSubStructureSearchThreaded(const std::string &smiles,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2, int numThreads = 1,
          std::size_t maxResults = 0) const;
 
      /**
       * Do a substructure isomorphism search for all molecules matching the
       * fingerprint. With a good fingerprint this is the fastest method to
       * do a substructure search.
       *
       * @param query The query molecule.
       * @param fingerprints The fingerprints to use. These must be in the
       *        database.
       *
       * @return The record indices for all molecules matching the query.
       */
      std::vector<std::size_t> quickSubStructureSearch(OpenBabel::OBMol *query,
          const std::vector<std::string> &fingerprints,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> quickSubStructureSearch(OpenBabel::OBMol *query,
          const std::string &fingerprint,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> quickSubStructureSearch(const std::string &smiles,
          const std::vector<std::string> &fingerprints,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0) const;
      /**
       * @overload
       */
      std::vector<std::size_t> quickSubStructureSearch(const std::string &smiles,
          const std::string &fingerprints,
          enum SubStructureSearchAlgorithm algorithm = MolDBReducedGraphVF2,
          std::size_t maxResults = 0) const;



    private:
      friend class DatabaseImpl;

      /**
       * Database data pointer. This ensures binary compatibility.
       */
      DatabasePrivate * const d;
  };

  class TitleCustomDataFunctor : public Database::CustomDataFunctor
  {
    public:
      void write(std::fstream &ofs, OpenBabel::OBMol *molecule);
      void read(std::ifstream &ifs);
      std::string title(std::size_t index) const;
    private:
      StringVector m_titles;
  };

  class MolecularWeightCustomDataFunctor : public Database::CustomDataFunctor
  {
    public:
      void write(std::fstream &ofs, OpenBabel::OBMol *molecule);
      void read(std::ifstream &ifs);
      double weight(std::size_t index) const;
    private:
      std::vector<double> m_weights;
  };

  class CombineCustomDataFunctor : public Database::CustomDataFunctor
  {
    public:
      CombineCustomDataFunctor(const std::vector<Database::CustomDataFunctor*> &functors);
      void write(std::fstream &ofs, OpenBabel::OBMol *molecule);
      void read(std::ifstream &ifs);
    private:
      std::vector<CustomDataFunctor*> m_functors;
  };





}

#endif
