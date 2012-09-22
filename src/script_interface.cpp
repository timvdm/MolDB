#include "script_interface.h"
#include "database.h"
#include "messagehandler.h"


namespace MolDB {

  std::stringstream* errorStream()
  {
    static std::stringstream *ss = 0;
    if (!ss)
      ss = new std::stringstream;
    return ss;
  }

  Database *database(const std::string &filename = std::string())
  {
    static Database *db = 0;
    if (!db)
      db = new Database(filename, new StdStreamMessageHandler(errorStream()));
    return db;
  }

  bool loadDatabase()
  {
    static bool loaded = false;
    if (!loaded) {
      loaded = database()->load();
      return loaded;
    }
    return true;
  }

  std::string errors()
  {
    std::string err = errorStream()->str();
    errorStream()->str("");
    return err;
  }
  
  bool open_database(const std::string &filename)
  {
    database(filename);
    return loadDatabase();
  }

  std::size_t exact_structure_search(const std::string &smiles)
  {
    return database()->exactStructureSearch(smiles);
  }
  
  std::vector<std::pair<std::size_t, double> > similarity_search(const std::string &smiles, double tanimotoThreshold)
  {
    return database()->similaritySearch(smiles, tanimotoThreshold, database()->databaseFingerprints());
  }

}

