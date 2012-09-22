#include <string>
#include <vector>
#include <limits>

namespace MolDB {

  std::string errors();
  bool open_database(const std::string &filename);
  std::size_t exact_structure_search(const std::string &smiles);
  std::vector<std::pair<std::size_t, double> > similarity_search(const std::string &smiles, double tanimotoThreshold);

  std::size_t invalidIndex()
  {
    return std::numeric_limits<std::size_t>::max();
  }

}

