#include "../src/timer.h"
#include "../src/openbabel.h"

#include <iostream>
#include <map>

using namespace MolDB;
using namespace OpenBabel;

typedef Graph<TypeTraits::EmptyType, unsigned char> GraphType;
typedef GraphType::VertexType VertexType;
typedef GraphType::EdgeType EdgeType;





int main(int argc, char **argv)
{
  OB::File file(argv[1]);

  OBMol *mol;
  unsigned long count = 0;

  std::map<unsigned int, unsigned long> vertexCounts;
  std::map<unsigned int, unsigned long> vertexDegreeCounts;
  while ((mol = file.next())) {
    count++;
    if ((count % 10000) == 0)
      std::cout << count << std::endl;

    vertexCounts[mol->NumAtoms()]++;
    FOR_ATOMS_OF_MOL (atom, mol)
      vertexDegreeCounts[atom->GetValence()]++;

    delete mol;
  }
  
  std::cout << "VertexCounts" << std::endl;
  for (std::map<unsigned int, unsigned long>::const_iterator i = vertexCounts.begin(); i != vertexCounts.end(); ++i)
    std::cout << i->first << "," << i->second << std::endl;
  std::cout << "VertexDegreeCounts" << std::endl;
  for (std::map<unsigned int, unsigned long>::const_iterator i = vertexDegreeCounts.begin(); i != vertexDegreeCounts.end(); ++i)
    std::cout << i->first << "," << i->second << std::endl;



  return 0;
}
