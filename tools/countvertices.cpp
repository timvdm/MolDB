#include "../src/vertex.h"
#include "../src/edge.h"
#include "../src/graph.h"
#include "../src/openbabel.h"

#include <iostream>

using namespace MolDB;
using namespace OpenBabel;

int main(int argc, char **argv)
{
  OB::File file(argv[1]);

  OBMol *mol;
  unsigned long numAtoms = 0, numVertices = 0, count = 0;


  while ((mol = file.next())) {

    unsigned long vertexCount = 0;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (atom->GetValence() > 2)
        vertexCount++;
    }

    numVertices += vertexCount;
    numAtoms += mol->NumAtoms();
    count++;

    std::cout << count << ": " << mol->NumAtoms() << " -> " << vertexCount << "  " << numVertices << "/" << numAtoms << std::endl;

    delete mol;
  }

}
