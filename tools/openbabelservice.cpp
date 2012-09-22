#include <boost/process.hpp>

#include "../src/openbabel.h"
#include "../src/database.h"

typedef boost::dynamic_bitset<> BitSet;

namespace bp = boost::process;

int main()
{
  OpenBabel::OBConversion conv;

  while (true) {
    std::string command;

    std::getline(std::cin, command);

    if (command == "EXIT")
      return 0;

    if (command.substr(0, 11) == "FINGERPRINT") {
      if (command.size() < 14) {
        std::cout << "ERROR Fingerprint name not found" << std::endl;
        continue;
      }

      std::size_t fingerprintNamePos = command.find(" ");
      if (fingerprintNamePos == std::string::npos) {
        std::cout << "ERROR Could not determine fingerprint name" << std::endl;
        continue;
      }
      
      std::size_t smilesPos = command.find(" ", fingerprintNamePos + 1);
      if (smilesPos == std::string::npos) {
        std::cout << "ERROR Could not determine smiles" << std::endl;
        continue;
      }
 
      std::string fingerprintName = command.substr(fingerprintNamePos + 1, smilesPos - fingerprintNamePos - 1);
      
      std::string smiles = command.substr(smilesPos + 1);

      if (fingerprintName.find("OpenBabel::") == std::string::npos) {
        std::cout << "ERROR Not a valid OpenBabel fingerprint name" << std::endl;
        continue;
      }

      OpenBabel::OBFingerprint *fp = OpenBabel::OBFingerprint::FindFingerprint(fingerprintName.substr(11, fingerprintName.size() - 11).c_str());
      if (!fp) {
        std::cout << "ERROR Fingerprint not found in OpenBabel" << std::endl;
        continue;
      }

      std::vector<unsigned int> fpbits;
      OpenBabel::OBMol *mol = MolDB::OB::File::readSmiles(smiles);
      fp->GetFingerprint(mol, fpbits);
      delete mol;

      BitSet::size_type numBits = fpbits.size() * fp->Getbitsperint();
      BitSet fpbitset(numBits);
      for (std::size_t i = 0; i < fpbits.size(); ++i)
        fpbitset |= (BitSet(numBits, fpbits[i]) << (fp->Getbitsperint() * i));

      std::cout << "FINGERPRINT " << fpbitset << std::endl;
    }

  }

  return 0;










  bp::context ctx;
  ctx.stdout_behavior = bp::capture_stream();

  std::string exec = "moldb_exe";
  std::vector<std::string> args;
  args.push_back("");

  bp::child c = bp::launch(exec, args, ctx);


  bp::pistream &is = c.get_stdout();
  std::string line;
  while (std::getline(is, line))
    std::cout << "Got line form child: " << line << std::endl;

  /*
  boost::process::status s = c.wait();

  if (s.exited() && s.exit_status() == EXIT_SUCCESS)
    std::cout << "Success" << std::endl;
  else
    std::cout << "Failed" << std::endl;
  std::cout << std::endl;
  */
}
