#include "../src/timer.h"
#include "../src/openbabel.h"
#include "../src/isomorphism.h"
#include "../src/typetraits.h"
#include "../src/database.h"
#include "../src/util.h"

#include "../src/graph.h"
//#include "../src/tinygraph.h"

#include <openbabel/parsmart.h>

#ifdef HAVE_CHEMKIT
#include <chemkit/molecule.h>
#include <chemkit/substructurequery.h>
#endif

#include <iostream>
#include <map>
#include <algorithm>

using namespace MolDB;
using namespace OpenBabel;

typedef Graph<TypeTraits::EmptyType, unsigned char> GraphType;
typedef GraphType::VertexType VertexType;
typedef GraphType::EdgeType EdgeType;


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

void argument_error(const std::string &arg, const std::string &opt)
{
  std::cout << "Error: no " << opt << " for " << arg << " provided" << std::endl;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    //            /*****************************************************************************/
    std::cout << "MolDB version 0.1, Copyright(c) 2012, Tim Vandermeersch" << std::endl
              << "Usage: " << argv[0] << " [options] database" << std::endl
              << "Options:" << std::endl
              << "    -insert <filename>     Insert molecules from file" << std::endl
              << "    -query <smiles>        Query the database" << std::endl
              << "    -fingerprint-only      Only perform a fingerprint search for the query" << std::endl
              << "    -no-fingerprint        Do not use fingerprints for the query and do a full" << std::endl
              << "                           substructure search" << std::endl
              << "    -fingerprint <fp>      Select the fingerprint to use for queries [FP2, FP3," << std::endl
              << "                           FP4, MACCS], default is FP2" << std::endl
              << "    -start <index>         Start from molecule with index (works for insert)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << " * Do a full substructure search in a smiles file:" << std::endl
              << "   moldb -query \"N#Cc1ccccc1C#N\" molecules.smi" << std::endl
              << std::endl
              << " * Insert molecules from a smiles file in a MolDB database:" << std::endl
              << "   moldb -insert molecules.smi database.moldb" << std::endl
              << std::endl
              << "   This creates the MolDB database if needed." << std::endl
              << std::endl
              << " * Resume inserting molecules in a MolDB database starting with molecule 1001:" << std::endl
              << "   moldb -insert molecules.smi -start 1001 database.moldb" << std::endl
              << std::endl
              << " * Do a substructure search in a MolDB database:" << std::endl
              << "   moldb -query \"N#Cc1ccccc1C#N\" database.moldb" << std::endl
              << std::endl
              << "   A isomorphism substructure search is only performed for molecules that pass" << std::endl
              << "   the fingerprint test." << std::endl
              << std::endl
              << " * Do a substructure search in a MolDB database without using a fingerprint:" << std::endl
              << "   moldb -query \"N#Cc1ccccc1C#N\" -no-fingerprint database.moldb" << std::endl
              << std::endl
              ;
    return 0;
  }

  std::map<std::string, std::string> arguments;
  arguments["fingerprint"] = "OpenBabel::FP2";
  arguments["Smin"] = "0.8";
  arguments["no-fingerprint"] = "false";
  arguments["fingerprint-only"] = "false";
  //arguments["algorithm"] = "MolDBRegularGraphVF2";
  arguments["algorithm"] = "OpenBabelSMARTS";

  for (int i = 0; i < argc; ++i) {
    if (std::string(argv[i]) == "-insert") {
      if (i+1 >= argc) {
        argument_error("-insert", "filename");
        return 1;
      }
      arguments["insert"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-add-fingerprint") {
      if (i+1 >= argc) {
        argument_error("-add-fingerprint", "fingerprint");
        return 1;
      }
      arguments["add-fingerprint"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-query") {
      if (i+1 >= argc) {
        argument_error("-query", "smiles");
        return 1;
      }
      arguments["query"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-fingerprint") {
      if (i+1 >= argc) {
        argument_error("-fingerprint", "fp");
        return 1;
      }
      arguments["fingerprint"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-Smin") {
      if (i+1 >= argc) {
        argument_error("-Smin", "threshold");
        return 1;
      }
      arguments["Smin"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-start") {
      if (i+1 >= argc) {
        argument_error("-start", "index");
        return 1;
      }
      arguments["start"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-algorithm") {
      if (i+1 >= argc) {
        argument_error("-algorithm", "algorithm");
        return 1;
      }
      arguments["algorithm"] = argv[i+1];
      i++;
    } else if (std::string(argv[i]) == "-not-unique") {
      arguments["not-unique"] = "true";
    } else if (std::string(argv[i]) == "-fingerprint-only") {
      arguments["fingerprint-only"] = "true";
    } else if (std::string(argv[i]) == "-list-fingerprints") {
      arguments["list-fingerprints"] = "true";
    } else if (std::string(argv[i]) == "-create") {
      arguments["create"] = "true";
    } else if (std::string(argv[i]) == "-no-fingerprint") {
      arguments["no-fingerprint"] = "true";
    } else if (std::string(argv[i]) == "-info") {
      arguments["info"] = "true";
    } else {
      arguments["database"] = argv[i];
    }

  }


  //
  // list the fingerprints
  //
  if (arguments.find("list-fingerprints") != arguments.end()) {
    std::vector<std::string> fingerprints = Database::availableFingerprints();
    for (std::size_t i = 0; i < fingerprints.size(); ++i)
      std::cout << fingerprints[i] << std::endl;
    return 0;
  }
  
  Database database(arguments["database"]);

  //
  // create a new database
  //
  if (arguments.find("create") != arguments.end()) {
    std::vector<std::string> fingerprints = Impl::tokenize(arguments["fingerprint"], ",");
    database.create(fingerprints);
  }
 
  //
  // add new fingerprint
  //
  if (arguments.find("add-fingerprint") != arguments.end()) {
    database.load();
    database.addFingerprint(arguments["add-fingerprint"]);
  }
  
  //
  // print database information
  //
  if (arguments.find("info") != arguments.end()) {
    database.load();
    database.databaseInformation();
  }

  //
  // insert new molecules in the database
  //
  if (arguments.find("insert") != arguments.end()) {
    database.load();
    bool unique = (arguments.find("not-unique") == arguments.end()) ? true : false;
    database.insert(arguments["insert"], unique);
  }

  //
  // query the database
  //
  if (arguments.find("query") != arguments.end()) {
    std::vector<std::string> fingerprints = Impl::tokenize(arguments["fingerprint"], ",");
    database.load(fingerprints);

    enum Database::SubStructureSearchAlgorithm algorithm = Database::InvalidSubStructureSearchAlgorithm;
    if (arguments["algorithm"] == "MolDBReducedGraphVF2")
      algorithm = Database::MolDBReducedGraphVF2;
    else if (arguments["algorithm"] == "MolDBRegularGraphVF2")
      algorithm = Database::MolDBRegularGraphVF2;
    else if (arguments["algorithm"] == "OpenBabelVF2")
      algorithm = Database::OpenBabelVF2;
    else if (arguments["algorithm"] == "OpenBabelSMARTS")
      algorithm = Database::OpenBabelSMARTS;
    else if (arguments["algorithm"] == "ChemKitVF2")
      algorithm = Database::ChemKitVF2;
    else if (arguments["algorithm"] == "RDKitVF2")
      algorithm = Database::RDKitVF2;

    //std::vector<std::size_t> hits = database.fullSubStructureSearch(arguments["query"], algorithm);
    //std::vector<std::size_t> hits = database.fullSubStructureSearchThreaded(arguments["query"], algorithm, 4);
    //std::cout << hits.size() << " hits" << std::endl;
    
    database.databaseInformation();
    
    //database.exactStructureSearch(arguments["query"]);
    
    //database.fingerprintScreen(arguments["query"], "OpenBabel::FP2");


    /*
    const std::vector<std::string> &fingerprints = database.databaseFingerprints();
    if (std::find(fingerprints.begin(), fingerprints.end(), arguments["fingerprint"]) == fingerprints.end())
      arguments["fingerprint"] = fingerprints[0];
    */


    // find all hits
    /*
    { 
      //std::ofstream hitSmiles("hits.smi");

      int count = 0;
      std::vector<std::vector<std::size_t> > hits;
      std::ifstream ifs(arguments["query"].c_str());
      std::string line;
      while (std::getline(ifs, line)) {
        std::string smiles = line.substr(0, line.find(" "));
        if (smiles.size() > 40)
          continue;
        hits.push_back(database.quickSubStructureSearch(smiles, fingerprints, algorithm));
        //for (std::size_t i = 0; i < hits.back().size(); ++i)
        //  hitSmiles << database.record(hits.back()[i]).canonicalSmiles << "\t" << i << std::endl;

        if ((++count % 10000) == 0)
          break;
        std::cout << count << std::endl;
      }

      std::ofstream ofs("full.hits");
      for (std::size_t i = 0; i < hits.size(); ++i) {
        for (std::size_t j = 0; j < hits[i].size(); ++j)
          ofs << hits[i][j] << " ";
        ofs << std::endl;
      }
    
      return 0;
    }
    */

    // load hits
    /*
    std::size_t allHitsSize = 0;
    std::vector<std::vector<std::size_t> > allHits;
    {
      std::ifstream ifs("full.hits");
      std::string line;
      while (std::getline(ifs, line)) {
        allHits.resize(allHits.size() + 1);
        std::size_t pos = 0;
        while (true) {
          std::size_t next_pos = line.find(" ", pos);
          if (next_pos == std::string::npos)
            break;
          std::stringstream ss(line.substr(pos, next_pos - pos));
          std::size_t index;
          ss >> index;
          //std::cout << index << std::endl;
          allHits.back().push_back(index);
          pos = next_pos + 1;
        }
        std::sort(allHits.back().begin(), allHits.back().end());
        allHitsSize += allHits.back().size();
      }
    }

    std::ofstream queries("queries.smi");


    std::vector<int> screens, hits;

    std::size_t screenedSize = 0;
    int count = 0;
    std::ifstream ifs(arguments["query"].c_str());
    std::string line;
    while (std::getline(ifs, line)) {
      std::string smiles = line.substr(0, line.find(" "));
      if (smiles.size() > 40)
        continue;

      queries << smiles << std::endl;
      
      std::vector<std::size_t> screened = database.fingerprintScreen(smiles, fingerprints);
      screens.push_back(screened.size());

      screenedSize += screened.size();
      std::sort(screened.begin(), screened.end());
      assert(screened.size() >= allHits[count].size());

      std::vector<std::size_t> intersection;
      std::set_intersection(allHits[count].begin(), allHits[count].end(), screened.begin(), screened.end(), std::back_inserter(intersection));
      std::cout << smiles << std::endl;
      std::cout << "correct: " << allHits[count].size() << std::endl;
      std::cout << "current: " << intersection.size() << std::endl;
      if (intersection.size() != allHits[count].size()) {
        for (std::size_t i = 0; i < allHits[count].size(); ++i)
          if (!std::binary_search(screened.begin(), screened.end(), allHits[count][i]))
            std::cout << database.record(allHits[count][i]).canonicalSmiles << std::endl;
        assert(intersection.size() == allHits[count].size());
      }
      for (std::size_t i = 0; i < intersection.size(); ++i)
        assert(intersection[i] == allHits[count][i]);

      //if ((++count % 1000) == 0) {
      if ((++count % 17) == 0) {
        std::string names;
        for (std::size_t i = 0; i < fingerprints.size(); ++i)
          names += convertFingerprintName(fingerprints[i]);
        std::ofstream ofs(make_string(names, "_", count, ".screen").c_str());
        for (std::size_t i = 0; i < screens.size(); ++i)
          ofs << screens[i] << std::endl;
      }
      //if ((count % 10000) == 0) {
      if ((count % 17) == 0) {
        assert(allHitsSize <= screenedSize);
        std::cout << allHitsSize << " <= " << screenedSize << std::endl;
        return 0;
      }
      std::cout << count << std::endl;
    }
  */

    int quick = 0;
    int full = 0;
    int screen = 1;
    int similarity = 0;

    if (full) {
      int hits = 0;
      hits += database.fullSubStructureSearch("ONC1CC(C(O)C1O)[n]2cnc3c(NC4CC4)ncnc23", algorithm).size();
      hits += database.fullSubStructureSearch("Nc1ncnc2[n]cnc12", algorithm).size();
      hits += database.fullSubStructureSearch("CNc1ncnc2[n](C)cnc12", algorithm).size();
      hits += database.fullSubStructureSearch("Nc1ncnc2[n](cnc12)C3CCCC3", algorithm).size();
      hits += database.fullSubStructureSearch("CC12CCC3C(CCC4=CC(O)CCC34C)C1CCC2", algorithm).size();
      hits += database.fullSubStructureSearch("OC2=CC(=O)c1c(cccc1)O2", algorithm).size();
      hits += database.fullSubStructureSearch("Nc1nnc(S)s1", algorithm).size();
      hits += database.fullSubStructureSearch("C1C2SCCN2C1", algorithm).size();
      hits += database.fullSubStructureSearch("CP(O)(O)=O", algorithm).size();
      hits += database.fullSubStructureSearch("CCCCCP(O)(O)=O", algorithm).size();
      hits += database.fullSubStructureSearch("N2CCC13CCCCC1C2Cc4c3cccc4", algorithm).size();
      hits += database.fullSubStructureSearch("s1cncc1", algorithm).size();
      hits += database.fullSubStructureSearch("C34CCC1C(CCC2CC(=O)CCC12)C3CCC4", algorithm).size();
      hits += database.fullSubStructureSearch("CCCCCCCCCCCP(O)(O)=O", algorithm).size();
      hits += database.fullSubStructureSearch("CC1CCCC1", algorithm).size();
      hits += database.fullSubStructureSearch("CCC1CCCC1", algorithm).size();
      hits += database.fullSubStructureSearch("CCCC1CCCC1", algorithm).size();
      hits += database.quickSubStructureSearch("C(F)(F)F", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C[Mg]Cl", fingerprints, algorithm).size();
      std::cout << "# Hits: " << hits << std::endl;
    }

    if (quick) {
      int hits = 0;
      hits += database.quickSubStructureSearch("ONC1CC(C(O)C1O)[n]2cnc3c(NC4CC4)ncnc23", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1ncnc2[n]cnc12", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CNc1ncnc2[n](C)cnc12", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1ncnc2[n](cnc12)C3CCCC3", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CC12CCC3C(CCC4=CC(O)CCC34C)C1CCC2", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("OC2=CC(=O)c1c(cccc1)O2", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1nnc(S)s1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C1C2SCCN2C1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CP(O)(O)=O", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCCCCP(O)(O)=O", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("N2CCC13CCCCC1C2Cc4c3cccc4", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("s1cncc1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C34CCC1C(CCC2CC(=O)CCC12)C3CCC4", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCCCCCCCCCCP(O)(O)=O", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CC1CCCC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCC1CCCC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCCC1CCCC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C(F)(F)F", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C[Mg]Cl", fingerprints, algorithm).size();
      return 0;
      
      hits += database.quickSubStructureSearch("c1ncnc2c1[nH]cn2", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1ccc(cc1)Cl", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1ccc(cc1)C=O", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1(c2c(cnc1)c(cnn2)N1CC1)Cl", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1(c2c(nnc1)c(cnc2)c1ccsc1)N1CC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1(c2c(cc(c1)S(=O)(=O)O)cc(c(c2O)C=C)S(=O)(=O)O)N", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1(c2c(cc(c1N=N)S(=O)(=O)O)cc(cc2)S(=O)(=O)O)O", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1cc(ccc1)C1NOC(=O)N1c1ccccc1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("n1cnc(c2c1[nH]cn2)SCc1ccccc1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("c1c2c(ccc1OC)nc1c(c2)ccc(c1)Cl", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("n1c(cc(nc1NS(=O)(=O)c1ccccc1)C)C", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("n1c2c(c(cc1)N)ccc(c2)Cl", fingerprints, algorithm).size();

      hits += database.quickSubStructureSearch("CSC(=S)NNC(=O)c1ccccc1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCN(CC)c1cc(N)ccc1C", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Cc1nc2CCCc2c2CCCc12", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1cc(I)c(I)c(I)c1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Cc1cccc2cc3CNCCOc3nc12", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CC1CNCc2cc3ccc(C)c(C)c3nc2O1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("NCc1cc(O)nc(n1)N1CCCC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Cc1nn(Cc2ccccc2)c(N)c1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C=CCC(=O)N(C)CC1(CO)CCCC1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1nnc(SCCOc2ccccc2)[nH]1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Cc1ccc(C#N)c(S)n1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("O=C1OCC1(C)C", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("C[Sn](C)(F)F", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("Nc1ccnc(F)c1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCOc1ccncc1Br", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("O=Cc1c[nH]nc1", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("N#[Cr]", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("CCCCC#CCCC#C", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("ClCCC(C)CCCl", fingerprints, algorithm).size();
      hits += database.quickSubStructureSearch("O=C1S/C(=C/c2ccco2)/C(=O)N1C", fingerprints, algorithm).size();
      std::cout << "# Hits: " << hits << std::endl;
    }

    if (screen) {
      Timer timer;
      database.fingerprintScreen("ONC1CC(C(O)C1O)[n]2cnc3c(NC4CC4)ncnc23", fingerprints);
      database.fingerprintScreen("Nc1ncnc2[n]cnc12", fingerprints);
      database.fingerprintScreen("CNc1ncnc2[n](C)cnc12", fingerprints);
      database.fingerprintScreen("Nc1ncnc2[n](cnc12)C3CCCC3", fingerprints);
      database.fingerprintScreen("CC12CCC3C(CCC4=CC(O)CCC34C)C1CCC2", fingerprints);
      database.fingerprintScreen("OC2=CC(=O)c1c(cccc1)O2", fingerprints);
      database.fingerprintScreen("Nc1nnc(S)s1", fingerprints);
      database.fingerprintScreen("C1C2SCCN2C1", fingerprints);
      database.fingerprintScreen("CP(O)(O)=O", fingerprints);
      database.fingerprintScreen("CCCCCP(O)(O)=O", fingerprints);
      database.fingerprintScreen("N2CCC13CCCCC1C2Cc4c3cccc4", fingerprints);
      database.fingerprintScreen("s1cncc1", fingerprints);
      database.fingerprintScreen("C34CCC1C(CCC2CC(=O)CCC12)C3CCC4", fingerprints);
      database.fingerprintScreen("CCCCCCCCCCCP(O)(O)=O", fingerprints);
      database.fingerprintScreen("CC1CCCC1", fingerprints);
      database.fingerprintScreen("CCC1CCCC1", fingerprints);
      database.fingerprintScreen("CCCC1CCCC1", fingerprints);
      
      return 0;
      database.fingerprintScreen("C(F)(F)F", fingerprints);
      database.fingerprintScreen("C[Mg]Cl", fingerprints);

      database.fingerprintScreen("c1ncnc2c1[nH]cn2", fingerprints);
      database.fingerprintScreen("c1ccc(cc1)Cl", fingerprints);
      database.fingerprintScreen("c1ccc(cc1)C=O", fingerprints);
      database.fingerprintScreen("c1(c2c(cnc1)c(cnn2)N1CC1)Cl", fingerprints);
      database.fingerprintScreen("c1(c2c(nnc1)c(cnc2)c1ccsc1)N1CC1", fingerprints);
      database.fingerprintScreen("c1(c2c(cc(c1)S(=O)(=O)O)cc(c(c2O)C=C)S(=O)(=O)O)N", fingerprints);
      database.fingerprintScreen("c1(c2c(cc(c1N=N)S(=O)(=O)O)cc(cc2)S(=O)(=O)O)O", fingerprints);
      database.fingerprintScreen("c1cc(ccc1)C1NOC(=O)N1c1ccccc1", fingerprints);
      database.fingerprintScreen("n1cnc(c2c1[nH]cn2)SCc1ccccc1", fingerprints);
      database.fingerprintScreen("c1c2c(ccc1OC)nc1c(c2)ccc(c1)Cl", fingerprints);
      database.fingerprintScreen("n1c(cc(nc1NS(=O)(=O)c1ccccc1)C)C", fingerprints);
      database.fingerprintScreen("n1c2c(c(cc1)N)ccc(c2)Cl", fingerprints);
      database.fingerprintScreen("CSC(=S)NNC(=O)c1ccccc1", fingerprints);
      database.fingerprintScreen("CCN(CC)c1cc(N)ccc1C", fingerprints);
      database.fingerprintScreen("Cc1nc2CCCc2c2CCCc12", fingerprints);
      database.fingerprintScreen("Nc1cc(I)c(I)c(I)c1", fingerprints);
      database.fingerprintScreen("Cc1cccc2cc3CNCCOc3nc12", fingerprints);
      database.fingerprintScreen("CC1CNCc2cc3ccc(C)c(C)c3nc2O1", fingerprints);
      database.fingerprintScreen("NCc1cc(O)nc(n1)N1CCCC1", fingerprints);
      database.fingerprintScreen("Cc1nn(Cc2ccccc2)c(N)c1", fingerprints);
      database.fingerprintScreen("C=CCC(=O)N(C)CC1(CO)CCCC1", fingerprints);
      database.fingerprintScreen("Nc1nnc(SCCOc2ccccc2)[nH]1", fingerprints);
      database.fingerprintScreen("Cc1ccc(C#N)c(S)n1", fingerprints);
      database.fingerprintScreen("O=C1OCC1(C)C", fingerprints);
      database.fingerprintScreen("C[Sn](C)(F)F", fingerprints);
      database.fingerprintScreen("Nc1ccnc(F)c1", fingerprints);
      database.fingerprintScreen("CCOc1ccncc1Br", fingerprints);
      database.fingerprintScreen("O=Cc1c[nH]nc1", fingerprints);
      database.fingerprintScreen("N#[Cr]", fingerprints);
      database.fingerprintScreen("CCCCC#CCCC#C", fingerprints);
      database.fingerprintScreen("ClCCC(C)CCCl", fingerprints);
      database.fingerprintScreen("O=C1S/C(=C/c2ccco2)/C(=O)N1C", fingerprints);
      std::cout << "elapsed: " << timer.elapsed() << " seconds" << std::endl;
    }

    if (similarity) {
      std::stringstream ss(arguments["Smin"]);
      double Smin;
      ss >> Smin;

      Timer timer;
      database.similaritySearch("ONC1CC(C(O)C1O)[n]2cnc3c(NC4CC4)ncnc23", Smin, fingerprints);
      database.similaritySearch("Nc1ncnc2[n]cnc12", Smin, fingerprints);
      database.similaritySearch("CNc1ncnc2[n](C)cnc12", Smin, fingerprints);
      database.similaritySearch("Nc1ncnc2[n](cnc12)C3CCCC3", Smin, fingerprints);
      database.similaritySearch("CC12CCC3C(CCC4=CC(O)CCC34C)C1CCC2", Smin, fingerprints);
      database.similaritySearch("OC2=CC(=O)c1c(cccc1)O2", Smin, fingerprints);
      database.similaritySearch("Nc1nnc(S)s1", Smin, fingerprints);
      database.similaritySearch("C1C2SCCN2C1", Smin, fingerprints);
      database.similaritySearch("CP(O)(O)=O", Smin, fingerprints);
      database.similaritySearch("CCCCCP(O)(O)=O", Smin, fingerprints);
      database.similaritySearch("N2CCC13CCCCC1C2Cc4c3cccc4", Smin, fingerprints);
      database.similaritySearch("s1cncc1", Smin, fingerprints);
      database.similaritySearch("C34CCC1C(CCC2CC(=O)CCC12)C3CCC4", Smin, fingerprints);
      database.similaritySearch("CCCCCCCCCCCP(O)(O)=O", Smin, fingerprints);
      database.similaritySearch("CC1CCCC1", Smin, fingerprints);
      database.similaritySearch("CCC1CCCC1", Smin, fingerprints);
      database.similaritySearch("CCCC1CCCC1", Smin, fingerprints);
      
      database.similaritySearch("c1ncnc2c1[nH]cn2", Smin, fingerprints);
      database.similaritySearch("c1ccc(cc1)Cl", Smin, fingerprints);
      database.similaritySearch("c1ccc(cc1)C=O", Smin, fingerprints);
      database.similaritySearch("c1(c2c(cnc1)c(cnn2)N1CC1)Cl", Smin, fingerprints);
      database.similaritySearch("c1(c2c(nnc1)c(cnc2)c1ccsc1)N1CC1", Smin, fingerprints);
      database.similaritySearch("c1(c2c(cc(c1)S(=O)(=O)O)cc(c(c2O)C=C)S(=O)(=O)O)N", Smin, fingerprints);
      database.similaritySearch("c1(c2c(cc(c1N=N)S(=O)(=O)O)cc(cc2)S(=O)(=O)O)O", Smin, fingerprints);
      database.similaritySearch("c1cc(ccc1)C1NOC(=O)N1c1ccccc1", Smin, fingerprints);
      database.similaritySearch("n1cnc(c2c1[nH]cn2)SCc1ccccc1", Smin, fingerprints);
      database.similaritySearch("c1c2c(ccc1OC)nc1c(c2)ccc(c1)Cl", Smin, fingerprints);
      database.similaritySearch("n1c(cc(nc1NS(=O)(=O)c1ccccc1)C)C", Smin, fingerprints);
      database.similaritySearch("n1c2c(c(cc1)N)ccc(c2)Cl", Smin, fingerprints);
      database.similaritySearch("CSC(=S)NNC(=O)c1ccccc1", Smin, fingerprints);
      database.similaritySearch("CCN(CC)c1cc(N)ccc1C", Smin, fingerprints);
      database.similaritySearch("Cc1nc2CCCc2c2CCCc12", Smin, fingerprints);
      database.similaritySearch("Nc1cc(I)c(I)c(I)c1", Smin, fingerprints);
      database.similaritySearch("Cc1cccc2cc3CNCCOc3nc12", Smin, fingerprints);
      database.similaritySearch("CC1CNCc2cc3ccc(C)c(C)c3nc2O1", Smin, fingerprints);
      database.similaritySearch("NCc1cc(O)nc(n1)N1CCCC1", Smin, fingerprints);
      database.similaritySearch("Cc1nn(Cc2ccccc2)c(N)c1", Smin, fingerprints);
      database.similaritySearch("C=CCC(=O)N(C)CC1(CO)CCCC1", Smin, fingerprints);
      database.similaritySearch("Nc1nnc(SCCOc2ccccc2)[nH]1", Smin, fingerprints);
      database.similaritySearch("Cc1ccc(C#N)c(S)n1", Smin, fingerprints);
      database.similaritySearch("O=C1OCC1(C)C", Smin, fingerprints);
      database.similaritySearch("C[Sn](C)(F)F", Smin, fingerprints);
      database.similaritySearch("Nc1ccnc(F)c1", Smin, fingerprints);
      database.similaritySearch("CCOc1ccncc1Br", Smin, fingerprints);
      database.similaritySearch("O=Cc1c[nH]nc1", Smin, fingerprints);
      database.similaritySearch("N#[Cr]", Smin, fingerprints);
      database.similaritySearch("CCCCC#CCCC#C", Smin, fingerprints);
      database.similaritySearch("ClCCC(C)CCCl", Smin, fingerprints);
      database.similaritySearch("O=C1S/C(=C/c2ccco2)/C(=O)N1C", Smin, fingerprints);
      database.similaritySearch("C(F)(F)F", Smin, fingerprints);
      database.similaritySearch("C[Mg]Cl", Smin, fingerprints);
  
      std::cout << "elapsed: " << timer.elapsed() << " seconds" << std::endl;
    }



    return 0;














    // read the query smiles string
    OBMol *query = OB::File::readSmiles(arguments["query"]);

    if (query) {
    
      // create the query reduced graph
      GraphType *queryGraph = (arguments["no-reducedgraph"] == "true") ? OB::graph<GraphType>(query) : OB::reducedGraph<GraphType>(query);
      std::cout << queryGraph->toString() << std::endl;

      // initialize the fingerprints only if needed
      //bool useFP = (arguments["no-fingerprint"] == "false");
      //if (arguments["fingerprint"] == "FP2" && useFP)

      if (arguments["database"].find(".graphs") != std::string::npos) {
        std::cout << "Database " << arguments["database"] << " is a graph file" << std::endl;
        std::ifstream ifs(arguments["database"].c_str());
        if (!ifs) {
          std::cout << "Could not open file " << arguments["database"] << std::endl;
          return 1;
        }

        std::vector<GraphType*> graphs;

        std::string line;
        while (ifs.good()) {

          std::getline(ifs, line);
          if (line.empty())
            break;
          GraphType *queriedGraph = GraphType::fromString(line);

          GraphIsomorphismVF2<true, GraphType> iso;
          bool found = iso.isSubgraph(queryGraph, queriedGraph);

          graphs.push_back(queriedGraph);
        }

        double totalTime = 0.0;
        for (int i = 0; i < 20; ++i) {
          Timer timer;
          for (std::size_t i = 0; i < graphs.size(); ++i) {
            GraphIsomorphismVF2<true, GraphType> iso;
            bool found = iso.isSubgraph(queryGraph, graphs[i]);
          }
          totalTime += timer.elapsed();
          std::cout << timer.elapsed() << " seconds" << std::endl;
        }
        std::cout << "average over 20 runs: " << totalTime / 20 << " seconds" << std::endl;

        //std::cout << "File I/O " << ioTimer.elapsed() << " seconds, isomorphism " << isoTimer.elapsed() << " seconds" << std::endl;


      } else if (OB::File::isValidFileFormat(arguments["database"])) {
        // search in a molecule file format
        OB::File file(arguments["database"].c_str());

        OBMol *queried;
        unsigned long count = 0;

        std::map<unsigned int, std::pair<unsigned int, double> > times;
        double time = 0.0;
        while ((queried = file.next())) {
          count++;
          //if (count % 1000 == 0)
            std::cout << "Checking molecule #" << count << std::endl;

          GraphType *queriedGraph = (arguments["no-reducedgraph"] == "true") ? OB::graph<GraphType>(queried) : OB::reducedGraph<GraphType>(queried);

          std::vector<OpenBabel::OBAtom*> atoms;
          for (std::size_t i = 0; i < queried->NumAtoms(); ++i)
            atoms.push_back(queried->GetAtom(i + 1));
          std::random_shuffle(atoms.begin(), atoms.end());

          //std::cout << queryGraph->toString() << std::endl;
          //std::cout << queriedGraph->toString() << std::endl;

          Timer timer;
          if (arguments.find("smarts") != arguments.end()) {
            if (count == 1) {
              delete queryGraph;
              queryGraph = 0;
            }
            query = new OpenBabel::OBMol(*queried);
            queried->RenumberAtoms(atoms);
            OpenBabel::OBSmartsPattern sp;
//            sp.Init(arguments["query"]);
            std::string smilesQuery = OB::File::writeSmiles(query);
            smilesQuery = smilesQuery.substr(0, smilesQuery.find('\t'));
            sp.Init(smilesQuery);
            //std::cout << smilesQuery << std::endl;
            bool found = sp.Match(*queried);
            //std::cout << "MATCH = " << found << std::endl;
            delete query;
            query = 0;
          } else if (arguments.find("chemkit") != arguments.end()) {
#ifdef HAVE_CHEMKIT
            if (count == 1) {
              delete queryGraph;
              queryGraph = 0;
            }
            query = new OpenBabel::OBMol(*queried);
            queried->RenumberAtoms(atoms);

            boost::shared_ptr<chemkit::Molecule> chemkitQuery(chemKitMolecule(query));
            boost::shared_ptr<chemkit::Molecule> chemkitQueried(chemKitMolecule(queried));

            chemkit::SubstructureQuery q(chemkitQuery);
            bool found = q.matches(chemkitQueried.get());
            //std::cout << "MATCH = " << found << std::endl;
            delete query;
            query = 0;
#endif
          } else {
            queried->RenumberAtoms(atoms);
            delete queryGraph;
            //std::cout << queriedGraph->toString() << std::endl << std::endl;
            queryGraph = (arguments["no-reducedgraph"] == "true") ? OB::graph<GraphType>(queried) : OB::reducedGraph<GraphType>(queried);
            //std::cout << queryGraph->toString() << std::endl << std::endl;

            bool found;
            if (arguments["no-reducedgraph"] == "true") {
              GraphIsomorphismVF2<false, GraphType> iso;
              found = iso.isSubgraph(queryGraph, queriedGraph);
            } else {
              GraphIsomorphismVF2<true, GraphType> iso;
              found = iso.isSubgraph(queryGraph, queriedGraph);
            }
            //std::cout << "MATCH = " << found << std::endl;
          }

          times[queried->NumAtoms()].first++;
          times[queried->NumAtoms()].second += timer.elapsed();
          time += timer.elapsed();
          //std::cout << "found isomorphism = " << found << std::endl;

          delete queried;
          delete queriedGraph;
        }

        for (std::map<unsigned int, std::pair<unsigned int, double> >::const_iterator i = times.begin(); i != times.end(); ++i)
          std::cout << i->first << "," << i->second.first << "," << i->second.second << std::endl;
        std::cout << "time = " << time << " seconds" << std::endl;
      } else {
        // search in a MolDB database file

      }

      if (query)
        delete query;
      if (queryGraph)
        delete queryGraph;
    }
  }


  return 0;
}
