#include "../src/database.h"
#include "../src/openbabel.h"

#include "test.h"

using namespace MolDB;

void createInsertFile(const std::string &filename)
{
  std::ofstream ofs(filename.c_str());
  MOLDB_REQUIRE(ofs);

  ofs << "[Cl-].OC(=O)C(CS)[NH3+] 9976" << std::endl;
  ofs << "CC(C)CCC[C@@H](C)[C@H]1CCC2C3CC=C4CC(CC[C@@]4(C)C3CC[C@@]12C)OC(=O)c1ccccc1 9978" << std::endl;
  ofs << "CCCC(=O)OC1CC[C@]2(C)C(=CCC3C4CC[C@H]([C@H](C)CCCC(C)C)[C@]4(C)CCC23)C1 9980" << std::endl;
  ofs << "CCCCCCCCCC(=O)OC1CC[C@]2(C)C(=CCC3C2CC[C@@]2(C)[C@H](CCC32)[C@H](C)CCCC(C)C)C1 9982" << std::endl;
  /*
  ofs << "CCCCCC(=O)OC1CC[C@]2(C)C(=CCC3C2CC[C@@]2(C)[C@H](CCC32)[C@H](C)CCCC(C)C)C1 9984" << std::endl;
  ofs << "CCCCCCCC(=O)OC1CC[C@]2(C)C(=CCC3C2CC[C@@]2(C)[C@H](CCC32)[C@H](C)CCCC(C)C)C1 9986" << std::endl;
  ofs << "CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@H]3CC=C4CC(Cl)CC[C@@]4(C)[C@@H]3CC[C@@]12C 9988" << std::endl;
  ofs << "CC(C)CCC[C@@H](C)[C@H]1CCC2C3CC=C4CC(CC[C@@]4(C)C3CC[C@@]12C)OC(=O)Cl 9990" << std::endl;
  ofs << "CC(C)CCC[C@@H](C)[C@H]1CCC2C3CC=C4CC(CC[C@@]4(C)C3CC[C@@]12C)OC(=O)/C=C/c1ccccc1 9992" << std::endl;
  ofs << "O=COC1CC[C@]2(C)C(=CCC3C4CC[C@H]([C@H](C)CCCC(C)C)[C@]4(C)CCC23)C1 9994" << std::endl;
  */
}

unsigned int numRecords = 4;

int main()
{
  std::string databaseFilename = std::string(TESTDATADIR) + "test_database.moldb";
  std::string insertFilename = std::string(TESTDATADIR) + "test_database.smi";
  Database database1(databaseFilename);

  MOLDB_ASSERT(database1.availableFingerprints().size() != 0);

  // create the database
  database1.create(database1.availableFingerprints());
  MOLDB_COMPARE(database1.numRecords(), 0);
  database1.load();
  MOLDB_COMPARE(database1.databaseFingerprints().size(), database1.availableFingerprints().size());

  // insert molecules
  createInsertFile(insertFilename);
  database1.insert(insertFilename);
  MOLDB_COMPARE(database1.numRecords(), numRecords);
  bool delete_the_test_database_db_file = true;
  MOLDB_REQUIRE(database1.numRecords() == numRecords && delete_the_test_database_db_file);

  database1.databaseInformation();

//  for (std::size_t i = 0; i < database1.record(0).fingerprints.size(); ++i)
//    std::cout << database1.record(0).fingerprints[i].size() << std::endl;


  // creat a new database object for the same database file and load it
  Database database2(databaseFilename);
  database2.load();
  MOLDB_COMPARE(database1.numRecords(), database2.numRecords());
  // check a record
  MOLDB_COMPARE(database1.record(2).position, database2.record(2).position);
  MOLDB_COMPARE(database1.record(2).fingerprints.size(), database2.record(2).fingerprints.size());
  MOLDB_COMPARE(database1.record(2).bitcounts.size(), database2.record(2).bitcounts.size());
  MOLDB_COMPARE(database1.record(2).graph, database2.record(2).graph);
  MOLDB_COMPARE(database1.record(2).canonicalSmiles, database2.record(2).canonicalSmiles);

  MOLDB_COMPARE(database2.fingerprintScreen("CCC", "OpenBabel::FP2").size(), numRecords);
  MOLDB_COMPARE(database2.fingerprintScreen("CCC", "NoneExistantFingerprint").size(), 0);

  MOLDB_COMPARE(database2.fullSubStructureSearch(OB::File::readSmiles("C"), Database::MolDBRegularGraphVF2).size(), numRecords);
}
