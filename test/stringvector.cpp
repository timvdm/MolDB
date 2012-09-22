#include "../src/stringvector.h"

#include "test.h"

int main()
{
  int input;

  /*
  std::vector<std::string> stdsv;
  for (std::size_t i = 0; i < 20000000; ++i)
    stdsv.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

  std::cout << "std::vector<std::string> with 20M strings of 40 characters." << std::endl;
  std::cin >> input;

  stdsv.clear();
  std::cout << "clearing std::vector" << std::endl;
  std::cin >> input;

  MolDB::StringVector moldbsv;
  moldbsv.reserve(20000000, 50);
  for (std::size_t i = 0; i < 20000000; ++i)
    moldbsv.push_back("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

  std::cout << "MolDB::StringVector with 20M strings of 40 characters." << std::endl;
  std::cin >> input;

  return 0;
  */


  MolDB::StringVector sv1;
  sv1.push_back("AAA");
  sv1.push_back("BBB");
  sv1.push_back("CCC");


  MolDB::StringVector sv;
  MOLDB_COMPARE(sv.size(), 0);
  sv.push_back("test1");
  MOLDB_COMPARE(sv.size(), 1);
  MOLDB_COMPARE(std::string(sv[0]), std::string("test1"));
  sv.push_back("test2");
  MOLDB_COMPARE(std::string(sv[1]), std::string("test2"));
  sv.push_back("test3");
  MOLDB_COMPARE(std::string(sv[2]), std::string("test3"));
  sv.push_back("test4");
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(std::string(sv[3]), std::string("test4"));
  MOLDB_COMPARE(sv.size(), 4);
  sv.clear();
  MOLDB_COMPARE(sv.size(), 0);


  for (int i = 0; i < 10000; ++i) {
    sv.push_back("fdhsdfhsdlkghlksdhglkfdshglhs");
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
  }
  
  for (int i = 0; i < 10000; ++i) {
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
  }

  for (int i = 0; i < 10000; ++i) {
    sv.push_back("fdhsdfhsdlkghlksdhglkfdshglhs");
    MOLDB_REQUIRE(std::string(sv[i]) == std::string("fdhsdfhsdlkghlksdhglkfdshglhs"));
    const char * const p = sv[i];
  }

}
