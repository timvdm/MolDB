#include "../src/fingerprintindex.h"

#include <iostream>

#include "test.h"


using namespace MolDB;

unsigned long numNodes(int numBits)
{
  if (numBits == 1)
    return 3;

  unsigned long sum = 3;
  for (int i = 0; i < numBits - 1; ++i)
    sum = 2 * sum + 1;

  return sum;
};

unsigned long numLeafNodes(int numBits)
{
  unsigned long sum = 2;
  for (int i = 0; i < numBits - 1; ++i)
    sum *= 2;
  return sum;
};



int main()
{
  /*
  StdOutMessageHandler msg;
  InMemoryFingerprintStorage storage("test_database", "OpenBabel::FP2", 5);
  CompressedBinarySearchTreeFingerprintIndex<InMemoryFingerprintStorage> index("test_database", "OpenBabel::FP2", 5, storage);
  typedef CompressedBinarySearchTreeFingerprintIndex<InMemoryFingerprintStorage>::LeafNode LeafNode;


  Fingerprint fp00000(5, 0);
  index.createLeaf(fp00000).push_back(1);
  MOLDB_COMPARE(index.createLeaf(fp00000)[0], 1);

  Fingerprint fp01000(5, 8);
  index.createLeaf(fp01000).push_back(2);
  MOLDB_COMPARE(index.createLeaf(fp00000)[0], 1);
  MOLDB_COMPARE(index.createLeaf(fp01000)[0], 2);

  Fingerprint fp10111(5, 23);
  index.createLeaf(fp10111).push_back(3);
  MOLDB_COMPARE(index.createLeaf(fp00000)[0], 1);
  MOLDB_COMPARE(index.createLeaf(fp01000)[0], 2);
  MOLDB_COMPARE(index.createLeaf(fp10111)[0], 3);

  Fingerprint fp01111(5, 15);
  index.createLeaf(fp01111).push_back(4);
  MOLDB_COMPARE(index.createLeaf(fp00000)[0], 1);
  MOLDB_COMPARE(index.createLeaf(fp01000)[0], 2);
  MOLDB_COMPARE(index.createLeaf(fp10111)[0], 3);
  MOLDB_COMPARE(index.createLeaf(fp01111)[0], 4);
*/
  int depth = 0;
  /*
  std::cout << index.findLeaf(fp00000, 0, depth)->fingerprints[0] << std::endl;
  depth = 0;
  std::cout << index.findLeaf(fp01000, 0, depth)->fingerprints[0] << std::endl;
  depth = 0;
  std::cout << index.findLeaf(fp10111, 0, depth)->fingerprints[0] << std::endl;
  depth = 0;
  std::cout << index.findLeaf(fp01111, 0, depth)->fingerprints[0] << std::endl;
  
  depth = 0;
  Fingerprint fp11111(5, 31);
  std::cout << index.findLeaf(fp11111, 0, depth) << std::endl;

  depth = 0;
  Fingerprint fp10000(5, 16);
  std::cout << index.findLeaf(fp10000, 0, depth)->fingerprints[0] << std::endl;

  depth = 0;
  Fingerprint fp00011(5, 3);
  std::cout << index.findLeaf(fp00011, 0, depth)->fingerprints[0] << std::endl;

  */

  std::vector<std::size_t> hits;

  std::cout << "######################################################" << std::endl;

  depth = 0;
  //LeafNode *leaf = 0;
//  while ((leaf = index.findLeaf(fp00000, leaf, depth))) {
//    hits.push_back(leaf->fingerprints[0]);  
//  }

  std::cout << "hits: ";
  for (int i = 0; i < hits.size(); ++i)
    std::cout << hits[i] << " ";
  std::cout << std::endl;

  int N = 5;
  
  //std::vector<int> leafs, counts;
  //for (int i = 0; i < std::pow(2, N); ++i)
  //  leafs.push_back(i);
  
  //InMemoryFingerprintStorage storage3("test_database2", "OpenBabel::FP2", N);
  //CompressedBinarySearchTreeFingerprintIndex<InMemoryFingerprintStorage> index3("test_database2", "OpenBabel::FP2", N, storage);
  
  /*
  for (int i = 0; i < std::pow(2, N); ++i) {
    Fingerprint fp(N, i);
    index3.createLeaf(fp).push_back(i);  
  }*/

  /*
  for (int i = 0; i < std::pow(2, N); ++i) {
    depth = 0;
    leaf = 0;
    int count = 0;
    Fingerprint fp(N, i);

    std::vector<std::size_t> hits;
    index3.contains(fp, hits);
    std::cout << "hits " << i << ": ";
    for (int j = 0; j < hits.size(); ++j)
      std::cout << hits[j] << " ";
    std::cout << std::endl;
    //while ((leaf = index3.findLeaf(fp, leaf, depth)))
    //  count++;
    //counts.push_back(count);
  }
  */
  


  /*
  for (int j = 0; j < 1000; ++j) {

    std::random_shuffle(leafs.begin(), leafs.end());
  
    InMemoryFingerprintStorage storage2("test_database2", "OpenBabel::FP2", N);
    CompressedBinarySearchTreeFingerprintIndex<InMemoryFingerprintStorage> index2("test_database2", "OpenBabel::FP2", N, storage);

    for (int i = 0; i < std::pow(2, N); ++i) {
      Fingerprint fp(N, leafs[i]);
      index2.createLeaf(fp).push_back(leafs[i]);  
    }

    for (int i = 0; i < std::pow(2, N); ++i) {
      depth = 0;
      leaf = 0;
      Fingerprint fp(N, i);
      //std::cout << "leafs " << i << ": ";
      int count = 0;
      while ((leaf = index2.findLeaf(fp, leaf, depth)))
        count++;

      MOLDB_COMPARE(counts[i], count);
        //std::cout << leaf->fingerprints[0] << " ";
      //std::cout << std::endl;
    }

  }
  */
 


  /*
  for (int i = 0; i < 64; ++i)
    std::cout << "numLeafNodes(" << i + 1 << ") = " << numLeafNodes(i + 1) << std::endl;
  for (int i = 0; i < 64; ++i)
    std::cout << "numLeafNodes(" << i + 1 << ") = " << numLeafNodes(i + 1) << std::endl;

  std::cout << "divide 1024" << std::endl;
  for (int i = 0; i < 1024; ++i)
    if (1024 % (i + 1) == 0)
      std::cout << i + 1 << std::endl;

  Impl::SingleBitBinaryTree tree(4);
  MOLDB_COMPARE(tree.leafNodes().size(), 16);

  std::vector<BitSet> fingerprints;
  for (int i = 0; i < 16; ++i)
    fingerprints.push_back(BitSet(4, i));

  for (int i = 0; i < 16; ++i)
    tree.insert(fingerprints[i], i);

  for (int i = 0; i < 16; ++i)
    MOLDB_COMPARE(tree.findLeaf(fingerprints[i])[0], i);
*/

}
