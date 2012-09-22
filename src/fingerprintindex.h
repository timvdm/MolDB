#ifndef MOLDB_FINGERPRINTINDEX_H
#define MOLDB_FINGERPRINTINDEX_H

#include "fingerprintstorage.h"
#include "timer.h"

#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <boost/unordered_map.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_CUDA
#include "cuda.h"
#endif

namespace MolDB {

  class FingerprintIndex
  {
    public:
      virtual ~FingerprintIndex()
      {
      }

      virtual std::string name() const = 0;

      virtual std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const = 0;
      virtual std::vector<std::size_t> contains(MessageHandler *msg, const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        return contains(fingerprint, maxResults);
      }

      virtual std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const = 0;

      virtual void reserve(std::size_t size)
      {
      }

      virtual void clear()
      {
      }

      virtual bool create(MessageHandler *msg)
      {
        return true;
      }

      virtual bool load(MessageHandler *msg)
      {
        return true;
      }

      virtual bool startInsert(MessageHandler *msg)
      {
        return true;
      }

      virtual bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        return true;
      }

      virtual bool commit()
      {
        return true;
      }

      virtual bool stopInsert(MessageHandler *msg)
      {
        return true;
      }

      virtual std::size_t memoryUsage() const = 0;

  };


  template<typename FingerprintStorageType>
  class BruteForceFingerprintIndex : public FingerprintIndex
  {
    public:
      BruteForceFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_fingerprints(fingerprints)
      {
      }
      
      std::string name() const
      {
        return "Brute Force Fingerprint Index";
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        //std::cout << "Hits: " << fingerprint << std::endl;

        std::vector<std::size_t> hits;

#ifdef _OPENMP
        #pragma omp parallel
        {
          std::vector<std::size_t> priv_hits;
          #pragma omp for
          for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
            if (isSubsetSuperset(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks()))
              priv_hits.push_back(i);
          }

          #pragma omp critical
          hits.insert(hits.end(), priv_hits.begin(), priv_hits.end());
        }

        if (hits.size() > maxResults)
          hits.resize(maxResults);
#else
        hits.reserve(m_fingerprints.size());
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          if (isSubsetSuperset(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks()))
            hits.push_back(i);
          //if (hits.size() == maxResults)
          //  return hits;
        }
#endif
        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);
        const FingerprintWord *pointer = &query[0];

        std::vector<std::pair<std::size_t, double> > hits;

#ifdef _OPENMP
        #pragma omp parallel
        {
          std::vector<std::pair<std::size_t, double> > priv_hits;
          #pragma omp for
          for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
            double S = MolDB::tanimoto(pointer, m_fingerprints.pointer(i), fingerprint.num_blocks());
            if (S >= threshold)
              priv_hits.push_back(std::make_pair(i, S));
          }

          #pragma omp critical
          hits.insert(hits.end(), priv_hits.begin(), priv_hits.end());
        }

        if (hits.size() > maxResults)
          hits.resize(maxResults);
#else
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          double S = MolDB::tanimoto(pointer, m_fingerprints.pointer(i), fingerprint.num_blocks());
          if (S >= threshold)
            hits.push_back(std::make_pair(i, S));
          if (hits.size() == maxResults)
            return hits;
        } 
#endif
        return hits;
      }

      std::size_t memoryUsage() const
      {
        return sizeof(BruteForceFingerprintIndex<FingerprintStorageType>);
      }

    private:
      const FingerprintStorageType &m_fingerprints;
  };

  /**
   * @class InvertedFingerprintIndex fingerprintindex.h
   *
   * Similar to 1D kD-Grid with smaller memory usage.
   */
  template<typename FingerprintStorageType>
  class InvertedFingerprintIndex : public FingerprintIndex
  {
      void insert(const Fingerprint &fingerprint, std::size_t index)
      {
        std::vector<int> bits;
        
        Fingerprint::size_type pos = fingerprint.find_first();
        while (pos != Fingerprint::npos) {
          bits.push_back(pos);
          pos = fingerprint.find_next(pos);
        }
        
        for (int i = 0; i < bits.size(); ++i)
          m_setBits[bits[i]].push_back(index);
      }

    public:
      InvertedFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_fingerprints(fingerprints), m_setBits(numBits)
      {
      }

      std::string name() const
      {
        return "Random Access Inverted Fingerprint Index";
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<int> bits;
        std::vector<std::size_t> hits;

        Fingerprint::size_type pos = fingerprint.find_first();
        while (pos != Fingerprint::npos) {
          bits.push_back(pos);
          pos = fingerprint.find_next(pos);
        }

        if (bits.empty()) {
          for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
            hits.push_back(i);
            if (hits.size() == maxResults)
              return hits;
          }

          return hits;
        }
       
        // select the shortest leaf
        int bit = 0;
        for (int i = 0; i < bits.size(); ++i) {
          const std::vector<std::size_t> &bitLeaf = m_setBits[bits[bit]];
          const std::vector<std::size_t> &leaf = m_setBits[bits[i]];
          if (leaf.size() < bitLeaf.size())
            bit = i;
        }

        const std::vector<std::size_t> &bitLeaf = m_setBits[bits[bit]];
        for (std::size_t i = 0; i < bitLeaf.size(); ++i) {
          bool found = true;
          for (int j = 0; j < bits.size(); ++j) {
            if (j == bit)
              continue;
            
            const std::vector<std::size_t> &leaf = m_setBits[bits[j]];
            if (!std::binary_search(leaf.begin(), leaf.end(), bitLeaf[i])) {
              found = false;
              break;
            }
          }

          if (found)
            hits.push_back(bitLeaf[i]);
        }

        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::pair<std::size_t, double> > hits;
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          double S = MolDB::tanimoto(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks());
          if (S >= threshold)
            hits.push_back(std::make_pair(i, S));
          if (hits.size() == maxResults)
            return hits;
        }
 
        return hits;
      }

      void clear()
      {
        m_setBits.clear();
        m_setBits.resize(m_fingerprints.numBits());
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        insert(fingerprint, m_fingerprints.size() - 1);
        return true;
      }

      bool load(MessageHandler *msg)
      {
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          insert(m_fingerprints[i], i);
        }


        for (std::size_t i = 0; i < m_setBits.size(); ++i) {
          std::ofstream ofs(make_string("inverted_index_distribution", i + 1, ".csv").c_str());
          for (std::size_t j = 0; j < m_setBits[i].size(); ++j)
            ofs << m_setBits[i][j] << std::endl;
          ofs << std::endl;
        }

        return true;
      }

      std::size_t memoryUsage() const
      {
        std::size_t size = 0;
        for (std::size_t i = 0; i < m_setBits.size(); ++i)
          size += sizeof(std::vector<std::size_t>) + m_setBits[i].size() * sizeof(std::size_t);
        return sizeof(InvertedFingerprintIndex<FingerprintStorageType>) + size;
      }

    private:
      const FingerprintStorageType &m_fingerprints;
      std::vector<std::vector<std::size_t> > m_setBits;
  };

  template<typename FingerprintStorageType>
  class SequentialAccessInvertedFingerprintIndex : public FingerprintIndex
  {
      void insert(const Fingerprint &fingerprint, std::size_t index)
      {
        std::vector<int> bits;
        
        Fingerprint::size_type pos = fingerprint.find_first();
        while (pos != Fingerprint::npos) {
          bits.push_back(pos);
          pos = fingerprint.find_next(pos);
        }
        
        for (int i = 0; i < bits.size(); ++i)
          m_setBits[bits[i]].push_back(index);
      }

    public:
      SequentialAccessInvertedFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_fingerprints(fingerprints), m_setBits(numBits)
      {
        std::string sequentialFilename = filename + "." + convertFingerprintName(fingerprint) + ".saifi.moldb";
        m_sequentialFingerprints = FingerprintStorageType(sequentialFilename, fingerprint, numBits);
      }

      std::string name() const
      {
        return "Sequential Access Inverted Fingerprint Index";
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<int> bits;
        std::vector<std::size_t> hits;

        Fingerprint::size_type pos = fingerprint.find_first();
        while (pos != Fingerprint::npos) {
          bits.push_back(pos);
          pos = fingerprint.find_next(pos);
        }

        if (bits.empty()) {
          for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
            hits.push_back(i);
            if (hits.size() == maxResults)
              return hits;
          }

          return hits;
        }
       
        // select the shortest leaf
        int bit = 0;
        for (int i = 0; i < bits.size(); ++i) {
          const std::vector<std::size_t> &bitLeaf = m_setBits[bits[bit]];
          const std::vector<std::size_t> &leaf = m_setBits[bits[i]];
          if (leaf.size() < bitLeaf.size())
            bit = i;
        }

        const std::vector<std::size_t> &bitLeaf = m_setBits[bits[bit]];
        for (std::size_t i = 0; i < bitLeaf.size(); ++i) {
          bool found = true;
          for (int j = 0; j < bits.size(); ++j) {
            if (j == bit)
              continue;
            
            const std::vector<std::size_t> &leaf = m_setBits[bits[j]];
            if (!std::binary_search(leaf.begin(), leaf.end(), bitLeaf[i])) {
              found = false;
              break;
            }
          }

          if (found)
            hits.push_back(bitLeaf[i]);
        }

        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::pair<std::size_t, double> > hits;
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          double S = MolDB::tanimoto(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks());
          if (S >= threshold)
            hits.push_back(std::make_pair(i, S));
          if (hits.size() == maxResults)
            return hits;
        }
 
        return hits;
      }

      void clear()
      {
        m_setBits.clear();
        m_setBits.resize(m_fingerprints.numBits());
      }

      bool stopInsert(MessageHandler *msg)
      {
        load(msg);
        return true;
      }

      bool load(MessageHandler *msg)
      {
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          insert(m_fingerprints[i], i);
          
        m_sequentialMap.resize(m_setBits.size());
        m_sequentialFingerprints.create(msg);
        m_sequentialFingerprints.startInsert(msg);
        m_sequentialFingerprints.reserve(m_fingerprints.size());
 
        std::size_t sequentialIndex = 0;
        for (std::size_t i = 0; i < m_setBits.size(); ++i)
          for (std::size_t j = 0; j < m_setBits[i].size(); ++i) {
            m_sequentialFingerprints.insert(msg, m_fingerprints[m_setBits[i][j]]);
            m_sequentialMap[m_setBits[i][j]] = sequentialIndex;
            ++sequentialIndex;
          }
        
        m_sequentialFingerprints.stopInsert(msg);
        
        return true;
      }

      std::size_t memoryUsage() const
      {
        std::size_t size = m_sequentialFingerprints.memoryUsage() + m_sequentialMap.size() * sizeof(std::size_t);
        for (std::size_t i = 0; i < m_setBits.size(); ++i)
          size += sizeof(std::vector<std::size_t>) + m_setBits[i].size() * sizeof(std::size_t);
        return sizeof(SequentialAccessInvertedFingerprintIndex<FingerprintStorageType>) + size;
      }

    private:
      const FingerprintStorageType &m_fingerprints;
      FingerprintStorageType m_sequentialFingerprints;
      std::vector<std::vector<std::size_t> > m_setBits;
      /**
       * random -> sequential index
       */
      std::vector<std::size_t> m_sequentialMap;
  };
  
  namespace Impl {

    template<typename FingerprintStorageType>
    class FingerprintBinarySearchTree
    {
      public:
        FingerprintBinarySearchTree(unsigned int height) : m_height(height), m_numNodes(std::pow(2, height + 1) - 1),
            m_numLeafNodes(std::pow(2, height))
        {
        }

        unsigned long parentIndex(unsigned long index) const
        {
          return (index - 1) / 2;
        }

        unsigned long leftIndex(unsigned long index) const
        {
          return 2 * index + 1;
        }

        unsigned long rightIndex(unsigned long index) const
        {
          return 2 * index + 2;
        }

        bool isLeft(unsigned long index) const
        {
          return (index == 0) ? false : index % 2 == 1;
        }

        bool isRight(unsigned long index) const
        {
          return (index == 0) ? false : index % 2 == 0;
        }

        int height() const
        {
          return m_height;
        }

        std::size_t numNodes() const
        {
          return m_numNodes;
        }

        std::size_t numLeafNodes() const
        {
          return m_numLeafNodes;
        }

        const std::vector<std::size_t>& leafNode(unsigned long index) const
        {
          return m_leafs[index - m_numLeafNodes + 1];
        }

        std::vector<std::size_t>& leafNode(unsigned long index)
        {
          return m_leafs[index - m_numLeafNodes + 1];
        }

        const std::vector<std::vector<std::size_t> >& leafNodes() const
        {
          return m_leafs;
        }

        std::vector<std::vector<std::size_t> >& leafNodes()
        {
          return m_leafs;
        }

        template<typename GoRightFunctor>
        unsigned long findLeafNode(const Fingerprint &fingerprint, const GoRightFunctor &functor) const
        {
          int depth = 0;
          unsigned long index = 0;

          while (depth < m_height) {
            if (functor(fingerprint, depth, index))
              index = rightIndex(index);
            else
              index = leftIndex(index);
            ++depth;
          }

          return index;
        }

        template<typename GoRightFunctor>
        unsigned long findNextLeafNode(const Fingerprint &fingerprint, unsigned long lastLeafIndex, const GoRightFunctor &functor) const
        {
          if (lastLeafIndex == numNodes() - 1)
            return 0;

          int depth = m_height;
          unsigned long index = lastLeafIndex;

          while (isRight(index)) {
            index = parentIndex(index);
            --depth;
          }

          if (index)
            index = rightIndex(parentIndex(index));

          while (depth < m_height) {
            if (functor(fingerprint, depth, index))
              index = rightIndex(index);
            else
              index = leftIndex(index);
            ++depth;
          }

          return index;
        }

        void clear()
        {
          for (std::size_t i = 0; i < m_leafs.size(); ++i)
            m_leafs[i].clear();
        }

        std::size_t memoryUsage() const
        {
          std::size_t leafMemoryUsage = m_leafs.size() * sizeof(std::vector<std::size_t>);
          for (std::size_t i = 0; i < m_leafs.size(); ++i)
            leafMemoryUsage += m_leafs[i].size() * sizeof(std::size_t);

          return sizeof(FingerprintBinarySearchTree<FingerprintStorageType>) + leafMemoryUsage;
        }

      private:
        unsigned int m_height;
        unsigned long m_numNodes;
        unsigned long m_numLeafNodes;
        std::vector<std::vector<std::size_t> > m_leafs;
    };

  }

  template<typename FingerprintStorageType, int Height, int BitsPerNode> 
  class BitSelectionBinarySearchTreeFingerprintIndex : public FingerprintIndex
  {
      typedef Impl::FingerprintBinarySearchTree<FingerprintStorageType> BinarySearchTree;

      struct GoRightFunctor
      {
        GoRightFunctor(const std::vector<unsigned short> &bits_) : bits(bits_)
        {
        }

        bool operator()(const Fingerprint &fingerprint, int depth, unsigned long index) const
        {
          for (int i = 0; i < BitsPerNode; ++i)
            if (!fingerprint[bits[BitsPerNode * index + i]])
              return false;
          return true;
        }

        const std::vector<unsigned short> &bits;
      };

      struct GoRightFunctor2
      {
        GoRightFunctor2(const Fingerprint &bits_) : bits(bits_)
        {
        }

        bool operator()(const Fingerprint &fingerprint, int depth, unsigned long index) const
        {
          return bits[depth];
        }

        const Fingerprint &bits;
      };

      void findBits(MessageHandler *msg)
      {
        msg->report(Information, "Optimizing fingerprint index...");
        msg->report(Information, progressBar(0, Height));

        Timer timer;

        int depth = 0;
        std::vector<std::vector<std::size_t> > bins(1);
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          bins.back().push_back(i);

        while (true) {
          if (depth == Height)
            break;

          std::vector<std::vector<std::size_t> > nextBins;

          for (std::size_t i = 0; i < bins.size(); ++i) {
            const std::vector<std::size_t> &bin = bins[i];

            std::vector<int> bits;
            for (int l = 0; l < BitsPerNode; ++l) {

              std::vector<std::size_t> diffs;
              for (int j = 0; j < m_fingerprints.numBits(); ++j) {
                if (std::find(bits.begin(), bits.end(), j) != bits.end()) {
                  diffs.push_back(std::numeric_limits<std::size_t>::max());
                  continue;
                }

                std::size_t set0 = 0;
                std::size_t set1 = 0;
                for (std::size_t k = 0; k < bin.size(); ++k) {
                  bool inSet1 = true;
                  for (int m = 0; m < bits.size(); ++m) {
                    if (!m_fingerprints[bin[k]][bits[m]]) {
                      inSet1 = false;
                      break;
                    }
                  }

                  if (inSet1 && m_fingerprints[bin[k]][j])
                    ++set1;
                  else
                    ++set0;
                }
                std::size_t diff = std::labs(set0 - set1);
                diffs.push_back(diff);
                if (!diff)
                  break;
              }

              unsigned short bit = std::min_element(diffs.begin(), diffs.end()) - diffs.begin();
              m_bits.push_back(bit);
              bits.push_back(bit);
            }

            nextBins.resize(nextBins.size() + 2);
            for (std::size_t j = 0; j < bin.size(); ++j) {
              bool right = true;
              for (int k = 0; k < BitsPerNode; ++k) {
                if (!m_fingerprints[bin[j]][bits[k]]) {
                  right = false;
                  break;
                }
              }

              if (right)
                nextBins[nextBins.size() - 1].push_back(bin[j]);
              else
                nextBins[nextBins.size() - 2].push_back(bin[j]);
            }
          }

          std::swap(bins, nextBins);
          ++depth;

          msg->report(Information, progressBar(depth, Height));
        }

        std::swap(bins, m_tree.leafNodes());
        
        msg->report(Information, make_string("Index build in ", timer.elapsed(), " seconds"));
      }

    public:
      BitSelectionBinarySearchTreeFingerprintIndex(const std::string &filename, const std::string &fingerprint,
          unsigned int numBits, const FingerprintStorageType &fingerprints) : m_fingerprints(fingerprints), m_tree(Height)
      {
        // the filename for the node bits
        std::string fp = fingerprint;
        fp.replace(fingerprint.find("::"), 2, "_");
        m_filename = make_string(filename, ".", fp, ".ft_", m_tree.height(), "_", BitsPerNode, ".moldb");
      }

      std::string name() const
      {
        return make_string("Bit Selection Binary Search Tree Fingerprint Index (Height = ", m_tree.height(), ")");
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::size_t> hits;
        for (unsigned long index = m_tree.findLeafNode(fingerprint, GoRightFunctor(m_bits)); index; index = m_tree.findNextLeafNode(fingerprint, index, GoRightFunctor(m_bits))) {
          const std::vector<std::size_t> &leaf = m_tree.leafNode(index);
          for (std::size_t i = 0; i < leaf.size(); ++i) {
            if (isSubsetSuperset(&query[0], m_fingerprints.pointer(leaf[i]), fingerprint.num_blocks()))
              hits.push_back(leaf[i]);
            if (hits.size() == maxResults)
              return hits;
          }
        }

        /*
        unsigned long leafIndex = m_tree.findLeafNode(fingerprint, GoRightFunctor(m_bits));
        std::size_t index = leafIndex - m_tree.numLeafNodes();

        for (std::size_t i = 0; i < m_nextLeafs[index].size(); ++i) {
          const std::vector<std::size_t> &leaf = m_tree.leafNode(m_nextLeafs[index][i]);
          for (std::size_t j = 0; j < leaf.size(); ++j) {
            if (isSubsetSuperset(&query[0], m_fingerprints.pointer(leaf[j]), fingerprint.num_blocks()))
              hits.push_back(leaf[j]);
            if (hits.size() == maxResults)
              return hits;
          }
        }
        */

        
        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::pair<std::size_t, double> > hits;
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          double S = MolDB::tanimoto(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks());
          if (S >= threshold)
            hits.push_back(std::make_pair(i, S));
          if (hits.size() == maxResults)
            return hits;
        }

        return hits;
      }

      void clear()
      {
        m_bits.clear();
        m_tree.clear();
      }

      bool load(MessageHandler *msg)
      {
        // search for file with greater tree height if needed
        std::string filename = m_filename;
        if (!std::ifstream(m_filename.c_str(), std::ios_base::in | std::ios_base::binary)) {
          std::string extension = make_string(".ft_", m_tree.height(), "_", BitsPerNode, ".moldb");
          std::size_t pos = filename.rfind(extension);
          for (int i = Height; i < 50; ++i) {
            filename.replace(pos, std::string::npos, make_string(".ft_", i, "_", BitsPerNode, ".moldb"));
            if (std::ifstream(filename.c_str(), std::ios_base::in | std::ios_base::binary)) {
              msg->report(Information, make_string("Fingerprint binary search tree file ", m_filename, " not found, but ", filename, " can be used."));
              break;
            }
          }
        }

        // load the bits for the search tree
        std::fstream ifs(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        if (!ifs) {
          msg->report(Warning, make_string("Fingerprint binary search tree file ", m_filename, " not found, rebuilding index."));
          // build the binary search tree index
          stopInsert(msg);
          return true;
        }
        ifs.close();

        if (!IO::openFile(msg, ifs, filename, std::ios_base::in, IO::MagicNumber::FingerprintNodeBitsFileFormat))
          return false;

        // get the number of bits in the fingerprint
        unsigned long numBits = IO::read<unsigned long>(ifs);
        if (m_fingerprints.size() && numBits < BitsPerNode * std::pow(2, m_tree.height()) - BitsPerNode) {
          msg->report(Error, make_string("Fingerprint binary search file ", filename, " does not contain the correct number of bits."));
          return false;
        }

        // start loading
        m_bits.reserve(numBits);
        for (std::size_t i = 0; i < numBits; ++i) {
          unsigned short bit = IO::read<unsigned short>(ifs);
          m_bits.push_back(bit);
        }

        // load the tree
        m_tree.leafNodes().resize(m_tree.numLeafNodes());
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          m_tree.leafNodes()[m_tree.findLeafNode(m_fingerprints[i], GoRightFunctor(m_bits)) - m_tree.numLeafNodes() + 1].push_back(i);
       
        /*
        m_nextLeafs.resize(m_tree.numLeafNodes());
        for (std::size_t i = 0; i < m_tree.numLeafNodes(); ++i) {
          Fingerprint fp(m_tree.height(), i);
          for (unsigned long index = m_tree.findLeafNode(fp, GoRightFunctor2(fp)); index; index = m_tree.findNextLeafNode(fp, index, GoRightFunctor2(fp)))
            m_nextLeafs[i].push_back(index);
        }
        */
 
        return true;
      }

      bool stopInsert(MessageHandler *msg)
      {
        findBits(msg);

        // write the bits for the search tree
        std::fstream ofs;
        if (!IO::openFile(msg, ofs, m_filename, std::ios_base::out, IO::MagicNumber::FingerprintNodeBitsFileFormat))
          return false;

        // write the number of bits in the search tree
        unsigned long numBits = m_bits.size();
        IO::write(ofs, numBits);

        // start writing
        for (std::size_t i = 0; i < numBits; ++i)
          IO::write(ofs, m_bits[i]);

        return true;
      }

      std::size_t memoryUsage() const
      {
        return sizeof(BitSelectionBinarySearchTreeFingerprintIndex<FingerprintStorageType, Height, BitsPerNode>) + 
               m_tree.memoryUsage() + (m_tree.numLeafNodes() - 1) * sizeof(unsigned short);
      }

    private:
      std::string m_filename;
      const FingerprintStorageType &m_fingerprints;
      BinarySearchTree m_tree;
      std::vector<unsigned short> m_bits;
      std::vector<std::vector<std::size_t> > m_nextLeafs;
  };

  template<typename FingerprintStorageType, int Dimension>
  class KDBinarySearchTreeFingerprintIndex : public FingerprintIndex
  {
      typedef Impl::FingerprintBinarySearchTree<FingerprintStorageType> BinarySearchTree;

      struct GoRightFunctor
      {
        GoRightFunctor(const std::vector<unsigned short> &bitCounts_) : bitCounts(bitCounts_)
        {
        }

        bool operator()(const Fingerprint &fingerprint, int depth, unsigned long index) const
        {
          return bitCount(fingerprint, depth) > bitCounts[index];
        }

        const std::vector<unsigned short> &bitCounts;
      };

      static unsigned short bitCount(const Fingerprint &fingerprint, int depth)
      {
        int beginBit = depth * (fingerprint.size() / Dimension);
        int endBit = (depth + 1) * (fingerprint.size() / Dimension);
 
        int count = 0;
        for (int i = beginBit; i < endBit; ++i)
          if (fingerprint[i])
            ++count;

        return count;
      }

      void findBitCounts(MessageHandler *msg)
      {
        msg->report(Information, "Optimizing fingerprint index...");
        msg->report(Information, progressBar(0, Dimension));

        Timer timer;

        int depth = 0;
        std::vector<std::vector<std::size_t> > bins(1);
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          bins.back().push_back(i);

        // precompute the bit counts
        std::vector<std::vector<unsigned short> > bitCounts(Dimension);
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          for (unsigned short j = 0; j < Dimension; ++j)
            bitCounts[j].push_back(bitCount(m_fingerprints[i], j));

        while (true) {
          if (depth == Dimension)
            break;

          std::vector<std::vector<std::size_t> > nextBins;

          for (std::size_t i = 0; i < bins.size(); ++i) {
            const std::vector<std::size_t> &bin = bins[i];

            std::vector<std::size_t> diffs;
            for (int j = 0; j < m_fingerprints.numBits(); ++j) {
              std::size_t set0 = 0;
              std::size_t set1 = 0;
              for (std::size_t k = 0; k < bin.size(); ++k) {
                if (bitCounts[depth][bin[k]] > j)
                  ++set1;
                else
                  ++set0;
              }
              std::size_t diff = std::labs(set0 - set1);
              diffs.push_back(diff);
              if (!diff)
                break;
            }

            unsigned short count = std::min_element(diffs.begin(), diffs.end()) - diffs.begin();
            m_bitCounts.push_back(count);

            nextBins.resize(nextBins.size() + 2);
            for (std::size_t j = 0; j < bin.size(); ++j) {
              if (bitCounts[depth][bin[j]] > m_bitCounts.back())
                nextBins[nextBins.size() - 1].push_back(bin[j]);
              else
                nextBins[nextBins.size() - 2].push_back(bin[j]);
            }

          }

          std::swap(bins, nextBins);
          ++depth;

          msg->report(Information, progressBar(depth, Dimension));
        }

        std::swap(bins, m_tree.leafNodes());
        
        msg->report(Information, make_string("Index build in ", timer.elapsed(), " seconds"));
      }

    public:
      KDBinarySearchTreeFingerprintIndex(const std::string &filename, const std::string &fingerprint,
          unsigned int numBits, const FingerprintStorageType &fingerprints) : m_fingerprints(fingerprints),
          m_tree(Dimension)
      {
        // the filename for the node bits
        std::string fp = fingerprint;
        fp.replace(fingerprint.find("::"), 2, "_");
        m_filename = make_string(filename, ".", fp, ".fkt_", m_tree.height(), ".moldb");
      }

      std::string name() const
      {
        return make_string("kD Binary Search Tree Fingerprint Index (k = ", m_tree.height(), ")");
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::size_t> hits;
        for (unsigned long index = m_tree.findLeafNode(fingerprint, GoRightFunctor(m_bitCounts)); index; index = m_tree.findNextLeafNode(fingerprint, index, GoRightFunctor(m_bitCounts))) {
          const std::vector<std::size_t> &leaf = m_tree.leafNode(index);
          for (std::size_t i = 0; i < leaf.size(); ++i) {
            if (isSubsetSuperset(&query[0], m_fingerprints.pointer(leaf[i]), fingerprint.num_blocks()))
              hits.push_back(leaf[i]);
            if (hits.size() == maxResults)
              return hits;
          }
        }
        
        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::pair<std::size_t, double> > hits;
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          double S = MolDB::tanimoto(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks());
          if (S >= threshold)
            hits.push_back(std::make_pair(i, S));
          if (hits.size() == maxResults)
            return hits;
        }

        return hits;
      }

      void clear()
      {
        m_tree.clear();
      }

      bool load(MessageHandler *msg)
      {
        // search for file with greater tree height if needed
        std::string filename = m_filename;
        if (!std::ifstream(m_filename.c_str(), std::ios_base::in | std::ios_base::binary)) {
          std::string extension = make_string(".fkt_", m_tree.height(), ".moldb");
          std::size_t pos = filename.rfind(extension);
          for (int i = Dimension; i < 50; ++i) {
            filename.replace(pos, std::string::npos, make_string(".fkt_", i, ".moldb"));
            if (std::ifstream(filename.c_str(), std::ios_base::in | std::ios_base::binary)) {
              msg->report(Information, make_string("Fingerprint k-d binary search tree file ", m_filename, " not found, but ", filename, " can be used."));
              break;
            }
          }
        }

        // load the bit counts for the search tree
        std::fstream ifs(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        if (!ifs) {
          msg->report(Warning, make_string("Fingerprint k-d binary search tree file ", m_filename, " not found, rebuilding index."));
          // build the binary search tree index
          stopInsert(msg);
          return true;
        }
        ifs.close();

        if (!IO::openFile(msg, ifs, filename, std::ios_base::in, IO::MagicNumber::FingerprintNodeBitsFileFormat))
          return false;

        // get the number of bits in the fingerprint
        unsigned long numBitCounts = IO::read<unsigned long>(ifs);
        if (m_fingerprints.size() && numBitCounts < std::pow(2, m_tree.height()) - 1) {
          msg->report(Error, make_string("Fingerprint k-d binary search file ", filename, " does not contain the correct number of bit counts."));
          return false;
        }

        // start loading
        m_bitCounts.reserve(numBitCounts);
        for (std::size_t i = 0; i < numBitCounts; ++i) {
          unsigned short bit = IO::read<unsigned short>(ifs);
          m_bitCounts.push_back(bit);
        }

        // load the tree
        m_tree.leafNodes().resize(m_tree.numLeafNodes());
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          m_tree.leafNodes()[m_tree.findLeafNode(m_fingerprints[i], GoRightFunctor(m_bitCounts)) - m_tree.numLeafNodes() + 1].push_back(i);
 
        return true;
      }

      bool stopInsert(MessageHandler *msg)
      {
        findBitCounts(msg);

        // write the bit counts for the search tree
        std::fstream ofs;
        if (!IO::openFile(msg, ofs, m_filename, std::ios_base::out, IO::MagicNumber::FingerprintNodeBitsFileFormat))
          return false;

        // write the number of bit counts in the search tree
        unsigned long numBitCounts = m_bitCounts.size();
        IO::write(ofs, numBitCounts);

        // start writing
        for (std::size_t i = 0; i < numBitCounts; ++i)
          IO::write(ofs, m_bitCounts[i]);

        return true;
      }

      std::size_t memoryUsage() const
      {
        return sizeof(KDBinarySearchTreeFingerprintIndex<FingerprintStorageType, Dimension>) + m_tree.memoryUsage();
      }

    private:
      std::string m_filename;
      const FingerprintStorageType &m_fingerprints;
      BinarySearchTree m_tree;
      std::vector<unsigned short> m_bitCounts;
  };

  template<typename FingerprintStorageType, bool InMemoryTree = true>
  class CompressedBinarySearchTreeFingerprintIndex : public FingerprintIndex
  {
      struct Node
      {
        virtual ~Node() {}
      };

      struct TreeNode : public Node
      {
        TreeNode() : leftDistance(0), rightDistance(0), left(0), right(0)
        {
        }

        int leftDistance;
        int rightDistance;
        Node *left;
        Node *right;
      };

      struct LeafNode : public Node
      {
        std::vector<std::size_t> fingerprints;
      };

      template<typename NodeType>
      NodeType *allocateNode()
      {
        if (TypeTraits::IsSame<NodeType, TreeNode>::result)
          ++m_nodeCount;
        else
          ++m_leafCount;

        if (InMemoryTree)
          return new NodeType;

        assert(m_allocationPointer + sizeof(NodeType) - m_mappedFile.data() <= 5000000000);
        NodeType *node = new (m_allocationPointer) NodeType;
        m_allocationPointer += sizeof(NodeType);
        return node;
      }

      std::vector<std::size_t>& createLeaf(const Fingerprint &fingerprint)
      {
        int depth = 0;
        TreeNode *node = m_tree;

        while (depth < m_numBits) {
          bool right = fingerprint[depth];
            
          int distance = -1;
          while (fingerprint[depth] == right) {
            ++depth;
            ++distance;
            if (depth == m_numBits)
              break;
          }
          
          while (true) {
            Node *child = right ? node->right : node->left;
            int nodeDistance = right ? node->rightDistance : node->leftDistance;
            
            if (!nodeDistance) {
              // there is a right/left child at distance 0
              if (child) {
                if (distance) {
                  --distance;
                  node = static_cast<TreeNode*>(child);
                  continue;
                }
                // done
                if (depth == m_numBits)
                  return static_cast<LeafNode*>(child)->fingerprints;

                node = static_cast<TreeNode*>(child);
                break;
              }

              // insert a new child at distance distance
              child = (depth == m_numBits) ? static_cast<Node*>(allocateNode<LeafNode>()) : static_cast<Node*>(allocateNode<TreeNode>());
              if (right) {
                node->rightDistance = distance;
                node->right = child;
              } else {
                node->leftDistance = distance;
                node->left = child;
              }
              // done
              if (depth == m_numBits)
                return static_cast<LeafNode*>(child)->fingerprints;

              node = static_cast<TreeNode*>(child);
              break;
            }

            if (distance > nodeDistance) {
              // continue to child and decrement distance
              distance -= nodeDistance + 1;
              node = static_cast<TreeNode*>(child);
              continue;
            }
              
            if (distance == nodeDistance) {
              // continue to the child at same distance
              if (depth == m_numBits)
                return static_cast<LeafNode*>(child)->fingerprints;

              node = static_cast<TreeNode*>(child);
              break;
            }
              
            if (distance < nodeDistance) {
              // insert new node between node and child
              TreeNode *newNode = allocateNode<TreeNode>();
              
              if (right) {
                node->rightDistance = distance;
                node->right = newNode;
                newNode->rightDistance = nodeDistance - distance - 1;
                newNode->right = child;
              } else {
                node->leftDistance = distance;
                node->left = newNode;
                newNode->leftDistance = nodeDistance - distance - 1;
                newNode->left = child;
              }

              // done
              node = newNode;
              break;
            }
          }
        }

        assert(0);
      }

      void containsDFS(const Fingerprint &fingerprint, TreeNode *node, int depth, std::vector<std::size_t> &hits) const
      {
        bool right = fingerprint[depth];

        if (!right && node->left) {
          int distance = 0;
          if (depth + 1 < m_numBits) {
            while (fingerprint[depth + distance + 1] == right) {
              ++distance;
              if (depth + distance + 1 == m_numBits)
                break;
            }
          }
 
          if (distance >= node->leftDistance) {
            if (depth + node->leftDistance + 1 < m_numBits)
              containsDFS(fingerprint, static_cast<TreeNode*>(node->left), depth + node->leftDistance + 1, hits);
            else {
              LeafNode *leaf = static_cast<LeafNode*>(node->left);
              std::copy(leaf->fingerprints.begin(), leaf->fingerprints.end(), std::back_inserter(hits));
            }
          }
        }

        if (node->right) {
          if (depth + node->rightDistance + 1 < m_numBits)
            containsDFS(fingerprint, static_cast<TreeNode*>(node->right), depth + node->rightDistance + 1, hits);
          else {
            LeafNode *leaf = static_cast<LeafNode*>(node->right);
            std::copy(leaf->fingerprints.begin(), leaf->fingerprints.end(), std::back_inserter(hits));
          }
        }
      }
 
      void tanimotoDFS(const Fingerprint &fingerprint, double threshold, std::size_t maxResults, TreeNode *node,
          int depth, int runningUnionCount, int runningIntersectionCount, int bitCount,
          std::vector<std::pair<std::size_t, double> > &hits) const
      {
        if (node->left) {
          int leftRUC = runningUnionCount;
          int leftBitCount = bitCount;
          for (int i = depth; i < depth + node->leftDistance + 1; ++i) {
            if (fingerprint[i]) {
              ++leftRUC;
              --leftBitCount;
            }
          }

          double S = static_cast<double>(runningIntersectionCount + leftBitCount) / (leftRUC + leftBitCount);

          if (S > threshold) {
            if (depth + node->leftDistance + 1 < m_numBits)
              tanimotoDFS(fingerprint, threshold, maxResults, static_cast<TreeNode*>(node->left), 
                  depth + node->leftDistance + 1, leftRUC, runningIntersectionCount, leftBitCount, hits);
            else {
              LeafNode *leaf = static_cast<LeafNode*>(node->left);
              for (std::size_t i = 0; i < leaf->fingerprints.size(); ++i) {
                hits.push_back(std::make_pair(leaf->fingerprints[i], S));
                if (hits.size() == maxResults)
                  break;
              }
            }
          }
        }

        if (hits.size() == maxResults)
          return;

        if (node->right) {
          int rightRUC = runningUnionCount + node->rightDistance + 1;
          int rightRIC = runningIntersectionCount;
          int rightBitCount = bitCount;
          for (int i = depth; i < depth + node->rightDistance + 1; ++i) {
            if (fingerprint[i]) {
              ++rightRIC;
              --rightBitCount;
            }
          }
          
          double S = static_cast<double>(rightRIC + rightBitCount) / (rightRUC + rightBitCount);

          if (S > threshold) {
            if (depth + node->rightDistance + 1 < m_numBits)
              tanimotoDFS(fingerprint, threshold, maxResults, static_cast<TreeNode*>(node->right),
                  depth + node->rightDistance + 1, rightRUC, rightRIC, rightBitCount, hits);
            else {
              LeafNode *leaf = static_cast<LeafNode*>(node->right);
              for (std::size_t i = 0; i < leaf->fingerprints.size(); ++i) {
                hits.push_back(std::make_pair(leaf->fingerprints[i], S));
                if (hits.size() == maxResults)
                  break;
              }
            }
          }
        }

      }
     
      void clearDFS(TreeNode *node, int depth)
      {
        if (node->left) {
          if (depth + node->leftDistance + 1 < m_numBits)
            clearDFS(static_cast<TreeNode*>(node->left), depth + node->leftDistance + 1);
          else if (InMemoryTree)
            delete static_cast<LeafNode*>(node->left);
          else
            static_cast<LeafNode*>(node->left)->~LeafNode();
        }

        if (node->right) {
          if (depth + node->rightDistance + 1 < m_numBits)
            clearDFS(static_cast<TreeNode*>(node->right), depth + node->rightDistance + 1);
          else if (InMemoryTree)
            delete static_cast<LeafNode*>(node->right);
          else
            static_cast<LeafNode*>(node->right)->~LeafNode();
        }
       
        if (InMemoryTree)
          delete node;
        else
          node->~TreeNode();
      }

      std::size_t memoryUsageDFS(TreeNode *node, int depth) const
      {
        std::size_t size = sizeof(TreeNode);

        if (node->left) {
          if (depth + node->leftDistance + 1 < m_numBits)
            size += memoryUsageDFS(static_cast<TreeNode*>(node->left), depth + node->leftDistance + 1);
          else {
            LeafNode *leaf = static_cast<LeafNode*>(node->left);
            size += sizeof(LeafNode) + leaf->fingerprints.size() * sizeof(std::size_t);
          }
        }

        if (node->right) {
          if (depth + node->rightDistance + 1 < m_numBits)
            size += memoryUsageDFS(static_cast<TreeNode*>(node->right), depth + node->rightDistance + 1);
          else {
            LeafNode *leaf = static_cast<LeafNode*>(node->right);
            size += sizeof(LeafNode) + leaf->fingerprints.size() * sizeof(std::size_t);
          }
        }
 
        return size;        
      }

    public:
      CompressedBinarySearchTreeFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_fingerprints(fingerprints), m_numBits(numBits), m_nodeCount(0), m_leafCount(0)
      {
        if (!InMemoryTree) {
          // the filename for the node bits
          std::string fp = fingerprint;
          fp.replace(fingerprint.find("::"), 2, "_");
          std::string mappedFilename = make_string(filename, ".", fp, ".memory.moldb");
 
          std::ofstream ofs(mappedFilename.c_str());
          ofs.seekp(5000000000);
          ofs << "A";
          ofs.close();

          m_mappedFile.open(mappedFilename);
          m_allocationPointer = m_mappedFile.data();
        }
      
        m_tree = allocateNode<TreeNode>();
      }

      ~CompressedBinarySearchTreeFingerprintIndex()
      {
        if (InMemoryTree)
          clearDFS(m_tree, 0);
      }

      std::string name() const
      {
        return "Compressed Binary Search Tree Fingerprint Index";
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<std::size_t> hits;
        containsDFS(fingerprint, m_tree, 0, hits);

        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<std::pair<std::size_t, double> > hits;
        tanimotoDFS(fingerprint, threshold, maxResults, m_tree, 0, 0, 0, fingerprint.count(), hits);
        
        return hits;
      }

      void clear()
      {
        clearDFS(m_tree, 0);
        m_allocationPointer = m_mappedFile.data();
        m_tree = allocateNode<TreeNode>();
      }

      bool load(MessageHandler *msg)
      {
        // clear the index;
        clear();
       
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          createLeaf(m_fingerprints[i]).push_back(i);

        std::cout << "nodeCount: " << m_nodeCount << std::endl;
        std::cout << "leafCount: " << m_leafCount << std::endl;

        return true;
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        createLeaf(fingerprint).push_back(m_fingerprints.size() - 1);
        return true;
      }

      std::size_t memoryUsage() const
      {
        if (InMemoryTree)
          return sizeof(std::vector<Fingerprint>) + memoryUsageDFS(m_tree, 0);

        return sizeof(std::vector<Fingerprint>) + (m_allocationPointer - m_mappedFile.const_data()) +
               m_fingerprints.size() * sizeof(std::size_t);
      }

    private:
      const FingerprintStorageType &m_fingerprints;
      unsigned int m_numBits;
      TreeNode *m_tree;

      boost::iostreams::mapped_file m_mappedFile;
      char *m_allocationPointer;
      unsigned long m_nodeCount, m_leafCount;
  };

  template<typename FingerprintStorageType>
  class BitCounter
  {
    public:
      BitCounter(const FingerprintStorageType &fingerprints) : m_fingerprints(fingerprints)
      {
      }

      int count(std::size_t index) const
      {
        return bitCount(m_fingerprints, index);
      }

      void update()
      {
      }

    private:
      const FingerprintStorageType &m_fingerprints;
  };

  template<typename FingerprintStorageType>
  class CachedBitCounter
  {
    public:
      CachedBitCounter(const FingerprintStorageType &fingerprints) : m_fingerprints(fingerprints)
      {
      }

      int count(std::size_t index) const
      {
        return m_bitCounts[index];
      }

      void update()
      {
        m_bitCounts.clear();
        m_bitCounts.reserve(m_fingerprints.size());
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          m_bitCounts.push_back(bitCount(m_fingerprints, i));
      }

    private:
      const FingerprintStorageType &m_fingerprints;
      std::vector<unsigned short> m_bitCounts;
  };

  template<typename FingerprintStorageType, int Dimension>
  class KDGridFingerprintIndex : public FingerprintIndex
  {
      struct Node
      {
      };

      struct TreeNode : public Node
      {
        TreeNode(int size)
        {
          children.resize(size, 0);
        }

        std::vector<Node*> children;
      };

      struct LeafNode : public Node
      {
        std::vector<std::size_t> fingerprints;
      };

      int bitCount(const Fingerprint &fingerprint, int depth) const
      {
        int beginBit = depth * (m_numBits / Dimension);
        int endBit = (depth + 1) * (m_numBits / Dimension);
        
        if (depth == Dimension - 1)
          endBit = (depth + 1) * (m_numBits / Dimension) + (m_numBits % ((depth + 1) * (m_numBits / Dimension)));

        int bitCount = 0;
        for (int i = beginBit; i < endBit; ++i)
          if (fingerprint[i])
            ++bitCount;

        return bitCount;
      }

      int childSize() const
      {
        return m_numBits / Dimension + 1;
      }

      std::vector<std::size_t>& findLeaf(const Fingerprint &fingerprint)
      {
        int depth = 0;
        TreeNode *node = m_tree;

        while (depth < Dimension) {
          int count = bitCount(fingerprint, depth);

          if (depth == Dimension - 1) {
            LeafNode *leaf = static_cast<LeafNode*>(node->children[count]);
            if (!leaf) {
              leaf = new LeafNode;
              node->children[count] = leaf;
            }
            return leaf->fingerprints;            
          }

          TreeNode *child = static_cast<TreeNode*>(node->children[count]);
          if (!child) {
            child = new TreeNode(childSize());
            node->children[count] = child;
          }
        
          node = child;
          ++depth;
        }
      }

      void tanimotoDFS(const FingerprintWord *fingerprint, double threshold, std::size_t maxResults,
          TreeNode *node, int depth, int *n_j, int *bitCounts, int bitCount,
          std::vector<std::pair<std::size_t, double> > &hits) const
      {
        // bound on the number of 1-bits in logical and (first part of bitstring)
        int A_i_min = 0;
        for (int i = 0; i < depth; ++i)
          A_i_min += std::min(bitCounts[i], n_j[i]);

        // bound on the number of 1-bits in logical or (first part of bitstring)
        int A_i_max = 0;
        for (int i = 0; i < depth; ++i)
          A_i_max += std::max(bitCounts[i], n_j[i]);

        // bound on last part of bitstring
        int A_i_last = 0;
        for (int i = depth + 1; i < Dimension; ++i)
          A_i_last += bitCounts[i];

        // bitcount of current part
        int A_i = bitCounts[depth];
        // number of bins in the current part
        int N_i = childSize() - 1;

        int n_lower = std::max(static_cast<int>(std::ceil(threshold * (A_i_max + A_i + A_i_last) - (A_i_min + A_i_last))), 0);
        int n_upper = threshold == 0.0 ? N_i : std::min(static_cast<int>(std::floor((A_i_min + A_i + A_i_last - threshold * (A_i_max + A_i_last)) / threshold)), N_i);
        assert(n_lower >= 0);
        assert(n_upper <= N_i);

        ++depth;

        if (depth == Dimension) {
          int bitCountB = std::accumulate(n_j, n_j + Dimension, 0);

          for (int i = n_lower; i < n_upper + 1; ++i) {
            if (!node->children[i])
              continue;
            LeafNode *leaf = static_cast<LeafNode*>(node->children[i]);

            assert(leaf);
            for (std::size_t j = 0; j < leaf->fingerprints.size(); ++j) {
              double S = MolDB::tanimoto(fingerprint, m_fingerprints.pointer(leaf->fingerprints[j]), bitCount, bitCountB + i, m_numBlocks);
              if (S >= threshold)
                hits.push_back(std::make_pair(leaf->fingerprints[j], S));
              if (hits.size() == maxResults)
                return;
            }
          }
          return;
        }

        for (int i = n_lower; i < n_upper + 1; ++i) {
          if (!node->children[i])
            continue;
          n_j[depth - 1] = i;
          tanimotoDFS(fingerprint, threshold, maxResults, static_cast<TreeNode*>(node->children[i]), depth, n_j, bitCounts, bitCount, hits);
          if (hits.size() == maxResults)
            return;
        }
      }
        
      void clearDFS(TreeNode *node, int depth)
      {
        ++depth;

        for (int i = 0; i < node->children.size(); ++i) {
          if (!node->children[i])
            continue;
          if (depth < Dimension)
            clearDFS(static_cast<TreeNode*>(node->children[i]), depth); // recurse
          else
            delete node->children[i]; // delete leaf node
        }
        
        delete node;
      }

      std::size_t memoryUsageDFS(TreeNode *node, int depth) const
      {
        std::size_t size = sizeof(TreeNode) + node->children.size() * sizeof(Node*);
        ++depth;

        for (int i = 0; i < node->children.size(); ++i) {
          if (!node->children[i])
            continue;
          if (depth < Dimension) {
            size += memoryUsageDFS(static_cast<TreeNode*>(node->children[i]), depth);
          } else {
            LeafNode *leaf = static_cast<LeafNode*>(node->children[i]);
            size += sizeof(LeafNode) + leaf->fingerprints.size() * sizeof(std::size_t);
          }
        }

        return size;        
      }


    public:
      KDGridFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_fingerprints(fingerprints), m_numBits(numBits), m_numBlocks(Fingerprint(numBits).num_blocks())
      {
        m_tree = new TreeNode(m_numBits / Dimension + 1);
      }

      std::string name() const
      {
        return make_string("kD Grid Fingerprint Index (k = ", Dimension, ")");
      }

      ~KDGridFingerprintIndex()
      {
        clear();
        delete m_tree;
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        std::vector<std::size_t> hits;
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i) {
          if (isSubsetSuperset(&query[0], m_fingerprints.pointer(i), fingerprint.num_blocks()))
            hits.push_back(i);
          if (hits.size() == maxResults)
            return hits;
        }

        return hits;
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        if (!maxResults)
          maxResults = std::numeric_limits<std::size_t>::max();

        std::vector<FingerprintWord> query(fingerprint.num_blocks());
        boost::to_block_range(fingerprint, &query[0]);

        int count = 0;
        int n_j[Dimension];
        int bitCounts[Dimension];
        for (int i = 0; i < Dimension; ++i) {
          int part = bitCount(fingerprint, i);
          bitCounts[i] = part;
          count += part;
        }

        std::vector<std::pair<std::size_t, double> > hits;
        tanimotoDFS(&query[0], threshold, maxResults, m_tree, 0, n_j, bitCounts, count, hits);
        
        return hits;
      }

      void clear()
      {
        clearDFS(m_tree, 0);
        m_tree = new TreeNode(m_numBits / Dimension + 1);
      }

      bool load(MessageHandler *msg)
      {
        // clear the index;
        clear();
       
        for (std::size_t i = 0; i < m_fingerprints.size(); ++i)
          findLeaf(m_fingerprints[i]).push_back(i);

        return true;
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        findLeaf(fingerprint).push_back(m_fingerprints.size() - 1);
        return true;
      }

      std::size_t memoryUsage() const
      {
        return sizeof(std::vector<Fingerprint>) + memoryUsageDFS(m_tree, 0);
      }

    private:
      const FingerprintStorageType &m_fingerprints;
      unsigned int m_numBits;
      unsigned int m_numBlocks;
      TreeNode *m_tree;
      const std::vector<std::size_t> m_emptyLeaf;
  };

  template<typename FingerprintStorageType, typename SubsetIndexType, typename TanimotoIndexType>
  class MixedFingerprintIndex : public FingerprintIndex
  {
    public:
      MixedFingerprintIndex(const std::string &filename, const std::string &fingerprint, unsigned int numBits, const FingerprintStorageType &fingerprints)
          : m_tanimoto(filename, fingerprint, numBits, fingerprints), m_subset(filename, fingerprint, numBits, fingerprints)
      {
      }

      std::string name() const
      {
        return make_string("Mixed: ", m_subset.name(), " & ", m_tanimoto.name());
      }

      void reserve(std::size_t size)
      {
        m_tanimoto.reserve(size);
        m_subset.reserve(size);
      }

      std::vector<std::size_t> contains(const Fingerprint &fingerprint, std::size_t maxResults = 0) const
      {
        return m_subset.contains(fingerprint, maxResults);
      }

      std::vector<std::pair<std::size_t, double> > tanimoto(const Fingerprint &fingerprint, double threshold, std::size_t maxResults = 0) const
      {
        return m_tanimoto.tanimoto(fingerprint, threshold, maxResults);
      }

      void clear()
      {
        m_tanimoto.clear();
        m_subset.clear();
      }

      bool create(MessageHandler *msg)
      {
        return m_tanimoto.create(msg) && m_subset.create(msg);
      }

      bool load(MessageHandler *msg)
      {
        return m_tanimoto.load(msg) && m_subset.load(msg);
      }

      bool startInsert(MessageHandler *msg)
      {
        return m_tanimoto.startInsert(msg) && m_subset.startInsert(msg);
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        return m_tanimoto.insert(msg, fingerprint) && m_subset.insert(msg, fingerprint);
      }

      bool commit()
      {
        return m_tanimoto.commit() && m_subset.commit();
      }

      bool stopInsert(MessageHandler *msg)
      {
        return m_tanimoto.stopInsert(msg) && m_subset.stopInsert(msg);
      }

      std::size_t memoryUsage() const
      {
        return m_tanimoto.memoryUsage() + m_subset.memoryUsage();
      }

    private:
      TanimotoIndexType m_tanimoto;
      SubsetIndexType m_subset;
  };

}

#endif
