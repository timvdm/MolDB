#ifndef MOLDB_FINGERPRINT_H
#define MOLDB_FINGERPRINT_H

#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

#include <boost/dynamic_bitset.hpp>

namespace MolDB {

  //typedef unsigned long FingerprintWord;
  typedef unsigned long FingerprintWord;
  typedef boost::dynamic_bitset<FingerprintWord> Fingerprint;

  /**
   * OpenBabel::FP2 -> OpenBabel_FP2
   */
  inline std::string convertFingerprintName(const std::string &fingerprintName)
  {
    std::string fp = fingerprintName;
    fp.replace(fingerprintName.find("::"), 2, "_");
    return fp;
  }

  inline bool isSubsetSuperset(const FingerprintWord *fingerprint1, const FingerprintWord *fingerprint2, int numBlocks)
  {
    for (int i = 0; i < numBlocks; ++i)
      if (fingerprint1[i] & ~fingerprint2[i])
        return false;
    return true; 
  }

  inline double tanimoto(const Fingerprint &fingerprint1, const Fingerprint &fingerprint2)
  {
    return static_cast<double>((fingerprint1 & fingerprint2).count()) / (fingerprint1 | fingerprint2).count();
  }

  inline double tanimoto(const FingerprintWord *fingerprint1, const FingerprintWord *fingerprint2, int numBlocks)
  {
    int andCount = 0;
    int orCount = 0;
    
    for (int i = 0; i < numBlocks; ++i) {
      long andfp = fingerprint1[i] & fingerprint2[i];
      long orfp = fingerprint1[i] | fingerprint2[i];
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      //andCount += __builtin_popcountl(andfp);
      //orCount += __builtin_popcountl(orfp);
      andCount += __builtin_popcountll(andfp);
      orCount += __builtin_popcountll(orfp);
#else
      for(; andfp; andfp = andfp << 1)
        if(andfp < 0)
          ++andbits;
      for(; orfp; orfp = orfp << 1)
        if(orfp < 0)
          ++orbits;
#endif
    }

    return static_cast<double>(andCount) / orCount;
  }

  inline double tanimoto(const FingerprintWord *fingerprint1, const FingerprintWord *fingerprint2, int bitCount1, int bitCount2, int numBlocks)
  {
    int andCount = 0;
    
    for (int i = 0; i < numBlocks; ++i) {
      long andfp = fingerprint1[i] & fingerprint2[i];
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      andCount += __builtin_popcountl(andfp);
#else
      for(; andfp; andfp = andfp << 1)
        if(andfp < 0)
          ++andbits;
#endif
    }

    return static_cast<double>(andCount) / (bitCount1 + bitCount2 - andCount);
  }

  template<typename FingerprintStorageType>
  inline int bitCount(const FingerprintStorageType &fingerprints, std::size_t index)
  {
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
    // use popcount
    int count = 0;
    int numBlocks = fingerprints.numBlocks();
    const FingerprintWord *p = fingerprints.pointer(index);

    for (int i = 0; i < numBlocks; ++i)
      count += __builtin_popcountl(*(p + i));

    return count;
#else
    // use boost's 64-bit lookup table count method
    return m_fingerprints[index].count();
#endif
  }


  class FingerprintGenerator
  {
      struct Pattern
      {
        Pattern(int bit_, int count_, const std::string &smarts_) : bit(bit_), count(count_), smarts(smarts_)
        {
          sp.Init(smarts);
        }

        int bit;
        int count;
        std::string smarts;
        OpenBabel::OBSmartsPattern sp;
      };

      void load(const std::string &filename)
      {
        std::ifstream ifs(filename.c_str());

        int bit = 0;
        std::string line;
        while (std::getline(ifs, line)) {
          if (line.empty())
            continue;

          std::size_t pos = line.find(" ");
          std::string smarts = line.substr(0, pos);
          std::stringstream ss(line.substr(pos + 1, line.find(" ", pos + 1)));
          int count;
          ss >> count;

          m_patterns.push_back(Pattern(bit, count + 1, smarts));
          ++bit;
        }
      }
      FingerprintGenerator(const std::string &patternsFilename)
      {
        load(patternsFilename);
      }
    public:
      Fingerprint fingerprint(OpenBabel::OBMol *mol)
      {
        Fingerprint fp(m_patterns.size());
        for (std::size_t i = 0; i < m_patterns.size(); ++i) {
          Pattern &pattern = m_patterns[i];

          bool value;
          if (pattern.count == 1) {
            value = pattern.sp.Match(*mol, true);          
          } else {
            pattern.sp.Match(*mol);
            value = pattern.sp.GetUMapList().size() >= pattern.count;
          }

          fp.set(pattern.bit, value);
        }

        return fp;
      }

      static FingerprintGenerator* instance(const std::string &patternsFilename)
      {
        static std::map<std::string, FingerprintGenerator*> m_instances;
        if (m_instances[patternsFilename] == 0) {
          m_instances[patternsFilename] = new FingerprintGenerator(patternsFilename);
        }
        return m_instances[patternsFilename];
      }
    private:
      std::vector<Pattern> m_patterns;
  };


}

#endif
