#ifndef MOLDB_DATABASE_SMILESINDEX_H
#define MOLDB_DATABASE_SMILESINDEX_H

#include "stringvector.h"
#include "fileio.h"

#include <map>

namespace MolDB {
  
  /**
   * Compare function object for ordered std::set.
   */
  struct SmilesCompareLess
  {
    bool operator()(const char *smiles1, const char *smiles2) const
    {
      assert(smiles1 && smiles2);
      return strcmp(smiles1, smiles2) < 0;
    }
  };

  /**
   * Compare function object for unordered std::tr1::unordered_set and boost::unordered_set.
   */
  struct SmilesCompareEqual
  {
    bool operator()(const char *smiles1, const char *smiles2) const
    {
      return strcmp(smiles1, smiles2) == 0;
    }
  };


  /**
   *
   *
   * @section File Format
   *
   * File Header:
   * @verbatim
   * bytes
   *  1  -   4   magic number (0xABCDEF01)
   *  5  -  12   number of smiles
   * 13  -  20   end of file (eof)
   * 21  - eof   smiles entries
   * @endverbatim
   *
   * Smiles Entry:
   * @verbatim
   * bytes
   *  1  -   4   smiles size (number of characters)
   *  5  -   x   smiles string (x = 4 + smiles size)
   * @endverbatim
   *
   */
  template<typename MapType = std::map<const char*, unsigned long, SmilesCompareLess> >
  class InMemoryMapSmilesIndex
  {
    public:
      const static unsigned int magicNumber = IO::MagicNumber::InMemoryMapSmilesIndex;

      InMemoryMapSmilesIndex(const std::string &filename)
          : m_filename(filename)
      {
      }

      void reserve(std::size_t numberOfSmiles, std::size_t averageSmilesSize)
      {
        m_smiles.reserve(numberOfSmiles, averageSmilesSize);
      }

      /**
       * Check if the index contains @p smiles.
       */
      std::size_t contains(const std::string &smiles) const
      {
        typename MapType::const_iterator i = m_smilesMap.find(smiles.c_str());
        if (i == m_smilesMap.end())
          return std::numeric_limits<std::size_t>::max();
        return i->second;
      }

      /**
       * Get the number of smiles in the index.
       */
      std::size_t size() const
      {
        return m_smiles.size();
      }

      /**
       * Get the smiles with index @p index.
       */
      const char* get(std::size_t index) const
      {
        return m_smiles[index];
      }

      /**
       * Clear the index.
       */
      void clear()
      {
        m_smiles.clear();
        m_smilesMap.clear();
      }

      std::size_t memoryUsage() const
      {
        return sizeof(StringVector) + sizeof(MapType) +
               sizeof(unsigned char*) * size() + m_smiles.poolSize();
      }

      bool create(MessageHandler *msg)
      {
        if (!IO::openFile(msg, m_file, m_filename + ".imssi.moldb", std::ios_base::out, magicNumber))
          return false;

        m_file.seekp(0);

        unsigned long endOfFile = sizeof(unsigned int) + 2 * sizeof(unsigned long);
        
        // write the magic number
        IO::write(m_file, magicNumber);
        // write the database size (i.e. number of smiles)
        unsigned long size = 0;
        IO::write(m_file, size);
         // write the end of file position
        IO::write(m_file, endOfFile);
        
        m_file.close();
        return true; 
      }

      bool load(MessageHandler *msg)
      {
        if (!IO::openFile(msg, m_file, m_filename + ".imssi.moldb", std::ios_base::in, magicNumber))
          return false;

        // clear the index
        clear();

        // get the number of smiles
        unsigned long numSmiles = IO::read<unsigned long>(m_file);

        // skip the end of file positons
        IO::read<unsigned long>(m_file);
        
        unsigned int BUFFER_SIZE = 1000;
        char *buffer = new char[BUFFER_SIZE];
        
        // start loading
        for (std::size_t i = 0; i < numSmiles; ++i) {
          // read the smiles size
          unsigned int smilesSize = IO::read<unsigned int>(m_file);
          if (smilesSize > BUFFER_SIZE) {
            delete [] buffer;
            buffer = new char[smilesSize];
            BUFFER_SIZE = smilesSize;
          }
          // read the smiles
          m_file.read(buffer, smilesSize);
          // add the smiles to the index
          if (m_smiles.size()) {
            const char *p = m_smiles[0];
            m_smiles.push_back(std::string(buffer, smilesSize));
            if (p != m_smiles[0]) {
              m_smilesMap.clear();
              for (std::size_t j = 0; j < m_smiles.size(); ++j)
                m_smilesMap[m_smiles[j]] = j;
            }
          } else
            m_smiles.push_back(std::string(buffer, smilesSize));
          m_smilesMap[m_smiles.back()] = i;
        }
        
        delete [] buffer;
        m_file.close();
      }

      bool startInsert(MessageHandler *msg)
      {
        if (!IO::openFile(msg, m_file, m_filename + ".imssi.moldb", std::ios_base::in | std::ios_base::out, magicNumber))
          return false;

        // skip the number of smiles
        m_file.seekp(sizeof(unsigned long), std::ios_base::cur);

        // get the end of file positions and set the call seekp
        unsigned long endOfFile = IO::read<unsigned long>(m_file);
        m_file.seekp(endOfFile);

        return true;
      }

      bool insert(MessageHandler *msg, const std::string &smiles)
      {
        // write the smiles size
        unsigned int smilesSize = smiles.size();
        IO::write(m_file, smilesSize);
        // write the smiles itself
        m_file.write(smiles.c_str(), smilesSize);
        // add the smiles to the index
        if (m_smiles.size()) {
          const char *p = m_smiles[0];
          m_smiles.push_back(smiles);
          if (p != m_smiles[0]) {
            m_smilesMap.clear();
            for (std::size_t i = 0; i < m_smiles.size(); ++i)
              m_smilesMap[m_smiles[i]] = i;
          }
        } else
          m_smiles.push_back(smiles);

        m_smilesMap[m_smiles.back()] = m_smiles.size() - 1;
        return true;
      }

      bool commit()
      {
        unsigned long endOfFile = m_file.tellp();
        // skip the magic number
        unsigned long size = m_smiles.size();
        m_file.seekp(sizeof(unsigned int));
        // write the database size (i.e. number of smiles)
        IO::write(m_file, size);
        // write the end of file position
        IO::write(m_file, endOfFile);
        // reset position
        m_file.seekp(endOfFile);
        return true;
      }

      bool stopInsert(MessageHandler *msg)
      {
        // close the file stream
        m_file.close();
        return true;
      }

    private:

      /**
       * The base database filename.
       */
      std::string m_filename;
      /**
       * The smiles index file streams.
       */
      std::fstream m_file;
      /**
       * A StringVector containing the smiles strings.
       */
      StringVector m_smiles;
      /**
       * A std::map to quickly find a smiles. std::map is almost always
       * implemented as a red-black binary search tree and has an O(log(n))
       * run time complexity for finding a smiles.
       */
      MapType m_smilesMap;
  };























































































































  class MemoryMappedRedBlackTreeSmilesIndex
  {
    private:
      const static unsigned int magicNumber = IO::MagicNumber::InMemoryMapSmilesIndex;

      bool openFiles(MessageHandler *msg, std::ios_base::openmode mode, const std::string &function)
      {
        // open the positions file
        std::string positionsFilename = m_filename + ".sip.moldb";
        m_positionsFile.open(positionsFilename.c_str(), std::ios_base::binary | mode);
        if (!m_positionsFile) {
          msg->report(Error, make_string("Could not open database file ", positionsFilename, "."));
          return false;
        }

        // TODO check magic numbers
        if (mode & std::ios_base::in)
          m_positionsFile.seekg(sizeof(unsigned int));
        if (mode & std::ios_base::out)
          m_positionsFile.seekp(sizeof(unsigned int));

        // open the smiles file
        std::string smilesFilename = m_filename + ".sis.moldb";
        m_smilesFile.open(smilesFilename.c_str(), std::ios_base::binary | mode);
        if (!m_smilesFile) {
          msg->report(Error, make_string("Could not open database file ", smilesFilename, "."));
          return false;
        }

        // TODO check magic numbers
        if (mode & std::ios_base::in)
          m_smilesFile.seekg(sizeof(unsigned int));
        if (mode & std::ios_base::out)
          m_smilesFile.seekp(sizeof(unsigned int));

        return true;
      }

    public:
      MemoryMappedRedBlackTreeSmilesIndex(const std::string &filename)
          : m_filename(filename)
      {
      }

      void reserve(std::size_t numberOfSmiles, std::size_t averageSmilesSize)
      {
        m_smiles.reserve(numberOfSmiles, averageSmilesSize);
      }

      /**
       * Check if the index contains @p smiles.
       */
      bool contains(const std::string &smiles) const
      {
        return m_smilesMap.find(smiles.c_str()) != m_smilesMap.end();
      }

      /**
       * Get the number of smiles in the index.
       */
      std::size_t size() const
      {
        return m_smiles.size();
      }

      /**
       * Get the smiles with index @p index.
       */
      const char* get(std::size_t index) const
      {
        return m_smiles[index];
      }

      /**
       * Clear the index.
       */
      void clear()
      {
        m_smiles.clear();
        m_smilesMap.clear();
      }

      std::size_t memoryUsage() const
      {
        return sizeof(StringVector) + sizeof(std::map<const char*, bool, CompareSmiles>) +
               sizeof(unsigned char*) * size() + sizeof(bool) * size() + m_smiles.poolSize();
      }

      bool create(MessageHandler *msg)
      {
        if (!openFiles(msg, std::ios_base::out, "create"))
          return false;

        m_smilesFile.seekp(0);
        m_positionsFile.seekp(0);

        unsigned long endOfFile = sizeof(unsigned int) + 2 * sizeof(unsigned long);
        
        // write the magic number
        IO::write(m_positionsFile, magicNumber);
        // write the database size (i.e. number of smiles)
        unsigned long size = 0;
        IO::write(m_positionsFile, size);
        // write the end of file position
        IO::write(m_positionsFile, endOfFile);

        // write the magic number
        IO::write(m_smilesFile, magicNumber);
        // write the database size (i.e. number of smiles)
        IO::write(m_smilesFile, size);
         // write the end of file position
        IO::write(m_smilesFile, endOfFile);
        
        m_positionsFile.close();
        m_smilesFile.close();
        return true; 
      }

      bool load(MessageHandler *msg)
      {
        if (!openFiles(msg, std::ios_base::in | std::ios_base::out, "load"))
          return false;

        // get the sizes
        unsigned long numPositions = IO::read<unsigned long>(m_positionsFile);
        unsigned long numSmiles = IO::read<unsigned long>(m_smilesFile);

        // check if the sizes match
        if (numPositions != numSmiles) {
          msg->report(Error, make_string("Database inconsistency: positions file contains ", numPositions, " and the smiles file ", numSmiles, "."));
          return false;
        }

        // skip the end of file positons
        IO::read<unsigned long>(m_positionsFile);
        IO::read<unsigned long>(m_smilesFile);
        
        unsigned int BUFFER_SIZE = 1000;
        char *buffer = new char[BUFFER_SIZE];
        
        // start loading
        for (std::size_t i = 0; i < numSmiles; ++i) {
          // read position
          IO::read<unsigned long>(m_positionsFile);
          // read the smiles size
          unsigned int smilesSize = IO::read<unsigned int>(m_smilesFile);
          if (smilesSize > BUFFER_SIZE) {
            delete [] buffer;
            buffer = new char[smilesSize];
            BUFFER_SIZE = smilesSize;
          }
          // read the smiles
          m_smilesFile.read(buffer, smilesSize);
          // add the smiles to the index
          m_smiles.push_back(std::string(buffer, smilesSize));
          m_smilesMap[m_smiles.back()] = true;
        }
        
        delete [] buffer;
        m_positionsFile.close();
        m_smilesFile.close();
      }

      bool startInsert(MessageHandler *msg)
      {
        if (!openFiles(msg, std::ios_base::in | std::ios_base::out, "startInsert"))
          return false;

        // get the sizes
        unsigned long positionsSize = IO::read<unsigned long>(m_positionsFile);
        unsigned long smilesSize = IO::read<unsigned long>(m_smilesFile);

        // check if the sizes match
        if (positionsSize != smilesSize) {
          msg->report(Error, make_string("Database inconsistency: positions file contains ", positionsSize, " and the smiles file ", smilesSize, "."));
          return false;
        }

        // get the end of file positions and set the call seekp
        unsigned long endOfFile = IO::read<unsigned long>(m_positionsFile);
        m_positionsFile.seekp(endOfFile);
        endOfFile = IO::read<unsigned long>(m_smilesFile);
        m_smilesFile.seekp(endOfFile);

        return true;
      }

      bool insert(MessageHandler *msg, const std::string &smiles)
      {
        // write the position in the smiles file to the positions file
        unsigned long position = m_smilesFile.tellp();
        IO::write(m_positionsFile, position);
        // write the smiles size
        unsigned int smilesSize = smiles.size();
        IO::write(m_smilesFile, smilesSize);
        // write the smiles itself
        m_smilesFile.write(smiles.c_str(), smilesSize);
        // add the smiles to the index
        m_smiles.push_back(smiles);
        m_smilesMap[m_smiles.back()] = true;
        return true;
      }

      bool commit()
      {
        unsigned long endOfFile = m_positionsFile.tellp();
        // skip the magic number
        m_positionsFile.seekp(sizeof(unsigned int));
        // write the database size (i.e. number of smiles)
        unsigned long size = m_smiles.size();
        IO::write(m_positionsFile, size);
        // write the end of file position
        IO::write(m_positionsFile, endOfFile);
        // reset position
        m_positionsFile.seekp(endOfFile);

        endOfFile = m_smilesFile.tellp();
        // skip the magic number
        m_smilesFile.seekp(sizeof(unsigned int));
        // write the database size (i.e. number of smiles)
        IO::write(m_smilesFile, size);
        // write the end of file position
        IO::write(m_smilesFile, endOfFile);
        // reset position
        m_smilesFile.seekp(endOfFile);
        return true;
      }

      bool stopInsert(MessageHandler *msg)
      {
        // close the file streams
        m_positionsFile.close();
        m_smilesFile.close();
        return true;
      }

    private:
      /**
       * Compare function object for std::map.
       */
      struct CompareSmiles
      {
        bool operator()(const char *smiles1, const char *smiles2) const
        {
          //if (!smiles1 || !smiles2)
          //  return 0;
          return strcmp(smiles1, smiles2) < 0;
        }
      };

      /**
       * The base database filename.
       */
      std::string m_filename;
      /**
       * The smiles index file streams.
       */
      std::fstream m_positionsFile, m_smilesFile;
      /**
       * A StringVector containing the smiles strings.
       */
      StringVector m_smiles;
      /**
       * A std::map to quickly find a smiles. std::map is almost always
       * implemented as a red-black binary search tree and has an O(log(n))
       * run time complexity for finding a smiles.
       */
      std::map<const char*, bool, CompareSmiles> m_smilesMap;
  };

}

#endif
