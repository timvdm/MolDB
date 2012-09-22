#ifndef MOLDB_FINGERPRINTSTORAGE_H
#define MOLDB_FINGERPRINTSTORAGE_H

#include "fileio.h"
#include "messagehandler.h"
#include "fingerprint.h"

#include <boost/iostreams/device/mapped_file.hpp>

namespace MolDB {

  /**
   * @class FingerprintFileFormat fingerprintindex.h
   *
   * The FingerprintFileFormat represents a file on disk containing fingerprints.
   *
   * @section File Format
   *
   * File Header:
   * @verbatim
   * bytes
   *  1  -   4   magic number (0xABCDEF02)
   *  5  -   8   fingerprint size (number of bits)
   *  9  -  16   number of fingerprints
   * 17  -  24   end of file (eof)
   * 25  - eof   fingerprint entries
   * @endverbatim
   *
   * Fingerprint Entry:
   * @verbatim
   * bytes
   *  1  -   x   fingerprint (x = number of bits / sizeof(FingerprintWord))
   * @endverbatim
   */
  class FingerprintFileFormat
  {
    public:
      /**
       * The file magic number.
       */
      const static unsigned int magicNumber = IO::MagicNumber::FingerprintFileFormat;

      /**
       * Constructor.
       *
       * @param filename The base filename. (database filename minus .moldb)
       * @param fingerprint The fingerprint type.
       * @param numBits The number of bits in the fingerprint.
       */
      FingerprintFileFormat(const std::string &filename, const std::string &fingerprint, unsigned int numBits)
          : m_numBits(numBits), m_numFingerprints(0)
      {
        if (filename.rfind(".moldb") == std::string::npos)
          m_filename = filename + "." + convertFingerprintName(fingerprint) + ".fi.moldb";
        else
          m_filename = filename;
      }

      /**
       * Get the full filename. (e.g. myDatabase.OpenBabel_FP2.fi.moldb)
       *
       * @return The full filename.
       */
      const std::string& filename() const
      {
        return m_filename;
      }

      /**
       * Create the fingerprint file and write the file header.
       *
       * @param msg The massage handler to report errors.
       *
       * @return True if successful.
       */
      bool create(MessageHandler *msg)
      {
        if (!IO::openFile(msg, m_file, m_filename, std::ios_base::out, magicNumber))
          return false;

        m_file.seekp(0);

        unsigned long endOfFile = 2 * sizeof(unsigned int) + 2 * sizeof(unsigned long);
        
        // write the magic number
        IO::write(m_file, magicNumber);
        // write the fingerprint size (i.e. number of bits)
        unsigned long numBits = m_numBits;
        IO::write(m_file, m_numBits);
        // write the database size (i.e. number of fingerprints)
        unsigned long size = 0;
        IO::write(m_file, size);
         // write the end of file position
        IO::write(m_file, endOfFile);
        
        m_file.close();
        return true; 
      }

      /**
       * Load the fingerprints by reading the header only. The number of
       * fingerprints is returned in @p numFingerprints.
       *
       * @param msg The massage handler to report errors.
       * @param numFingerprints Output parameter for returning the number of
       *        fingerprints in the file.
       *
       * @return True if successful.
       */
      bool load(MessageHandler *msg, unsigned long &numFingerprints)
      {
        if (!IO::openFile(msg, m_file, m_filename, std::ios_base::in, magicNumber))
          return false;

        // get the number of bits in the fingerprint
        int numBits = IO::read<unsigned int>(m_file);
        if (numBits != m_numBits) {
          msg->report(Error, make_string("Fingerprint file ", m_filename, " does not contain the correct fingerprints (number of bits incorrect)."));
          return false;
        }

        // get the number of fingerprints
        numFingerprints = IO::read<unsigned long>(m_file);
        m_numFingerprints = numFingerprints;

        m_file.close();

        return true;
      }

      /**
       * Load the fingerprints. All fingerprints in the file are added to @p
       * fingerprints.
       *
       * @param msg The massage handler to report errors.
       * @param fingerprints Output parameter for returning the loaded
       *        fingerprints.
       *
       * @return True if successful.
       */
       bool load(MessageHandler *msg, std::vector<FingerprintWord> &fingerprints)
      {
        if (!IO::openFile(msg, m_file, m_filename, std::ios_base::in, magicNumber))
          return false;

        // get the number of bits in the fingerprint
        int numBits = IO::read<unsigned int>(m_file);
        if (numBits != m_numBits) {
          msg->report(Error, make_string("Fingerprint file ", m_filename, " does not contain the correct fingerprints (number of bits incorrect)."));
          return false;
        }


        // get the number of fingerprints
        unsigned long numFingerprints = IO::read<unsigned long>(m_file);

        // skip the end of file positons
        IO::read<unsigned long>(m_file);
       
        // start loading
        int numBlocks = Fingerprint(numBits).num_blocks();
        fingerprints.reserve(numBlocks * numFingerprints);
        for (std::size_t i = 0; i < numFingerprints; ++i) {
          // read the fingerprint
          for (int j = 0; j < numBlocks; ++j)
            fingerprints.push_back(IO::read<FingerprintWord>(m_file));
        }
        assert(fingerprints.size() == numBlocks * numFingerprints);
        
        m_numFingerprints = numFingerprints;

        m_file.close();
      }

      bool startInsert(MessageHandler *msg)
      {
        if (!IO::openFile(msg, m_file, m_filename, std::ios_base::in | std::ios_base::out, magicNumber))
          return false;

        // skip the number of bits and number of fingerprints
        m_file.seekp(sizeof(unsigned int) + sizeof(unsigned long), std::ios_base::cur);

        // get the end of file positions and set the call seekp
        unsigned long endOfFile = IO::read<unsigned long>(m_file);
        m_file.seekp(endOfFile);

        return true;
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        // write the fingerprint
        FingerprintWord *blocks = new Fingerprint::block_type[fingerprint.num_blocks()];
        boost::to_block_range(fingerprint, blocks);
        for (int i = 0; i < fingerprint.num_blocks(); ++i) {
          FingerprintWord bits = blocks[i];
          // write the bits
          IO::write(m_file, bits);
        }
        delete [] blocks;
        ++m_numFingerprints;

        return true;
      }

      bool commit()
      {
        unsigned long endOfFile = m_file.tellp();
        // skip the magic number & number of bits
        m_file.seekp(2 * sizeof(unsigned int));
        // write the database size (i.e. number of fingerprints)
        unsigned long size = m_numFingerprints;
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
       * The full filename.
       */
      std::string m_filename;
      /**
       * The std::fstream for the file.
       */
      std::fstream m_file;
      /**
       * The number of bits in the fingerprint.
       */
      int m_numBits;
      /**
       * The number of fingerprints in the file.
       */
      std::size_t m_numFingerprints;
  };

  class InMemoryFingerprintStorage
  {
    public:
      InMemoryFingerprintStorage(const std::string &filename, const std::string &fingerprint, int numBits)
          : m_file(filename, fingerprint, numBits), m_fingerprintName(fingerprint), m_numBits(numBits),
          m_numBlocks(Fingerprint(numBits).num_blocks()), m_fingerprint(numBits)
      {
      }

      const std::string& fingerprint() const
      {
        return m_fingerprintName;
      }

      void reserve(std::size_t size)
      {
        m_fingerprints.reserve(size * m_numBlocks);
      }

      int numBits() const
      {
        return m_numBits;
      }

      int numBlocks() const
      {
        return m_numBlocks;
      }

      std::size_t size() const
      {
        return m_fingerprints.size() / m_numBlocks;
      }

      const FingerprintWord* pointer(std::size_t index) const
      {
        std::size_t offset = index * m_numBlocks;
        return &m_fingerprints[offset];
      }

      const Fingerprint& operator[](std::size_t index) const
      {
        const FingerprintWord *data = pointer(index);
        boost::from_block_range(data, data + m_numBlocks, m_fingerprint);
        return m_fingerprint;
      }

      void clear()
      {
        m_fingerprints.clear();
      }

      bool create(MessageHandler *msg)
      {
        return m_file.create(msg);
      }

      bool load(MessageHandler *msg)
      {
        // clear the index;
        clear();
        return m_file.load(msg, m_fingerprints);
      }

      bool startInsert(MessageHandler *msg)
      {
        return m_file.startInsert(msg);
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        // add the fingerprint to the index
        m_fingerprints.resize(m_fingerprints.size() + m_numBlocks);
        boost::to_block_range(fingerprint, &m_fingerprints[m_fingerprints.size() - m_numBlocks]);
        return m_file.insert(msg, fingerprint);
      }

      bool commit()
      {
        return m_file.commit();
      }

      bool stopInsert(MessageHandler *msg)
      {
        return m_file.stopInsert(msg);
      }

      std::size_t memoryUsage() const
      {
        return sizeof(InMemoryFingerprintStorage) + m_fingerprints.size() * sizeof(FingerprintWord);
      }

    private:
      FingerprintFileFormat m_file;
      std::string m_fingerprintName;
      std::vector<FingerprintWord> m_fingerprints;
      int m_numBits;
      int m_numBlocks;
      mutable Fingerprint m_fingerprint;
  };

  class MemoryMappedFingerprintStorage
  {
    public:
      MemoryMappedFingerprintStorage(const std::string &filename, const std::string &fingerprint, int numBits)
          : m_file(filename, fingerprint, numBits), m_fingerprintName(fingerprint), m_numFingerprints(0),
          m_numBits(numBits), m_numBlocks(Fingerprint(numBits).num_blocks()), m_fingerprint(numBits)
      {
      }

      const std::string& fingerprint() const
      {
        return m_fingerprintName;
      }

      void reserve(std::size_t size)
      {
      }

      int numBits() const
      {
        return m_numBits;
      }

      int numBlocks() const
      {
        return m_numBlocks;
      }
   
      std::size_t size() const
      {
        return m_numFingerprints;
      }

      const FingerprintWord* pointer(std::size_t index) const
      {
        std::size_t offset = 2 * sizeof(unsigned int) + 2 * sizeof(unsigned long) +
                             index * m_fingerprint.num_blocks() * sizeof(FingerprintWord);
      
        const FingerprintWord *data = reinterpret_cast<const FingerprintWord*>(m_mappedFile.data() + offset);
        return data;
      }

      const Fingerprint& operator[](std::size_t index) const
      {
        const FingerprintWord *data = pointer(index);
        boost::from_block_range(data, data + m_fingerprint.num_blocks(), m_fingerprint);
        return m_fingerprint;
      }

      void clear()
      {
      }

      std::size_t memoryUsage() const
      {
        return 0;
      }

      bool create(MessageHandler *msg)
      {
        return m_file.create(msg);
      }

      bool load(MessageHandler *msg)
      {
        bool result = m_file.load(msg, m_numFingerprints);
        m_mappedFile.open(m_file.filename());
        return result;
      }

      bool startInsert(MessageHandler *msg)
      {
        m_mappedFile.close();
        return m_file.startInsert(msg);
      }

      bool insert(MessageHandler *msg, const Fingerprint &fingerprint)
      {
        return m_file.insert(msg, fingerprint);
      }

      bool commit()
      {
        return m_file.commit();
      }

      bool stopInsert(MessageHandler *msg)
      {
        m_mappedFile.open(m_file.filename());
        return m_file.stopInsert(msg);
      }

    private:
      FingerprintFileFormat m_file;
      std::string m_fingerprintName;
      unsigned long m_numFingerprints;
      int m_numBits;
      int m_numBlocks;
      boost::iostreams::mapped_file_source m_mappedFile;
      mutable Fingerprint m_fingerprint;
  };

}

#endif
