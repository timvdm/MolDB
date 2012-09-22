#ifndef MOLDB_FILEIO_H
#define MOLDB_FILEIO_H

#include "typetraits.h"
#include "messagehandler.h"
#include "util.h"

#include <fstream>
#include <cassert>

namespace MolDB {

  namespace IO {

    struct MagicNumber
    {
      static const unsigned int Database = 0xABCDEF00;
      static const unsigned int InMemoryMapSmilesIndex = 0xABCDEF01;
      static const unsigned int FingerprintFileFormat = 0xABCDEF02;
      static const unsigned int FingerprintNodeBitsFileFormat = 0xABCDEF03;
    };

    template<typename FileStream>
    bool openFile(MessageHandler *msg, FileStream &fs, const std::string &filename, std::ios_base::openmode mode, unsigned int magicNumber)
    {
      // open the smiles file
      fs.open(filename.c_str(), std::ios_base::binary | mode);
      if (!fs) {
        msg->report(Error, make_string("Could not open database file ", filename, "."));
        return false;
      }

      // TODO check magic numbers
      if (mode & std::ios_base::in)
        fs.seekg(sizeof(unsigned int));
      if (mode & std::ios_base::out)
        fs.seekp(sizeof(unsigned int));

      return true;
    }


    template<typename OutStream, typename T>
    void write(OutStream &os, T t)
    {
      os.write(reinterpret_cast<char*>(&t), sizeof(T));
      assert(!os.bad());
    }

    template<typename T, typename InStream>
    T read(InStream &is)
    {
      T t;
      is.read(reinterpret_cast<char*>(&t), sizeof(T));
      return t;
    }

  }

}

#endif
