#ifndef MOLDB_STRINGVECTOR_H
#define MOLDB_STRINGVECTOR_H

#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <cassert>

namespace MolDB {

  class StringVector
  {
    public:

      void reserve(std::size_t size, std::size_t averageStringSize)
      {
        m_pool.reserve(size * (averageStringSize + 1));
        m_strings.reserve(size);
      }

      std::size_t size() const
      {
        return m_strings.size();
      }

      std::size_t poolSize() const
      {
        return m_pool.size();
      }

      void push_back(const std::string &str)
      {
        std::size_t index = m_pool.size();
        m_pool.resize(index + str.size() + 1);
        const char *c_str = str.c_str();
        char *p = &m_pool[index];
        strcpy(p, c_str);
        assert(p[str.size()] == '\0');
        assert(strlen(p) == str.size());
        m_strings.push_back(index);
      }

      const char * const operator[](std::size_t index) const
      {
        assert(index < m_strings.size());
        assert(m_strings[index] < m_pool.size());
        assert(&m_pool[m_strings[index]]);
        return &m_pool[m_strings[index]];
      }

      void clear()
      {
        m_pool.clear();
        m_strings.clear();
      }

      const char* back() const
      {
        assert(&m_pool[m_strings.back()]);
        return &m_pool[m_strings.back()];
      }

    private:
      std::vector<char> m_pool;
      std::vector<std::size_t> m_strings;
  };

}

#endif
