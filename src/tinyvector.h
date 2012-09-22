#ifndef MOLDB_TINYVECTOR_H
#define MOLDB_TINYVECTOR_H

#include <vector>
#include <limits>
#include <cassert>

#include <iostream> // FIXME

#include "typetraits.h"

//#define DEBUG_MOLDB_TINYVECTOR_H

namespace MolDB {

  namespace Impl {

    template<typename T, typename UsageType>
    class TinyVectorMemoryPool
    {
        typedef TinyVectorMemoryPool<T, UsageType> MemoryPool;
        typedef typename TypeTraits::Select<TypeTraits::IsPointer<T>::result, T, T*>::result PointerType;
      public:
        std::size_t size() const
        {
          return m_pool.size();
        }

        std::size_t capacity() const
        {
          return m_pool.capacity();
        }

        T* allocate(std::size_t size)
        //PointerType allocate(std::size_t size)
        {
          std::size_t poolSize = m_pool.size();
          //std::cout << "TinyVectorMemoryPool::allocate(" << size << ")" << std::endl;
          //std::cout << " pool size = " << m_pool.size() << ", size = " << size << std::endl;
          m_pool.resize(poolSize + size);
          //std::cout << "pool resized" << std::endl;
          //T *p = new (&m_pool[poolSize]) T;
          T *p = &m_pool[poolSize];
          //std::cout << "pointer obtained" << std::endl;
          //PointerType p = &m_pool[m_pool.size()];
          for (std::size_t i = 0; i < size; ++i)
            p[i] = invalidT();
          //std::cout << "pointers initialized" << std::endl;
          return p;
        }

        void deallocate()
        {
          m_pool.clear();
        }

        static MemoryPool& instance(std::size_t size = 0)
        {
          static MemoryPool *p = 0;
          if (!p)
            p = new MemoryPool;
          if (size)
            p->m_pool.reserve(size);
          return *p;
        }

        const std::vector<T>& pool() const
        {
          return m_pool;
        }

      private:
        TinyVectorMemoryPool()
        {
        }

        T invalidT() const
        {
          // return max value for integers or 0 for pointers
          return std::numeric_limits<T>::max();
        }

        std::vector<T> m_pool;
        TypeTraits::StaticAssert<TypeTraits::IsInteger<T>::result || TypeTraits::IsPointer<T>::result> TinyVector_T_must_be_an_integer_or_pointer;
    };

  }

  template<typename T, typename UsageType, typename MemoryPoolType = Impl::TinyVectorMemoryPool<T, UsageType> >
  class TinyVector
  {
      typedef TinyVector<T, UsageType, MemoryPoolType> VectorType;
      typedef typename TypeTraits::Select<TypeTraits::IsPointer<T>::result, T, T*>::result PointerType;
    public:
      TinyVector(std::size_t size) : m_elements(MemoryPoolType::instance().allocate(size + 1))
      {
        TypeTraits::StaticAssert<TypeTraits::IsInteger<T>::result || TypeTraits::IsPointer<T>::result>();
#ifdef DEBUG_MOLDB_TINYVECTOR_H
        m_size = vectorSize;
#endif
      }

      void push_back(T element)
      {
#ifdef DEBUG_MOLDB_TINYVECTOR_H
        assert(size() < m_size + 1);
#endif
        std::size_t index = 0;
        while (m_elements[index] != invalidT())
          ++index;
        assert(m_elements[index + 1] == invalidT());
        m_elements[index] = element;
      }

      std::size_t size() const
      {
        assert(m_elements);
        std::size_t index = 0;
        while (m_elements[index] != invalidT())
          ++index;
        return index;
      }

      const T& operator[](std::size_t index) const
      {
#ifdef DEBUG_MOLDB_TINYVECTOR_H
        assert(index < m_size);
#endif
        return m_elements[index];
      }

      T& operator[](std::size_t index)
      {
#ifdef DEBUG_MOLDB_TINYVECTOR_H
        assert(index < m_size);
#endif
        return m_elements[index];
      }

      const T& back() const
      {
        return m_elements[size() - 1];
      }

      T& back()
      {
        return m_elements[size() - 1];
      }

    private:
      T invalidT() const
      {
        // return max value for integers or 0 for pointers
        return std::numeric_limits<T>::max();
      }

      T* m_elements;
      //PointerType m_elements;
#ifdef DEBUG_MOLDB_TINYVECTOR_H
      std::size_t m_size;
#endif
  };

}

#endif
