#ifndef MOLDB_MEMORY_H
#define MOLDB_MEMORY_H

#include <vector>
#include <climits>

namespace MolDB {

  namespace Impl {

    struct Chunk
    {
      void init(std::size_t blockSize, unsigned char blocks)
      {
        m_data = new unsigned char[blockSize * blocks];
        m_firstAvailableBlock = 0;
        m_blocksAvailable = blocks;
        unsigned char *p = m_data;
        // store pointer to next available block inside the blocks
        for (unsigned char i = 0; i != blocks; p += blockSize)
          *p = ++i;
      }

      void* allocate(std::size_t blockSize)
      {
        if (!m_blocksAvailable)
          return 0;
        unsigned char *result = m_data + (m_firstAvailableBlock * blockSize);
        // update m_firstAvailableBlock to point to the next block
        m_firstAvailableBlock = *result;
        --m_blocksAvailable;
        return result;
      }

      void deallocate(void *p, std::size_t blockSize)
      {
        assert(p >= m_data);
        unsigned char *toRelease = static_cast<unsigned char*>(p);
        // alignment check
        assert((toRelease - m_data) % blockSize == 0);
        *toRelease = m_firstAvailableBlock;
        m_firstAvailableBlock = static_cast<unsigned char>((toRelease - m_data) / blockSize);
        // turncation check
        assert(m_firstAvailableBlock == (toRelease - m_data) / blockSize);
        ++m_blocksAvailable;
      }

      unsigned char *m_data;
      unsigned char m_firstAvailableBlock;
      unsigned char m_blocksAvailable;
    };

    template<std::size_t BlockSize, unsigned char NumBlocks>
    class FixedAllocator
    {
      public:
        FixedAllocator() : m_allocChunk(0), m_deallocChunk(0)
        {
        }

        void* allocate()
        {
          if (m_allocChunk == 0 || m_allocChunk->m_blocksAvailable == 0) {
            //std::cout << "no blocks available..." << std::endl;
            // no available memory in this chunk, try to find one
            for (std::vector<Chunk>::iterator i = m_chunks.begin(); ; ++i) {
              if (i == m_chunks.end()) {
                //std::cout << "creating a chunke..." << std::endl;
                // all filled up, add a new chunk
                m_chunks.reserve(m_chunks.size() + 1);
                Chunk newChunk;
                newChunk.init(BlockSize, NumBlocks);
                m_chunks.push_back(newChunk);
                m_allocChunk = &m_chunks.back();
                m_deallocChunk = &m_chunks.back();
                break;
              }
              if (i->m_blocksAvailable > 0) {
                //std::cout << "found a chunke..." << std::endl;
                // found a chunk
                m_allocChunk = &*i;
                break;
              }
            }
          }
          assert(m_allocChunk != 0);
          assert(m_allocChunk->m_blocksAvailable > 0);
          return m_allocChunk->allocate(BlockSize);
        }

        void deallocate(void *p)
        {
          if (m_deallocChunk) {
            if (isInChunk(m_deallocChunk, p)) {
              m_deallocChunk->deallocate(p, BlockSize);
              m_allocChunk = m_deallocChunk;
              return;
            }
          }
          for (std::vector<Chunk>::iterator i = m_chunks.begin(); i != m_chunks.end(); ++i) {
            if (isInChunk(&*i, p)) {
              i->deallocate(p, BlockSize);
              m_deallocChunk = &*i;
            }
          }
        }

        static FixedAllocator<BlockSize, NumBlocks>* instance()
        {
          static FixedAllocator<BlockSize, NumBlocks> *p = 0;
          if (!p)
            p = new FixedAllocator<BlockSize, NumBlocks>;
          return p;
        }

      private:
        bool isInChunk(const Chunk *chunk, void *p) const
        {
          return (p >= chunk->m_data && p < chunk->m_data + BlockSize * NumBlocks);
        }

        std::vector<Chunk> m_chunks;
        Chunk *m_allocChunk;
        Chunk *m_deallocChunk;
    };

    template<std::size_t BlockSize, std::size_t NumBlocks>
    class FixedAllocator2
    {
      public:
        FixedAllocator2() : m_allocIndex(0)/*, m_deallocIndex(0)*/
        {
          m_chunks.resize(NumBlocks);
          /*
          for (std::size_t i = 0; i < NumBlocks; ++i)
            m_chunks[i].index = i + 1;
          */
        }

        void* allocate()
        {
          assert(m_allocIndex < m_chunks.size());
          //++m_allocIndex = m_cunks[m_allocIndex].index;
          return &m_chunks[++m_allocIndex];
        }

        void deallocate(void *p)
        {
          /*
          for (std::size_t i = 0; i < m_chunks.size(); ++i) {
            if (p == static_cast<void*>(&m_chunks[i])) {
              static_cast<Chunk2*>(&m_chunks[i])->index = i

          }
          */
        }

        static FixedAllocator2<BlockSize, NumBlocks>* instance()
        {
          static FixedAllocator2<BlockSize, NumBlocks> *p = 0;
          if (!p)
            p = new FixedAllocator2<BlockSize, NumBlocks>;
          return p;
        }

      private:
        struct Chunk2
        {
          //union {
          //  std::size_t index;
            char data[BlockSize];
          //}
        };

        std::vector<Chunk2> m_chunks;
        std::size_t m_allocIndex;
        //std::szie_t m_deallocIndex;
    };


  }

  template<typename Object>
  class SmallObject
  {
    public:
      static void* operator new(std::size_t size)
      {
        //struct Helper { unsigned char data[sizeof(Object)]; };
        //return static_cast<void*>(new Helper);
        Impl::FixedAllocator<sizeof(Object), static_cast<unsigned char>(UCHAR_MAX)>::instance()->allocate();
        //Impl::FixedAllocator2<sizeof(Object), static_cast<std::size_t>(5000000)>::instance()->allocate();
      }

      static void operator delete(void *p, std::size_t size)
      {
        Impl::FixedAllocator<sizeof(Object), static_cast<unsigned char>(UCHAR_MAX)>::instance()->deallocate(p);
      }

      /*
      virtual ~SmallObject()
      {
      }
      */
  };
    

}

#endif
