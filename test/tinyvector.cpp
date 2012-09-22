#include "../src/tinyvector.h"

#include "test.h"

using namespace MolDB;

struct SomeUsage {};
struct SomeType {};

int main()
{
  //
  // primitive type
  //
  Impl::TinyVectorMemoryPool<int, SomeUsage>::instance(1000);
  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int, SomeUsage>::instance().capacity()), 1000);
  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int, SomeUsage>::instance().size()), 0);

  TinyVector<int, SomeUsage> v(5);
  MOLDB_COMPARE(v.size(), 0);
  v.push_back(42);
  MOLDB_COMPARE(v.size(), 1);
  MOLDB_COMPARE(v[0], 42);
  MOLDB_COMPARE(v[1], std::numeric_limits<int>::max());
  MOLDB_COMPARE(v[5], std::numeric_limits<int>::max());
  
  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int, SomeUsage>::instance().size()), 6);
  
  //
  // pointer
  // 
  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int*, SomeUsage>::instance().capacity()), 0);
  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int*, SomeUsage>::instance().size()), 0);
  const std::vector<int*> &pool = Impl::TinyVectorMemoryPool<int*, SomeUsage>::instance().pool();

  TinyVector<int*, SomeUsage> v2(5);;
  MOLDB_COMPARE(v2.size(), 0);
  v2.push_back(reinterpret_cast<int*>(42));
  v2.push_back(reinterpret_cast<int*>(16));
  MOLDB_COMPARE(v2.size(), 2);
  MOLDB_COMPARE(v2[0], reinterpret_cast<int*>(42));
  MOLDB_COMPARE(v2[1], reinterpret_cast<int*>(16));
  MOLDB_COMPARE(v2[5], reinterpret_cast<int*>(0)); 
 
  MOLDB_COMPARE(pool[0], reinterpret_cast<int*>(42));
  MOLDB_COMPARE(pool[1], reinterpret_cast<int*>(16));
  MOLDB_COMPARE(pool[2], reinterpret_cast<int*>(0));
  MOLDB_COMPARE(pool[3], reinterpret_cast<int*>(0));
  MOLDB_COMPARE(pool[4], reinterpret_cast<int*>(0));
  MOLDB_COMPARE(pool[5], reinterpret_cast<int*>(0));

  MOLDB_COMPARE((Impl::TinyVectorMemoryPool<int*, SomeUsage>::instance().size()), 6);

  std::cout << "pool for TinyVector<int*>: ";
  for (int i = 0; i < pool.size(); ++i)
    std::cout << pool[i] << " ";
  std::cout << std::endl;



  // correctly produces compile error
  //TinyVector<SomeType, SomeUsage> v3(5);;

}
