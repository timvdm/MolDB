#ifndef MOLDB_UTIL_H
#define MOLDB_UTIL_H

#include <sstream>

namespace MolDB {

  template<typename T1>
  inline std::string make_string(const T1 &t1)
  {
    std::stringstream ss;
    ss << t1;
    return ss.str();
  }

  template<typename T1, typename T2>
  inline std::string make_string(const T1 &t1, const T2 &t2)
  {
    std::stringstream ss;
    ss << t1 << t2;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const T10 &t10)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << t10;
    return ss.str();
  }

  inline std::string progressBar(std::size_t value, std::size_t maximum)
  {
    int percent = (maximum > 100) ? value / (maximum / 100) : value * (100 / maximum);
    if (value == maximum)
      percent = 100;
    std::string bar = "[";
    // the dots
    for (int i = 0; i < percent / 2; ++i)
      bar += ".";
    // empty space
    for (int i = 0; i < 50 - percent / 2; ++i)
      bar += " ";
    bar += "] ";

    // percentage
    if (percent < 10)
      bar += " ";
    if (percent < 100)
      bar += " ";
    bar += make_string(percent, "%");

    return bar;
  }


  inline std::vector<std::size_t> linear_intersection(const std::vector<std::vector<std::size_t> > &v, const std::vector<std::size_t> &indexes)
  {

  }

}

#endif
