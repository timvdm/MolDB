/**
 * Copyright (c) 2012, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MOLDB_SEMANTICS_H
#define MOLDB_SEMANTICS_H

#include <string>
#include <sstream>

#include "typetraits.h"

namespace MolDB {

  namespace Impl {

    template<typename T>
    inline std::string semanticsToStringHelper(T s)
    {
      std::stringstream ss;
      ss << s;
      return ss.str();
    }

    template<typename T>
    inline std::string semanticsToStringHelper(T *s)
    {
      std::stringstream ss;
      ss << *s;
      return ss.str();
    }

    template<>
    inline std::string semanticsToStringHelper(char s)
    {
      return semanticsToStringHelper(static_cast<int>(s));
    }

    template<>
    inline std::string semanticsToStringHelper(unsigned char s)
    {
      return semanticsToStringHelper(static_cast<int>(s));
    }

    template<>
    inline std::string semanticsToStringHelper(signed char s)
    {
      return semanticsToStringHelper(static_cast<int>(s));
    }

  }

  template<typename T>
  inline std::string semanticsToString(T s)
  {
    return Impl::semanticsToStringHelper(s);
  }

  namespace Impl {

    template<typename T>
    inline bool semanticsMatchHelper(T v, T w, TypeTraits::Int2Type<false>)
    {
      // T is not a pointer
      return v == w;
    }

    template<typename T>
    inline bool semanticsMatchHelper(T v, T w, TypeTraits::Int2Type<true>)
    {
      // T is a pointer
      // dereference pointers
      return *v == *w;
    }

  }

  template<typename T>
  inline bool semanticsMatch(T v, T w)
  {
    // dispatch to right helper (is T a pointer?)
    return Impl::semanticsMatchHelper(v, w, TypeTraits::Int2Type<TypeTraits::IsPointer<T>::result>());
  }


}

#endif
