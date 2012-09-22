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

#ifndef MOLDB_TYPETRAITS_H
#define MOLDB_TYPETRAITS_H

namespace MolDB {

  namespace Impl {

    template<typename IfTrueType, typename IfFalseType, int Expr>
    struct SelectHelper { typedef IfTrueType result; };

    template<typename IfTrueType, typename IfFalseType>
    struct SelectHelper<IfTrueType, IfFalseType, false> { typedef IfFalseType result; };

  }

  namespace TypeTraits {

    struct NullType {};

    struct EmptyType {};

    template<int Expr>
    struct StaticAssert;

    template<>
    struct StaticAssert<true> {};


    template<typename T, typename U>
    struct IsSame { enum { result = false }; };

    template<typename T>
    struct IsSame<T, T> { enum { result = true }; };


    template<typename T>
    struct IsPointer { enum { result = false }; };

    template<typename T>
    struct IsPointer<T*> { enum { result = true }; };

    template<typename T>
    struct RemovePointer { typedef T result; };

    template<typename T>
    struct RemovePointer<T*> { typedef T result; };

    /**
     *
     */
    template<typename T>
    struct IsInteger { enum { result = false }; };

    /// @cond IMPL
    template<> struct IsInteger<char> { enum { result = true }; };
    template<> struct IsInteger<signed char> { enum { result = true }; };
    template<> struct IsInteger<unsigned char> { enum { result = true }; };
    template<> struct IsInteger<short> { enum { result = true }; };
    //template<> struct IsInteger<signed short> { enum { result = true }; };
    template<> struct IsInteger<unsigned short> { enum { result = true }; };
    template<> struct IsInteger<int> { enum { result = true }; };
    //template<> struct IsInteger<signed int> { enum { result = true }; };
    template<> struct IsInteger<unsigned int> { enum { result = true }; };
    template<> struct IsInteger<long> { enum { result = true }; };
    //template<> struct IsInteger<signed long> { enum { result = true }; };
    template<> struct IsInteger<unsigned long> { enum { result = true }; };
    /// @endcond

    /**
     * Integer to type template. This is used to turn an integer to type and is
     * mainly used to dispatch to helper functions.
     *
     * Example:
     * @code
     * template<typename T>
     * void fooHelper(T v, Int2Type<true>)
     * {
     *   // v is a pointer
     *   v->bar();
     * }
     *
     * template<typename T>
     * void fooHelper(T v, Int2Type<false>)
     * {
     *   // v is a reference
     *   v.bar();
     * }
     *
     * template<typename T>
     * void foo(T v)
     * {
     *   fooHelper(v, Int2Type<IsPointer<T>::result>());
     * }
     * @endcode
     */
    template<int N>
    struct Int2Type { enum { result = N }; };

    template<typename T>
    struct Type2Type { typedef T result; };


    /**
     * Select a type based on an value of an expression. This works like the
     * expr ? true_value : false_value statement in C++ but operates on types
     * instead.
     *
     * Example:
     * @code
     * // ? : expression with values
     * bool foo = true;
     * int bar = foo ? 0 : 42; // bar = 0;
     *
     * Select<false, int, int*>::result v; // v is of type int*
     * @endcode
     */
    template<int Expr, typename IfTrueType, typename IfFalseType>
    struct Select { typedef typename Impl::SelectHelper<IfTrueType, IfFalseType, Expr>::result result; };

  }

}

#endif
