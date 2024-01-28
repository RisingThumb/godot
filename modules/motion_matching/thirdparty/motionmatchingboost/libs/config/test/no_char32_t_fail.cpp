//  This file was automatically generated on Thu May 29 09:24:04 2008
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-4.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id$
//


// Test file for macro BOOST_NO_CXX11_CHAR32_T
// This file should not compile, if it does then
// BOOST_NO_CXX11_CHAR32_T should not be defined.
// See file boost_no_char32_t.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifdef BOOST_NO_CXX11_CHAR32_T
#include "boost_no_char32_t.ipp"
#else
#error "this file should not compile"
#endif

int main( int, char *[] )
{
   return motionmatchingboost_no_cxx11_char32_t::test();
}

