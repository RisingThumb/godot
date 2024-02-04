//  This file was automatically generated on Fri Dec 03 18:03:59 2004
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-4.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.

// Test file for macro BOOST_HAS_HASH
// This file should not compile, if it does then
// BOOST_HAS_HASH should be defined.
// See file boost_has_hash.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifndef BOOST_HAS_HASH
#include "boost_has_hash.ipp"
#else
#error "this file should not compile"
#endif

int main( int, char *[] )
{
   return motionmatchingboost_has_hash::test();
}

