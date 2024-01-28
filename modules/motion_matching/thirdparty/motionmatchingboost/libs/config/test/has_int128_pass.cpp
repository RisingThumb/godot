//  This file was automatically generated on Thu Oct 25 10:14:36 2012
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-4.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id$
//


// Test file for macro BOOST_HAS_INT128
// This file should compile, if it does not then
// BOOST_HAS_INT128 should not be defined.
// See file boost_has_int128.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifdef BOOST_HAS_INT128
#include "boost_has_int128.ipp"
#else
namespace motionmatchingboost_has_int128 = empty_boost;
#endif

int main( int, char *[] )
{
   return motionmatchingboost_has_int128::test();
}

