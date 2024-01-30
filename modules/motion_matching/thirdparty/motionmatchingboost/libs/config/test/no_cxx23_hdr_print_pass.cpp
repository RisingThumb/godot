//  This file was automatically generated on Sat Sep  9 18:46:50 2023
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-21.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id$
//


// Test file for macro BOOST_NO_CXX23_HDR_PRINT
// This file should compile, if it does not then
// BOOST_NO_CXX23_HDR_PRINT should be defined.
// See file boost_no_cxx23_hdr_print.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifndef BOOST_NO_CXX23_HDR_PRINT
#include "boost_no_cxx23_hdr_print.ipp"
#else
namespace motionmatchingboost_no_cxx23_hdr_print = empty_boost;
#endif

int main( int, char *[] )
{
   return motionmatchingboost_no_cxx23_hdr_print::test();
}

