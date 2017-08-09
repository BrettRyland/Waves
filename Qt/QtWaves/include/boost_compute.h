#pragma once

#define BOOST_COMPUTE_USE_OFFLINE_CACHE
//#define BOOST_COMPUTE_HAVE_THREAD_LOCAL
//#define BOOST_COMPUTE_THREAD_SAFE
#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION
// The above #defines need to be defined before including boost.compute, so include this file before any other boost.compute includes.
#include <boost/compute.hpp>