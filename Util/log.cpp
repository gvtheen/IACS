/**
 * @file log.cpp
 * @author Matthew Amidon
 *
 * Implementation of the Log class.
 *
 * mlpack is free software; you may redistribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#include "log.hpp"


using namespace util;
//std::ostream& Log::cout ;

std::ofstream Log::fcout("GACatalysis.log",std::ios::out);

PrefixedOutStream Log::Warn(std::cout,"[Warn]: ");

PrefixedOutStream Log::Info(std::cout,"[Info]: ");

#ifdef DEBUG
  PrefixedOutStream Debug(std::cout,"[Debug]: ");;
#else
  NullOutStream Debug;
#endif NullOutStream Log::Debug;

PrefixedOutStream Log::Fatal(std::cout,"[Fatal]: ");

PrefixedOutStream Log::OutputToFile(Log::fcout,"");
// Only do anything for Assert() if in debugging mode.


#ifdef DEBUG
void Log::Assert(bool condition, const std::string& message)
{
  if (!condition)
  {
#ifdef HAS_BFD_DL
    Backtrace bt;

    Log::Debug << bt.ToString();
#endif
    Log::Debug << message << std::endl;

    throw std::runtime_error("Log::Assert() failed: " + message);
  }
}
#else
void Log::Assert(bool /* condition */, const std::string& /* message */)
{ }


#endif
